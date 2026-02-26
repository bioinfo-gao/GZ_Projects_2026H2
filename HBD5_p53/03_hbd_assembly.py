# 03_hbd_assembly.py
import pysam
from collections import Counter, defaultdict

def get_consensus_with_quality(sequences):
    """改进的共识组装算法：返回共识序列及每位置的碱基频率信息"""
    if not sequences: return "", []

    consensus_seq = ""
    base_frequencies = []  # 记录每个位置的碱基频率

    for i in range(len(sequences[0])):
        column = [s[i] for s in sequences if i < len(s)]
        base_counts = Counter(column)

        # 获取最常见碱基及其频率
        most_common_base, count = base_counts.most_common(1)[0]
        freq_info = {base: count/base_counts.total() for base, count in base_counts.items()}

        consensus_seq += most_common_base
        base_frequencies.append(freq_info)

    return consensus_seq, base_frequencies

def calculate_family_confidence(base_frequencies):
    """
    计算家族一致性（Confidence）分数
    - 基于每个位置上主要碱基的频率
    - 返回0-1之间的置信度分数，1表示完全一致，0表示完全不一致
    """
    if not base_frequencies:
        return 0.0
    
    # 计算每个位置的主要碱基频率的平均值
    major_base_freqs = []
    for freq_info in base_frequencies:
        if freq_info:
            max_freq = max(freq_info.values())
            major_base_freqs.append(max_freq)
    
    if not major_base_freqs:
        return 0.0
    
    average_major_freq = sum(major_base_freqs) / len(major_base_freqs)
    return average_major_freq

def classify_confidence_level(confidence_score):
    """
    根据置信度分数分类家族质量
    """
    if confidence_score >= 0.95:
        return "HIGH"
    elif confidence_score >= 0.85:
        return "MEDIUM"
    elif confidence_score >= 0.70:
        return "LOW"
    else:
        return "VERY_LOW"

def get_duplex_consensus(strand1_consensus, strand1_freqs, strand2_consensus, strand2_freqs):
    """
    Duplex consensus: 只有当两条互补链都支持同一个突变时，才认定为真实变异
    strand2_consensus 应该是 strand1_consensus 的反向互补序列
    """
    if not strand1_consensus or not strand2_consensus:
        return None, []
    
    # 首先验证strand2是否确实是strand1的反向互补
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    expected_strand2 = ''.join([complement[base] for base in strand1_consensus])[::-1]
    
    # 如果长度不匹配，返回None
    if len(strand1_consensus) != len(strand2_consensus):
        return None, []
    
    duplex_consensus = ""
    duplex_frequencies = []
    
    # 对于每个位置，检查两条链是否一致支持相同的碱基
    for i in range(len(strand1_consensus)):
        base1 = strand1_consensus[i]
        freq1 = strand1_freqs[i].get(base1, 0) if i < len(strand1_freqs) else 0
        
        # strand2的对应位置（从末尾开始）
        j = len(strand2_consensus) - 1 - i
        if j >= 0 and j < len(strand2_consensus):
            base2 = strand2_consensus[j]
            freq2 = strand2_freqs[j].get(base2, 0) if j < len(strand2_freqs) else 0
            
            # 检查base2是否是base1的互补碱基
            expected_complement = complement.get(base1, 'N')
            if base2 == expected_complement and freq1 >= 0.7 and freq2 >= 0.7:
                # 两条链都高置信度支持，使用strand1的碱基
                duplex_consensus += base1
                # 合并频率信息
                combined_freq = {base1: min(freq1, freq2)}  # 使用较低的置信度
                duplex_frequencies.append(combined_freq)
            else:
                # 不一致或置信度不够，标记为不确定
                duplex_consensus += "N"
                duplex_frequencies.append({'N': 1.0})
        else:
            duplex_consensus += "N"
            duplex_frequencies.append({'N': 1.0})
    
    return duplex_consensus, duplex_frequencies

def classify_error_source(mutation_freq, family_size):
    """
    根据突变频率和家族大小分类错误来源：
    - 低频自发突变：在多个家族中出现，频率较低但一致
    - PCR第一次扩增错误：在单个家族中高频出现（>50%），因为早期错误会被大量复制
    - PCR后续扩增错误：在单个家族中低频出现（<30%），因为晚期错误只影响少数拷贝
    """
    if mutation_freq >= 0.5:
        return "PCR_first_amplification_error"
    elif mutation_freq >= 0.1:
        return "PCR_later_amplification_error"
    else:
        return "low_frequency_spontaneous_mutation"

def analyze_mutations_detailed(consensus_seq, base_frequencies, reference_seq=None):
    """
    详细分析突变类型和可能来源，区分不同类型的PCR错误
    """
    mutations = []
    for i, (cons_base, freq_info) in enumerate(zip(consensus_seq, base_frequencies)):
        # 检查是否有其他碱基存在（非共识碱基）
        other_bases = {base: freq for base, freq in freq_info.items() if base != cons_base and base != 'N'}

        for base, freq in other_bases.items():
            # 根据频率和上下文判断突变类型
            mut_type = classify_error_source(freq, len(base_frequencies))
            
            mutations.append({
                'position': i,
                'original_base': cons_base,
                'variant_base': base,
                'frequency': freq,
                'type': mut_type
            })

    return mutations

def deduplicate_molecules_duplex(bam_path):
    """
    Duplex HBD深度去重：利用"二级样本Index + UMI + Read2起始坐标 + strand"作为分子唯一标识
    在duplex模式下，我们需要分别处理来自Watson链(alpha UMI)和Crick链(beta UMI)的reads
    """
    # 存储唯一的分子: { (sample_index, umi, read2_start, strand) : [seq1, seq2, ...] }
    unique_molecules = defaultdict(list)
    # 同时存储每个分子的详细信息
    molecule_details = defaultdict(list)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped: continue

            # 从查询名称中提取完整信息
            query_parts = r.query_name.split("_")
            if len(query_parts) >= 6:
                # 格式: READ_x_watson_alpha_UMI_primaryIdx_sampleIdx_umi_used
                # 或者: READ_x_crick_beta_UMI_primaryIdx_sampleIdx_umi_used
                strand = None
                umi_type = None
                umi_used = None
                primary_index = None
                sample_index = None
                
                # 查找strand和umi_type信息
                for i, part in enumerate(query_parts):
                    if part in ['watson', 'crick']:
                        strand = part
                        if i + 1 < len(query_parts) and query_parts[i + 1] in ['alpha', 'beta']:
                            umi_type = query_parts[i + 1]
                            if i + 2 < len(query_parts):
                                umi_used = query_parts[i + 2]
                        break
                
                # 如果没有找到完整信息，尝试从最后提取UMI
                if strand is None or umi_used is None:
                    umi_used = query_parts[-1]
                    strand = "unknown"
                    umi_type = "unknown"
                
                # 尝试提取sample_index和primary_index
                if len(query_parts) >= 8:
                    try:
                        primary_index = query_parts[-3]  # 倒数第三个
                        sample_index = query_parts[-2]   # 倒数第二个  
                    except:
                        sample_index = "unknown"
                        primary_index = "unknown"
                else:
                    sample_index = "unknown"
                    primary_index = "unknown"
            else:
                # 备用方案：如果格式不匹配
                umi_used = query_parts[-1] if query_parts else "unknown"
                strand = "unknown"
                sample_index = "unknown"
                primary_index = "unknown"

            # 获取Read2的起始比对位置（随机断裂点）
            if r.is_read1:
                read2_start = r.next_reference_start if r.has_tag('MC') else r.reference_start
            else:
                read2_start = r.reference_start

            # 四位一体键：二级样本Index + UMI + Read2起始坐标 + strand
            key = (sample_index, umi_used, read2_start, strand)

            unique_molecules[key].append(r.query_sequence)
            molecule_details[key].append({
                'query_name': r.query_name,
                'query_sequence': r.query_sequence,
                'reference_start': r.reference_start,
                'reference_end': r.reference_end,
                'mapping_quality': r.mapping_quality,
                'cigar': r.cigarstring,
                'primary_index': primary_index,
                'sample_index': sample_index,
                'umi': umi_used,
                'strand': strand,
                'umi_type': umi_type,  # 明确记录是alpha还是beta UMI
                'read2_start': read2_start
            })

    return unique_molecules, molecule_details

def group_duplex_families(unique_molecules, molecule_details):
    """
    将来自同一条原始双链DNA的Watson链(alpha)和Crick链(beta)家族组合成duplex家族
    """
    duplex_families = {}
    
    # 创建一个映射：(sample_index, read2_start) -> {watson_data, crick_data}
    family_groups = defaultdict(lambda: {'watson': None, 'crick': None})
    
    for key, sequences in unique_molecules.items():
        sample_index, umi, read2_start, strand = key
        
        # 使用(sample_index, read2_start)作为组键，忽略UMI差异
        group_key = (sample_index, read2_start)
        
        consensus_seq, base_freqs = get_consensus_with_quality(sequences)
        family_info = {
            'sequences': sequences,
            'consensus': consensus_seq,
            'base_frequencies': base_freqs,
            'umi': umi,
            'details': molecule_details[key]
        }
        
        if strand == 'watson':
            family_groups[group_key]['watson'] = family_info  # alpha UMI 家族
        elif strand == 'crick':
            family_groups[group_key]['crick'] = family_info    # beta UMI 家族
    
    # 创建duplex家族
    for group_key, strands in family_groups.items():
        if strands['watson'] is not None and strands['crick'] is not None:
            # 两条链都有数据，可以进行duplex consensus
            duplex_families[group_key] = strands
        elif strands['watson'] is not None:
            # 只有Watson链 (alpha UMI)
            duplex_families[group_key] = {'watson': strands['watson'], 'crick': None}
        elif strands['crick'] is not None:
            # 只有Crick链 (beta UMI)
            duplex_families[group_key] = {'watson': None, 'crick': strands['crick']}
    
    return duplex_families

def count_total_reads_from_fastq(fastq_file):
    """从FASTQ文件计算总reads数"""
    count = 0
    try:
        with open(fastq_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    count += 1
    except FileNotFoundError:
        count = 0
    return count

def generate_read_statistics_report(unique_molecules, raw_r1_file="raw_R1.fastq", clean_r1_file="clean_R1.fq"):
    """生成详细的reads统计报告"""
    
    # 1. 原始reads总数（来自图片模拟的数据）
    original_reads = count_total_reads_from_fastq(raw_r1_file)
    
    # 2. 清洗后的reads数（成功提取UMI和索引的reads）
    cleaned_reads = count_total_reads_from_fastq(clean_r1_file)
    
    # 3. 映射到参考基因组的reads数
    mapped_reads = sum(len(seqs) for seqs in unique_molecules.values())
    
    # 4. 分子家族数量（unique molecules after deduplication）
    unique_molecules_count = len(unique_molecules)
    
    # 5. PCR重复率计算
    pcr_duplicates = mapped_reads - unique_molecules_count
    
    # 6. 最终真实reads数目（共识序列数）
    final_consensus_reads = unique_molecules_count
    
    # 创建统计报告
    stats_report = {
        'original_reads': original_reads,
        'cleaned_reads': cleaned_reads,
        'mapped_reads': mapped_reads,
        'unique_molecules': unique_molecules_count,
        'pcr_duplicates': pcr_duplicates,
        'final_consensus_reads': final_consensus_reads,
        'pcr_duplication_rate': (pcr_duplicates / mapped_reads * 100) if mapped_reads > 0 else 0,
        'cleaning_efficiency': (cleaned_reads / original_reads * 100) if original_reads > 0 else 0,
        'mapping_rate': (mapped_reads / cleaned_reads * 100) if cleaned_reads > 0 else 0
    }
    
    return stats_report

def write_statistics_report(stats_report, output_file="read_analysis_report.txt"):
    """写入详细的统计报告"""
    with open(output_file, "w") as f:
        f.write("HBD5_p53 Duplex Read Analysis Statistics Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("1. Original Reads (from simulation/images):\n")
        f.write(f"   Total original reads: {stats_report['original_reads']:,}\n")
        f.write("   - These represent the raw simulated data including all PCR duplicates\n")
        f.write("   - Generated from the TP53 simulation with introduced errors\n\n")
        
        f.write("2. Cleaned Reads (after UMI/index extraction):\n")
        f.write(f"   Successfully processed reads: {stats_report['cleaned_reads']:,}\n")
        f.write(f"   Cleaning efficiency: {stats_report['cleaning_efficiency']:.2f}%\n")
        f.write("   - Reads that successfully passed UMI and index extraction\n")
        f.write("   - Failed extractions are typically due to missing GGG anchor or short sequences\n\n")
        
        f.write("3. Mapped Reads (aligned to reference genome):\n")
        f.write(f"   Successfully mapped reads: {stats_report['mapped_reads']:,}\n")
        f.write(f"   Mapping rate: {stats_report['mapping_rate']:.2f}%\n")
        f.write("   - Reads that successfully aligned to the hg38 reference genome\n\n")
        
        f.write("4. Molecular Deduplication Analysis:\n")
        f.write(f"   Unique molecular families: {stats_report['unique_molecules']:,}\n")
        f.write(f"   PCR duplicates removed: {stats_report['pcr_duplicates']:,}\n")
        f.write(f"   PCR duplication rate: {stats_report['pcr_duplication_rate']:.2f}%\n")
        f.write("   - Each unique family represents one original DNA molecule\n")
        f.write("   - PCR duplicates are multiple reads from the same original molecule\n\n")
        
        f.write("5. Final Consensus Reads:\n")
        f.write(f"   Final high-fidelity consensus sequences: {stats_report['final_consensus_reads']:,}\n")
        f.write("   - These represent the true biological molecules after error correction\n")
        f.write("   - In duplex mode, only variants supported by both strands are retained\n")
        f.write("   - Family confidence scores indicate internal consistency of each molecular family\n\n")
        
        f.write("Summary:\n")
        f.write(f"   From {stats_report['original_reads']:,} original reads,\n")
        f.write(f"   we obtained {stats_report['final_consensus_reads']:,} high-confidence consensus sequences.\n")
        f.write(f"   This represents a {stats_report['final_consensus_reads']/stats_report['original_reads']*100:.2f}% yield\n")
        f.write("   of true biological molecules after removing PCR duplicates and sequencing errors.\n")

def assemble_families_duplex(bam_path):
    # 获取去重后的分子家族（duplex模式）
    unique_molecules, molecule_details = deduplicate_molecules_duplex(bam_path)
    
    # 生成详细的统计报告
    stats_report = generate_read_statistics_report(unique_molecules)
    write_statistics_report(stats_report)
    
    print(f"{'SampleIdx':<12} {'UMI':<8} {'UMIType':<8} {'Read2Start':<12} {'Strand':<8} {'Reads':<6} {'Confidence':<10} {'Consensus_Seq_Snippet'}")
    print("-" * 130)

    # 组合duplex家族
    duplex_families = group_duplex_families(unique_molecules, molecule_details)
    
    # 用于统计所有家族的突变信息
    all_mutations = []
    duplex_validated_mutations = []

    with open("hbd_final_assembly_report.csv", "w") as o:
        # 添加更多列来描述详细的错误分类和duplex验证，包括家族一致性
        o.write("Gene,Chrom,Pos,SampleIndex,Read2Start,WatsonReads,CrickReads,WatsonConsensus,CrickConsensus,DuplexConsensus,MutationInfo,VariantType,DuplexValidated,FamilyConfidence,ConfidenceLevel,WildTypeRatio\n")

        for group_key, strands in duplex_families.items():
            sample_index, read2_start = group_key
            
            watson_data = strands.get('watson')
            crick_data = strands.get('crick')
            
            watson_consensus = watson_data['consensus'] if watson_data else ""
            watson_freqs = watson_data['base_frequencies'] if watson_data else []
            watson_reads = len(watson_data['sequences']) if watson_data else 0
            
            crick_consensus = crick_data['consensus'] if crick_data else ""
            crick_freqs = crick_data['base_frequencies'] if crick_data else []
            crick_reads = len(crick_data['sequences']) if crick_data else 0
            
            # 计算家族置信度
            family_confidence = 0.0
            confidence_level = "UNKNOWN"
            
            if watson_data and crick_data:
                # Duplex模式：使用duplex consensus的置信度
                duplex_consensus, duplex_freqs = get_duplex_consensus(
                    watson_consensus, watson_freqs, 
                    crick_consensus, crick_freqs
                )
                if duplex_consensus:
                    family_confidence = calculate_family_confidence(duplex_freqs)
                    confidence_level = classify_confidence_level(family_confidence)
                else:
                    # 如果duplex consensus失败，使用单链中较好的那个
                    conf1 = calculate_family_confidence(watson_freqs)
                    conf2 = calculate_family_confidence(crick_freqs)
                    family_confidence = max(conf1, conf2)
                    confidence_level = classify_confidence_level(family_confidence)
                    duplex_consensus = None
            elif watson_data:
                # 只有Watson链 (alpha UMI)
                family_confidence = calculate_family_confidence(watson_freqs)
                confidence_level = classify_confidence_level(family_confidence)
                duplex_consensus = None
            elif crick_data:
                # 只有Crick链 (beta UMI)
                family_confidence = calculate_family_confidence(crick_freqs)
                confidence_level = classify_confidence_level(family_confidence)
                duplex_consensus = None
            else:
                family_confidence = 0.0
                confidence_level = "NONE"
                duplex_consensus = None
            
            # 分析突变
            mutations = []
            if duplex_consensus:
                mutations = analyze_mutations_detailed(duplex_consensus, duplex_freqs)
                # 标记为duplex validated
                for mut in mutations:
                    mut['duplex_validated'] = True
                    duplex_validated_mutations.append(mut)
            else:
                # 如果没有有效的duplex consensus，使用单链数据
                if watson_consensus:
                    mutations = analyze_mutations_detailed(watson_consensus, watson_freqs)
                    for mut in mutations:
                        mut['duplex_validated'] = False
                elif crick_consensus:
                    mutations = analyze_mutations_detailed(crick_consensus, crick_freqs)
                    for mut in mutations:
                        mut['duplex_validated'] = False
            
            # 获取染色体和位置信息
            chrom = "chr17"  # TP53位于chr17
            pos = read2_start if read2_start is not None else 0
            gene = "TP53" if chrom == "chr17" else "Unknown"

            # 计算野生型比例
            wildtype_ratio = 0
            if duplex_consensus:
                wildtype_ratio = sum([freq_info.get(base, 0) for base, freq_info in zip(duplex_consensus, duplex_freqs)]) / len(duplex_freqs) if duplex_freqs else 0
            elif watson_consensus:
                wildtype_ratio = sum([freq_info.get(base, 0) for base, freq_info in zip(watson_consensus, watson_freqs)]) / len(watson_freqs) if watson_freqs else 0

            # 格式化突变信息
            mut_info = ";".join([f"{m['position']}:{m['original_base']}->{m['variant_base']}({m['frequency']:.2f},{m['type']})" for m in mutations])

            # 打印部分结果到屏幕
            umi_display = watson_data['umi'] if watson_data else (crick_data['umi'] if crick_data else 'N/A')
            umi_type_display = 'alpha' if watson_data else ('beta' if crick_data else 'N/A')
            strand_display = 'both' if watson_data and crick_data else ('watson' if watson_data else ('crick' if crick_data else 'unknown'))
            total_reads = watson_reads + crick_reads
            consensus_display = duplex_consensus or watson_consensus or crick_consensus
            
            print(f"{sample_index:<12} {umi_display:<8} {umi_type_display:<8} {read2_start:<12} {strand_display:<8} {total_reads:<6} {family_confidence:.3f} ({confidence_level:<8}) {consensus_display[:30]}...")

            # 写入完整报表
            o.write(f"{gene},{chrom},{pos},{sample_index},{read2_start},{watson_reads},{crick_reads},{watson_consensus},{crick_consensus},{duplex_consensus if duplex_consensus else 'N/A'},{mut_info if mut_info else 'None'},family_specific,{len([m for m in mutations if m.get('duplex_validated', False)])>0},{family_confidence:.3f},{confidence_level},{wildtype_ratio:.3f}\n")

            # 记录此家族的所有突变信息
            for mut in mutations:
                all_mutations.append({
                    'chrom': chrom,
                    'pos': pos + mut['position'],  # 实际基因组位置
                    'sample_index': sample_index,
                    'mutation': mut,
                    'duplex_validated': mut.get('duplex_validated', False),
                    'family_confidence': family_confidence,
                    'confidence_level': confidence_level
                })

    # 统计分子家族大小分布
    family_sizes = []
    family_confidences = []
    for strands in duplex_families.values():
        size = 0
        if strands.get('watson'):
            size += len(strands['watson']['sequences'])
        if strands.get('crick'):
            size += len(strands['crick']['sequences'])
        family_sizes.append(size)
        
        # 计算置信度
        if strands.get('watson') and strands.get('crick'):
            conf1 = calculate_family_confidence(strands['watson']['base_frequencies'])
            conf2 = calculate_family_confidence(strands['crick']['base_frequencies'])
            family_confidences.append(max(conf1, conf2))
        elif strands.get('watson'):
            family_confidences.append(calculate_family_confidence(strands['watson']['base_frequencies']))
        elif strands.get('crick'):
            family_confidences.append(calculate_family_confidence(strands['crick']['base_frequencies']))
        else:
            family_confidences.append(0.0)
    
    with open("molecule_family_distribution.txt", "w") as f:
        f.write("FamilySize\tCount\n")
        size_counts = {}
        for size in family_sizes:
            size_counts[size] = size_counts.get(size, 0) + 1
        for size, count in sorted(size_counts.items()):
            f.write(f"{size}\t{count}\n")
    
    # 添加家族置信度分布
    with open("family_confidence_distribution.txt", "w") as f:
        f.write("ConfidenceLevel\tCount\tAverageConfidence\n")
        confidence_levels = {'HIGH': [], 'MEDIUM': [], 'LOW': [], 'VERY_LOW': [], 'NONE': []}
        for conf in family_confidences:
            level = classify_confidence_level(conf)
            confidence_levels[level].append(conf)
        
        for level, conf_list in confidence_levels.items():
            if conf_list:
                avg_conf = sum(conf_list) / len(conf_list)
                f.write(f"{level}\t{len(conf_list)}\t{avg_conf:.3f}\n")
            else:
                f.write(f"{level}\t0\t0.000\n")

    # 生成详细的突变报告
    with open("mutation_analysis_report.txt", "w") as f:
        f.write("Chrom\tPos\tSampleIndex\tPositionInSeq\tOriginalBase\tVariantBase\tFrequency\tType\tDuplexValidated\tFamilyConfidence\tConfidenceLevel\tFamilySize\n")
        for mut_record in all_mutations:
            mut = mut_record['mutation']
            # 计算家族大小
            family_size = 0
            for strands in duplex_families.values():
                if strands.get('strand1'):
                    family_size += len(strands['strand1']['sequences'])
                if strands.get('strand2'):
                    family_size += len(strands['strand2']['sequences'])
            f.write(f"{mut_record['chrom']}\t{mut_record['pos']}\t{mut_record['sample_index']}\t{mut['position']}\t{mut['original_base']}\t{mut['variant_base']}\t{mut['frequency']:.3f}\t{mut['type']}\t{mut_record['duplex_validated']}\t{mut_record['family_confidence']:.3f}\t{mut_record['confidence_level']}\t{family_size}\n")

    # 统计总体突变情况
    total_mutations = len(all_mutations)
    duplex_validated_count = len([m for m in all_mutations if m['duplex_validated']])
    non_duplex_count = total_mutations - duplex_validated_count
    
    # 统计置信度分布
    high_conf_count = len([m for m in all_mutations if m['confidence_level'] == 'HIGH'])
    medium_conf_count = len([m for m in all_mutations if m['confidence_level'] == 'MEDIUM'])
    low_conf_count = len([m for m in all_mutations if m['confidence_level'] in ['LOW', 'VERY_LOW']])

    print(f"\n突变分析总结:")
    print(f"总突变数: {total_mutations}")
    print(f"Duplex验证突变数: {duplex_validated_count}")
    print(f"单链突变数: {non_duplex_count}")
    print(f"Duplex验证率: {duplex_validated_count/total_mutations*100:.2f}%" if total_mutations > 0 else "Duplex验证率: 0%")
    print(f"\n家族置信度分布:")
    print(f"HIGH置信度家族: {high_conf_count}")
    print(f"MEDIUM置信度家族: {medium_conf_count}")
    print(f"LOW/VERY_LOW置信度家族: {low_conf_count}")
    
    # 打印统计摘要
    print(f"\nRead Statistics Summary:")
    print(f"原始reads总数: {stats_report['original_reads']:,}")
    print(f"清洗后reads数: {stats_report['cleaned_reads']:,}")
    print(f"映射reads数: {stats_report['mapped_reads']:,}")
    print(f"唯一分子家族数: {stats_report['unique_molecules']:,}")
    print(f"最终共识reads数: {stats_report['final_consensus_reads']:,}")
    print(f"PCR重复率: {stats_report['pcr_duplication_rate']:.2f}%")

assemble_families_duplex("sorted.bam")