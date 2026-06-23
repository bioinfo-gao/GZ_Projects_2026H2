ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:SORTMERNA (4)'

Caused by:
  Can't stage file https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-arc-23s-id98.fasta -- reason: Network is unreachable


Command executed:

  sortmerna \
      --ref rfam-5.8s-database-id98.fasta --ref rfam-5s-database-id98.fasta --ref silva-arc-16s-id95.fasta --ref silva-arc-23s-id98.fasta --ref silva-bac-16s-id90.fasta --ref silva-bac-23s-id98.fasta --ref silva-euk-18s-id95.fasta --ref silva-euk-28s-id98.fasta \
       \
      --reads C3_C_trimmed_1_val_1.fq.gz --reads C3_C_trimmed_2_val_2.fq.gz \
      --threads 7 \
      --workdir . \
      --aligned rRNA_reads --fastx --other non_rRNA_reads \
      --paired_in \
      --out2 \
      --num_alignments 1 -v --index 0
  
  
          mv non_rRNA_reads_fwd.f*q.gz C3_C_1.non_rRNA.fastq.gz
          mv non_rRNA_reads_rev.f*q.gz C3_C_2.non_rRNA.fastq.gz
          mv rRNA_reads.log C3_C.sortmerna.log
  
  
  cat <<-END_VERSIONS > versions.yml
  "NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:SORTMERNA":
      sortmerna: $(echo $(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*$//')
  END_VERSIONS

Command exit status:
  -

Command output:
  (empty)

Container:
  /home/gao/.singularity/nf-core/depot.galaxyproject.org-singularity-sortmerna-4.3.6--h9ee0642_0.img

Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

 -- Check '.nextflow.log' file for details
ERROR ~ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting

 -- Check '.nextflow.log' file for details