mkdir -p /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery/qc
mkdir -p /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery/results

cp /home/gao/projects_2026H2/3_cho_wgs_species_confirm/qc/wt1_fastp.html /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery/qc/

cp /home/gao/projects_2026H2/3_cho_wgs_species_confirm/qc/multiqc/multiqc_report.html /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery/qc/


cp /home/gao/projects_2026H2/3_cho_wgs_species_confirm/results/cho_flagstat.txt /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery/results/

cp /home/gao/projects_2026H2/3_cho_wgs_species_confirm/results/dhfr_gff.txt /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery/results/

tree /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery/

mv /home/gao/projects_2026H2/3_cho_wgs_species_confirm/delivery /home/gao/projects_2026H2/3_cho_wgs_species_confirm/Research_Report