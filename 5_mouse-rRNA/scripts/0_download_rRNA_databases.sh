
cd cd 2026_Item12_rRNA/
ln -s /Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/rRNA_databases   rRNA_databases
# #!/bin/bash
# # Script to download all required SortMeRNA rRNA databases
# # Run this on a machine with internet access, then transfer the files to your server

# set -e

# # Create directory for rRNA databases
# mkdir -p rRNA_databases
# cd rRNA_databases

# echo "Downloading SortMeRNA rRNA databases (v4.3.6)..."

# # Download all 8 required database files using v4.3.6 URLs
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/rfam-5.8s-database-id98.fasta
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/rfam-5s-database-id98.fasta  
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/silva-arc-16s-id95.fasta
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/silva-arc-23s-id98.fasta
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/silva-bac-16s-id90.fasta
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/silva-bac-23s-id98.fasta
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/silva-euk-18s-id95.fasta
# wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.6/data/rRNA_databases/silva-euk-28s-id98.fasta

# # Create the manifest file with local paths
# cat > sortmerna_database_manifest.txt << EOF
# \$(pwd)/rfam-5.8s-database-id98.fasta
# \$(pwd)/rfam-5s-database-id98.fasta
# \$(pwd)/silva-arc-16s-id95.fasta
# \$(pwd)/silva-arc-23s-id98.fasta
# \$(pwd)/silva-bac-16s-id90.fasta
# \$(pwd)/silva-bac-23s-id98.fasta
# \$(pwd)/silva-euk-18s-id95.fasta
# \$(pwd)/silva-euk-28s-id98.fasta
# EOF

# echo "Download completed!"
# echo "Transfer this 'rRNA_databases' directory to your server at:"
# echo "/Work_bio/gao/projects/2026_Item10_rRNA_removal/scripts/rRNA_databases/"

# # Verify downloads
# echo "Verifying file sizes..."
# ls -lh *.fasta