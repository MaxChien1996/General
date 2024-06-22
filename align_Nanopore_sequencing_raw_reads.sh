#!/bin/bash




# Make the script executable -> "chmod +x align_Nanopore_sequencing_raw_reads.sh"
# Run the script with appropriate arguments -> "./align_Nanopore_sequencing_raw_reads.sh reference.fasta query.fastq.gz output.bam"





reference_file="$1"
query_file="$2"
output_file="$3"

minimap2 -a "$reference_file" "$query_file" > "${output_file%.bam}".sam # use "%"" to remove suffix ".bam"

samtools view -bS "${output_file%.bam}".sam > "${output_file%.bam}".bam
rm "${output_file%.bam}".sam

samtools sort "${output_file%.bam}".bam > "${output_file%.bam}".sorted.bam
rm "${output_file%.bam}".bam

samtools index "${output_file%.bam}".sorted.bam