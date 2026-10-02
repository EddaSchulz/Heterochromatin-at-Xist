#!/bin/bash
# Get allelic info from SNPs
# Nelly Kanata
# OWL Schulz
# Created: 13.04.2026
# Modified: 24.04.2026

date
DIR='/your/working/directory/' # replace this with your working directory
bam_dir=${DIR}bam/
SNP_dir=${DIR}/SNP_reads/
# SNP_bed='TX1072_SNPs.bed' # file with SNPs between B6/Cast mouse strains

mkdir ${SNP_dir}/
cd ${SNP_dir}


# filter SNPs around RE57:
# $SNP_bed was used, here I provide the resulting RE57_SNPs.bed file

## RE57 proximal chrX:103,479,041-103,482,895
# bedtools intersect -a $SNP_bed -b <(echo -e "chrX\t103479041\t103482895") > RE57_SNPs.bed


# Get nucleotide info in each SNP position
    
for sample in $(basename -a $bam_dir/*.mapped.sort.blacklisted.bam | sed -E 's/_Rep.*//' | sort -u) ; do
    echo -e "Reading $sample"

    # For samples with 3 replicates:
    echo -e 'Chr\tStart\tRef\tReads_R1\tNucleotides_R1\tQualities_R1\tReads_R2\tNucleotides_R2\tQualities_R2\tReads_R3\tNucleotides_R3\tQualities_R3' > ${sample}_SNPs.txt

    samtools mpileup -l ${SNP_dir}RE57_SNPs.bed \
    ${bam_dir}/${sample}_Rep*_S*.mapped.sort.blacklisted.bam >> ${sample}_SNPs.txt

done



