#!/bin/bash
## Nelly Kanata
## OWL Schulz
## Created: 27.11.2025
## Modified: 14.01.2026

DIR='/your/working/directory/' # replace this with your working directory
bam_dir=${DIR}/bam/
fastq_dir=${DIR}fastq/
epic2_dir=${DIR}/epic2_peaks/

mkdir -p $epic2_dir

cd $DIR

for f in $(ls ${fastq_dir}*R1.fastq.gz | sed 's/_R1.*//' | uniq| xargs -n1 basename); do

    treatment_file=${bam_dir}/${f}.mapped.sort.blacklisted.bam
    output_file="${epic2_dir}/${f}_epic2.bed"

    # Skip if the epic2 output already exists
    if [[ -f "$output_file" ]]; then
        echo "Skipping $f — output already exists: $output_file"
        continue
    fi

    echo "Running epic2 for $f"
    epic2 --treatment $treatment_file \
        --guess-bampe \
        --genome mm10 \
        --keep-duplicates \
        --gaps-allowed 1 \
        --output ${output_file}

done

      