#!/bin/bash
# Count normalized reads at RE57
# Nelly Kanata
# OWL Schulz
# Created: 11.08.2023
# Modified: 22.10.2024

date
DIR='/your/working/directory/' # replace this with your working directory

ucsc_tools=~/Tools/ucsc_tools/
fastq_dir=${DIR}fastq/
bw_dir=${DIR}bigwig/
comp_dir=${DIR}/comparison/

mkdir $comp_dir

cd $DIR
for f in $(ls ${fastq_dir}*R1.fastq.gz | sed 's/_R1.*//' | uniq| xargs -n1 basename); do

    echo $f
    
    ${ucsc_tools}/bigWigAverageOverBed $bw_dir/$f'.bw' \
    RE57.bed \
    $comp_dir/$f.tab
    sed -i '1s/^/name\tsize\tcovered\tsum\tmean0\tmean \n/' $comp_dir/$f.tab

done

echo "Done!"