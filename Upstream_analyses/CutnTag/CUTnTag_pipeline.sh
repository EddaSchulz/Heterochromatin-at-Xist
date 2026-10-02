#!/bin/bash

# CUT&Tag-seq pipeline
# after: https://yezhengstat.github.io/CUTTag_tutorial/ (Zheng Y et al (2020). Protocol.io)

# Run as:
# bash ./CUTnTag_pipeline.sh > std.2.out 2> std.2.err & disown

date
DIR='/your/working/directory/' # replace this with your working directory

#reference genome with AmpR sequence (spike-in)
mm10=/path/to/your/genome_index/bowtie2/mm10_AmpR/mm10_AmpR # replace this with the path to your bowtie2 genome index

# blacklisted regions from (Amemiya et al, 2019)
mm10_bl=files/mm10-blacklist.v2.bed

fastq_dir=${DIR}fastq/
trim_dir=${DIR}trimmed/
bam_dir=${DIR}bam/
bw_dir=${DIR}bigwig/
qc_dir=${DIR}/qc/

# Create directories
mkdir -p $fastq_dir
mkdir -p $trim_dir
mkdir -p $bam_dir
mkdir -p $bw_dir
mkdir -p $bed_dir
mkdir -p $qc_dir


# pipeline for PE reads
cd $DIR

# for every sample, run the pipeline
for f in $(ls ${fastq_dir}*R1.fastq.gz | sed 's/_R1.*//' | uniq| xargs -n1 basename); do
  
  echo $f


  # Generate fastqc for all files
  fastqc -o ${fastq_dir} -f fastq ${fastq_dir}${f}_R1.fastq.gz
  fastqc -o ${fastq_dir} -f fastq ${fastq_dir}${f}_R2.fastq.gz

  # Trimming adapter sequences
  echo -e "Trimming $f with trim_galore\n"
  (trim_galore -j 8 --paired -o trimmed ${fastq_dir}${f}_R1.fastq.gz ${fastq_dir}${f}_R2.fastq.gz)  2> ${trim_dir}${f}_trimmingStats.txt

  # fastqc for trimmed files
  fastqc -o ${trim_dir} --threads 2 -f fastq ${trim_dir}${f}_R*_val_*.fq.gz

  # Alignment to reference genome
  echo -e "Alignment for $f with bowtie2\n"

  bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 2000 -p 16 -x $mm10 \
  -1 ${trim_dir}${f}_R1_val_1.fq.gz -2 ${trim_dir}${f}_R2_val_2.fq.gz \
  -S ${bam_dir}${f}.sam &> ${bam_dir}${f}_mappingStats.txt


  # Filter and keep only the mapped read pairs  
  samtools view -bS -F 0x04 ${bam_dir}${f}.sam > ${bam_dir}${f}.mapped.bam
  samtools sort -@ 4 ${bam_dir}${f}.mapped.bam -o ${bam_dir}${f}.mapped.sort.bam
  
  # Removing blacklisted region with a new blacklist file, important for Tsix repeat region
  echo -e "Remove blacklisted regions for $f"
  bedtools intersect -v -a ${bam_dir}${f}.mapped.sort.bam -b ${mm10_bl} > ${bam_dir}${f}.mapped.sort.blacklisted.bam
  samtools index ${bam_dir}${f}.mapped.sort.blacklisted.bam

  # Create bigwig tracks
  echo "creating RPKM normalized coverage track for $f"
  bamCoverage -b ${bam_dir}${f}.mapped.sort.blacklisted.bam -o ${bw_dir}${f}.bw --extendReads -bs 1 --normalizeUsing RPKM -p 20
  
done



##### QC metrics ####

# multiqc -d $fastq_dir -o $qc_dir
multiqc $fastq_dir*.zip -d $fastq_dir -o $qc_dir
multiqc $trim_dir*.zip -d $trim_dir -o $qc_dir'/trimmed' 

# Generating qc_file for the data
echo -e 'Sample\tTotal\tMapped\tMultimappers\tFiltered' > ${qc_dir}QC_metrics.txt 

for f in $(ls ${fastq_dir}*H3K9*R1.fastq.gz | sed 's/_R1.*//' | uniq| xargs -n1 basename); do

  echo -e 'Returning qc info for' ${f}
  fragments=$(grep -Po -m 1 'Total reads processed.*' ${trim_dir}${f}_trimmingStats.txt | grep -Po '[0-9,]*' | tr -d ,)
  map=$(grep -Po '[0-9, /.]*% overall alignment rate' ${bam_dir}${f}_mappingStats.txt| grep -Po '[0-9,/.]*')
  mapFrac=$(printf '%.4f' $(echo $map / 100 | bc -l))
  multimappers=$(grep -Po '[0-9, /.]*%\) aligned concordantly >1 times' ${bam_dir}${f}_mappingStats.txt| grep -Po '[0-9,/.]*' | head -n 1)
  filtered=$(expr $(samtools view -c ${bam_dir}${f}.mapped.sort.blacklisted.bam) / 2)
  echo -e "$f\t$fragments\t$mapFrac\t$multimappers\t$filtered" >> ${qc_dir}QC_metrics.txt 
done

#### Counting spike-in reads ####

files=${DIR}spikein_files/
mkdir $files

cd ${bam_dir}

# Count the number of AmpR reads
# all files with the same antibody (ab)

ab="H3K9me3"
touch ${files}$ab\_norm_reads.txt
echo -e "Sample\tTotal_reads\tAmpR_reads\tAmpR_normalized" >> ${files}$ab\_norm_reads.txt

for f in $(ls *$ab*.bam | cut -d'.' -f-1 | uniq) 
do
    echo -e "Counting Reads in $f"
    total=$(expr $(samtools view -c -@ 20 $f\.mapped.sort.blacklisted.bam) / 2)
    AmpR=$(expr $(samtools view -c -@ 20  $f\.mapped.sort.blacklisted.bam AMPR) / 2 )

    AmpR_norm=$(echo "$AmpR * 1000 / $total " | bc -l) #RPKM

    echo -e "$f\t$total\t$AmpR\t$AmpR_norm" >> ${files}$ab\_norm_reads.txt
done




##### Merging replicates for visualization ####
merged_dir=${bw_dir}/merged_replicates/
chrom_sizes=mm10_AMPR.chrom.sizes
ucsc_tools=~/Tools/ucsc_tools/

mkdir ${merged_dir}/

# Calculate mean of Replicates
for sample in $(basename -a $fastq_dir/*.fastq.gz | sed -E 's/_Rep.*//' | sort -u) ; do
    echo -e "Merging $sample"
    wiggletools mean \
        ${bw_dir}/${sample}*.bw \
        > ${merged_dir}/${sample}_mean.wig 
    
    ${ucsc_tools}/wigToBigWig ${merged_dir}/${sample}_mean.wig ${chrom_sizes} ${merged_dir}/${sample}_mean.bw
done


# removing unnecessary files
  rm -f ${fastq_dir}*fastqc.zip
  rm -rf ${trim_dir}*fq.gz
  rm -rf ${trim_dir}*.fastq.gz_trimming_report.txt
  rm -rf ${bam_dir}*.sam
  rm -rf ${bam_dir}*.mapped.bam
  rm -rf ${bam_dir}*.mapped.sort.bam
