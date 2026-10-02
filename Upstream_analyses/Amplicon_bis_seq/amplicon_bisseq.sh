#!/bin/bash
## Nelly Kanata
## OWL Schulz
## Created: 10.04.2024
## Modified: 28.02.2026

# amplicon bisulfite sequencing reads pre-processing and mapping
# USAGE: ./amplicon_biseq.sh

# Define variables
DIR='./' # replace this with your working directory
bsmapz=/your/tool/directory/bsmapz #change this for the directory where you installed bsmapz
PICARD=/your/tool/directory/picard.jar #change this to the directory where you installed PICARD.
REF=$DIR/ref/mm10_RE57-58.fa

# new directories
fastq_dir=${DIR}fastq/
trim_dir=${DIR}trimmed/
bam_dir=${DIR}bam/
qc_dir=${DIR}/qc/
fastqc_dir=$qc_dir/fastqc/
mcall_dir=${DIR}/mcall/

# Create directories
mkdir -p $trim_dir
mkdir -p $bam_dir/stats
mkdir -p $bam_dir/dedup/
mkdir -p $fastqc_dir/trimmed
mkdir -p $fastqc_dir/untrimmed
mkdir -p $mcall_dir

cd $DIR

# 1. Run FastQC on all fastq files in indicated directory
fastqc $fastq_dir/*.gz -t 2 -o=$fastqc_dir'/untrimmed/'
multiqc -d $fastqc_dir/untrimmed -o $qc_dir

for f in $(ls ${fastq_dir}*R1.fastq.gz | sed 's/_R1.*//' | uniq| xargs -n1 basename); do
  echo $f
  # Extract the last two characters as the sample number
  sample_number="${f##*S}"
    
    
    # 2. Trimming adapter sequences
    echo -e "Trimming $f with trim_galore\n"
   
    trim_galore --fastqc_args "--outdir ${fastqc_dir}/trimmed/" -o $trim_dir --clip_R1 10 --three_prime_clip_R1 5 --clip_R2 15 \
    --three_prime_clip_R2 5 --paired --cores 4 --basename $f ${fastq_dir}${f}_R1.fastq.gz ${fastq_dir}${f}_R2.fastq.gz 2> ${trim_dir}${f}_trimmingStats.txt



    # 3. Alignment to reference genome 
    echo -e "Alignment for $f with BSMapz\n"


    $bsmapz -a "${trim_dir}${f}_val_1.fq.gz" -b "${trim_dir}${f}_val_2.fq.gz" \
    -d $REF -m 0 -g 3 -n 1 -v 0.1 -w 100 -o $bam_dir/${f}.bam &> ${bam_dir}'/stats/'${f}_mappingStats.txt
    
    # options notes:
    # -n 1: map SE or PE reads to all 4 strands
    # -m 0: minimal insert size allowed 0
    # -g 3: gap size allowed
    # -v 0.1: mismatch rate
    # -w 100: maximum number of equal best hits to count 100

    samtools sort -T $bam_dir/temp --threads 4 -o "$bam_dir/${f}.sorted.bam" "$bam_dir/${f}.bam"
		

    # index
    samtools index "$bam_dir/${f}.sorted.bam"
    

    # Duplication rate for stats (optional)
    java -jar $PICARD MarkDuplicates I=${bam_dir}${f}.sorted.bam O=${bam_dir}/dedup/${f}.sorted.dupMarked.bam METRICS_FILE=${bam_dir}/dedup/${f}_dedupMetric.txt
  

done

# check if samples need to be merged (multiple sequencing runs)
for f in $(ls ${fastq_dir}*R1.fastq.gz | sed 's/_R1.*//' | uniq | xargs -n1 basename); do
  # Define the sample names
  sample="${f}"
  sample_b="${f}b"
  
  # Define the BAM file paths
  bam_file="${bam_dir}/${sample}.sorted.bam"
  bam_file_b="${bam_dir}/${sample_b}.sorted.bam"
  merged_bam="${bam_dir}/${sample}_merged.sorted.bam"

  # Check if both BAM files exist
  if [[ -f "$bam_file" && -f "$bam_file_b" ]]; then
    echo "Merging ${bam_file} and ${bam_file_b} into ${merged_bam}"
    samtools merge -f "$merged_bam" "$bam_file" "$bam_file_b"
  else
    echo "One or both of the files ${bam_file} and ${bam_file_b} do not exist."
  fi
done

# methylation calling
for f in $(ls ${fastq_dir}*R1.fastq.gz | sed 's/_R1.*//' | uniq | xargs -n1 basename); do
    
    echo "Running Methylation Calling CpG"
    cd $mcall_dir

  # Determine which BAM file to use
  if [[ -f "${bam_dir}/${f}_merged.sorted.bam" ]]; then
    echo "Using merged BAM file."
    bam_to_use="${bam_dir}/${f}_merged.sorted.bam"
  elif [[ -f "${bam_dir}/${f}.sorted.bam" ]]; then
    echo "Using standard BAM file"
    bam_to_use="${bam_dir}/${f}.sorted.bam"
  else
    echo "Error: Neither BAM file exists for ${f}"
    continue  # Skip to the next sample if neither BAM file exists
  fi

    mcall \
    --threads 4 \
    --reference $REF \
    --sampleName $f \
    --mappedFiles $bam_to_use \
    --outputDir $mcall_dir \
    --webOutputDir $mcall_dir \
    --trimWGBSEndRepairPE2Seq 0 \
    --trimWGBSEndRepairPE1Seq 0 \
    --reportCHX 0

done


# Quality control
multiqc -d $fastqc_dir/trimmed -o $qc_dir


output_file=${DIR}/qc/"readcount_table.txt"
# Initialize the table with headers
echo -e "Sample\tTotal\tMapped\tMultimappers\tduplicates\tRE57_reads\tRE58_reads\ttotal_Xist_reads\ttotal_mapped_reads_primary\tbisulfite_conversion_rate" > $output_file

for f in $(ls ${fastq_dir}*R1.fastq.gz | sed 's/_R1.*//' | uniq| xargs -n1 basename); do
  echo $f

  fragments=$(grep -Po -m 1 'Total reads processed.*' ${trim_dir}${f}_trimmingStats.txt | grep -Po '[0-9,]*' | tr -d ,)

  map=$(grep -Po 'aligned pairs: \d+ \(\K[0-9.]+(?=%)' ${bam_dir}/stats/${f}_mappingStats.txt)

  multimappers=$(grep -Po 'non-unique pairs: \d+ \(\K[0-9.]+(?=%)' ${bam_dir}/stats/${f}_mappingStats.txt)

  dedup_file="${bam_dir}/dedup/${f}_dedupMetric.txt"
  if [[ -f $dedup_file ]]; then
    dedup=$(grep -Po 'Unknown Library.*' $dedup_file | awk '{print $8}')
  else
    dedup="NA"  # Set a default value or handle it as needed
  fi

  RE57_reads=$(samtools view -c -F 256 "$bam_dir/${f}.sorted.bam" "FJ14-FJ15.FJ10-FJ11::chrX:103481327-103482781" )
  # -F 256 keeps out secondary alignments

  RE58_reads=$(samtools view -c -F 256 "$bam_dir/${f}.sorted.bam" "FJ28-FJ29.FJ26-FJ27::chrX:103482701-103483891" )

  total_Xist_reads=$((RE57_reads + RE58_reads))

  total_mapped_reads_primary=$(samtools view -c -F 256 "$bam_dir/${f}.sorted.bam")

  bisulfite_conversion_rate=$(grep 'bisulfite conversion ratio' ${bam_dir}/${f}".sorted.bam_stat.txt" | head -n 1 | awk -F ' = ' '{print $2}')

  # Append the values to the output file
  echo -e "${f}\t${fragments}\t${map}\t${multimappers}\t${dedup}\t${RE57_reads}\t${RE58_reads}\t${total_Xist_reads}\t${total_mapped_reads_primary}\t${bisulfite_conversion_rate}" >> $output_file

done


# Preparation for epiallele analysis (optional)
# epialleleR (package for epiallele analysis) needs name-sorted bam files

for f in $(ls ${fastq_dir}/*R1.fastq.gz | sed 's/_R1.*//' | uniq| xargs -n1 basename); do
  
    echo $f

    # sort merged bam file on name
    samtools sort -n -o "${bam_dir}/${f}.name-sorted.bam" "${bam_dir}/${f}.sorted.bam"

done

echo "Done!"

