# Transcription-dependent heterochromatin at the Xist promoter shapes the random choice of the inactive X chromosome

Eleni Kanata, Ingrid Pelaez-Conde, Gemma Noviello, Ilona Dunkel, Lena Milanowska, Till Schwämmle, Melissa Bothe, Rutger A.F. Gjaltema, Edda G. Schulz

https://www.biorxiv.org/content/10.64898/2026.07.22.737064v1

## Abstract
In female mammals, Xist, the master regulator of X-chromosome inactivation (XCI), is expressed monoallelically. This pattern is established during early embryonic development, when the active Xist allele is chosen at random in each cell. How this choice is made remains incompletely understood. Combining knock-down and overexpression strategies in differentiating mouse embryonic stem cells, which recapitulate the onset of random XCI, we identify a role for the repressive chromatin mark H3K9me3 in the XCI initiation. We show that H3K9me3 accumulates at the promoter-proximal region of the silent Xist allele in female cells as monoallelic expression is established. Unexpectedly, this accumulation requires prior transcription of Xist itself, likely during the initial phase of upregulation, when Xist is frequently transcribed in male cells and from both X chromosomes in females. A repressive function of Xist-dependent H3K9me3 accumulation is supported by our finding that premature, transient Xist overexpression primes an allele for future silencing and skews the choice of the inactive X. Xist-dependent H3K9me3 recruitment does not require its antisense transcript Tsix, which can nonetheless enhance subsequent maintenance of the mark. In addition, the X-linked Xist activator RNF12 counteracts H3K9me3 formation, independently of its known target REX1. Our results thus point to facultative heterochromatin formation as a key contributor to choice at the onset of XCI, where activating and repressing mechanisms are intertwined to establish monoallelic Xist expression.

## Description
This repository contains the code used to produce main and supplementary figures of the publication. The data are either provided in this repository, or they are uploaded in GEO under the accession number GSE339369.

## Data analysis
### Whole-Genome Bisulfite Sequencing (WGBS)
Analysis is performed exactly as described here: https://github.com/EddaSchulz/Antisense_paper/tree/main/Align_methylation

### Amplicon-Bisulfite Sequencing
1. Download (or create your own) Amplicons.bed file, that contains coordinates of used amplicons. Major Satellite Repeat coordinates were retrieved from RepeatMasker.

2. Install necessary packages:

FastQC (v0.12.1) quality control tool for high throughput sequencing data from "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"   
MultiQC (v1.26) quality control report tool from "https://github.com/multiqc/multiqc"   
trim_galore (v0.6.4) perl wrapper around FastQC and Cutadapt from "https://github.com/FelixKrueger/TrimGalore"   
samtools (v1.22) collection of C scripts from "http://www.htslib.org/"   
Picard tools (v2.18.25) collection of JAVA scripts from "https://broadinstitute.github.io/picard/"   

bsmapz alignment tool for bisulfite converted DNA from "https://github.com/zyndagj/bsmapz"   
MOABS for methylation calling from "https://github.com/sunnyisgalaxy/moabs"   

R (v4.4.2) from "https://www.r-project.org/"

R-packages:   
epialleleR (v1.14.0)

3. Run the scripts in the following order:
`create_ref_genome.sh`: Given a bed file of amplicon coordinates, it creates a fasta file to use for mapping.   
`amplicon_bisseq.sh`: Pre-processing, mapping of reads and methylation calling. Also outputs qualilty control table and name-sorted bam files to be used for epiallele analysis.   
`amplicon_bis_seq_analysis_v2.R`: extract methylation tables for plotting.   

4. To plot the epiallele distributions use the script in `Fig1_heterochromatin_dynamics`:
`epiallele_densities.R`

### CUT&Tag
Before running the pipeline, you will need the following:
1. Append the AmpR (ampicillin resistance gene) sequence to the mm10 reference genome to be able to map spike-in reads.
```
>AMPR
GGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATC
```
2. Create bowtie2 index of that genome.
3. Download the ENCODE Blacklisted regions (mm10-v2) (Amemiya et al, 2019).
4. Install necessary packages:

FastQC (v0.12.1) quality control tool for high throughput sequencing data from "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"   
MultiQC (v1.26) quality control report tool from "https://github.com/multiqc/multiqc"   
trim_galore (v0.6.4) perl wrapper around FastQC and Cutadapt from "https://github.com/FelixKrueger/TrimGalore"   
bedtools (v2.30.0) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"   
samtools (v1.22) collection of C scripts from "http://www.htslib.org/"   
Bowtie2 (v2.5.0) aligner from "https://github.com/BenLangmead/bowtie2"   
Deeptools2 (v3.5.6) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"   
WiggleTools (v1.2.11) from "https://github.com/Ensembl/Wiggletools"   
USCS tools from "https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads"   
epic2 peak caller from "https://github.com/biocore-ntnu/epic2"   

R (v4.4.2) from "https://www.r-project.org/"

R-packages:   
DiffBind (v3.16.0)



5. Run the scripts in the following order:
`CUTnTag_pipeline.sh` : Quality control, trimming and alignment of reads. Merging of replicates and creation of bigwig tracks.   
`average_norm_reads_over_RE57.sh` : Counts how many normalized reads cover RE57.   
`epic2.sh` : Peak calling.   
`diffbind.R` : Differential peak analysis with epic2 peaks. Example script for RNF12-KO and DKO vs WT.   
`mpileup_SNPs.sh` : Read count over SNPs.

### RT-qPCR

1. Run qPCR in Quant-Studio 7 Flex real-time PCR machine (Thermo Fisher Scientific). 

2. Save results as `.eds` file, analyze with Design & Analysis 2.6.0 software (Thermo Fisher Scientific) and export result table as a single `.xlsx` file. 

3. Run `qPCR_analysis.R` script to calculate relative gene expression of target genes (normalized to the geometric mean of housekeeping Rplp0 (also called Arpo) and Rrm2).

## Figures
### Figure 1
`WGBS_at_RE57.R`: Figure 1d   
`dna_meth_timecourse.R`: Figure 1e   

`cnt_reads_at_RE57.R`: Supplementary Figure 1a   
`RNA-seq.R`: Supplementary Figure 1c   
`wgbs_vs_amplicon.R`: Supplementary Figure 1d   
`epiallele_densities.R`: Supplementary Figure 1e   

### Figure 2
`flowjo_workspace.R`: Figure 2b, Supplementary Figure 2g   
`fish_counts.R`: Figure 2c, Supplementary Figure 2f   
`sorted_total_reads_at_RE57.R`: Figure 2e   
`allelic_reads.R`: Figure 2f   

`karyotyping.R`: Supplementary Figure 2d   
`qPCR.R`: Supplementary Figure 2e   

### Figure 3
`qPCR_data.R`: Figure 3b, Supplementary Figure 3b,c   
`reads_at_RE57.R`: Figure 3d   
`dna_meth.R`: Figure 3e   
`pyrosequencing.R`: Figure 3h, Supplementary Figure 3d   
`SNP-FISH.R`: Figure 3k   

`karyotyping.R`: Supplementary Figure 3a

### Figure 4
`qPCR_data.R`: Figure 4d,h, Supplementary Figure 4c   
`reads_at_RE57.R`: Figure 4f,j, Supplementary Figure 4e   
`dna_meth.R`: Figure 4k, Supplementary Figure 4f   

`allelic_reads.R`: Supplementary Figure 4a

### Figure 5
`qPCR_data.R`: Figure 5c, Supplementary Figure 5c,d,f,g,h   
`reads_at_RE57.R`: Figure 5e,j, Supplementary Figure 6e   
`dna_meth.R`: Figure 5f,k, Supplementary Figure 6c,f   
`fish.R`: Figure 5l, Supplementary Figure 5e   
`volcano_plots.R`: Figure 5n, Supplementary Figure 7d   

`spike-in_norm.R`: Supplementary Figure 6a-b, Supplementary Figure 7a

`qPCR_data_DKO.R`: Supplementary Figure 7b,c   
`heatmap_RNF12KO_vs_DKO.R`: Supplementary Figure 7e   
`new_peaks.R`: Supplementary Figure 7f   
