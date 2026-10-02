#!/bin/bash
## Nelly Kanata
## OWL Schulz
## Created: 15.04.2024
## Modified: 15.04.2024

# create reference genome for amplicon bisulfite seq analysis
# script from MethPanel tool "https://github.com/thinhong/MethPanel"
# USAGE: ./create_ref_genome.sh

# Define variables
DIR='./' # replace this with your working directory
inAmp=$DIR/ref/Amplicons.bed"
build_name="mm10"
panel_name="RE57-58"
add=200

# new directories
REF=${DIR}/ref/ # where to place the index of amplicons sequences

mkdir -p $REF 

cd $DIR
  
### 1. prepare files
# name the ref fasta
build_name_panel="${build_name}_${panel_name}"
refs="$REF/${build_name_panel}.fa"
# name the amplicon coordinate with .extend.merged
inAmp_extend_merged_bed="$REF/${panel_name}.extend.merged.bed"
# name the amplicon coordinate with _sorted_mapped
inAmp_sorted_mapped="$REF/${panel_name}_sorted_mapped.bed"

# 2. Extend the regions by +-200 and check amplicon overlapping given a bed file of amplicon starts and ends.  
# If there is an overlapping, these overlapped amplicons are merged into one amplicon
bedtools merge -o collapse -delim "." -c 4 -i <( awk -v add="$add" '{OFS="\t"}NR>1{print $1, $2-add, $3+add, $5}' "$inAmp" \
| sort -k1,1 -k2,2n ) | \
sort -k1,1 -k2,2n > "$inAmp_extend_merged_bed"



### 2. Create a map of mm10 coordinate vs amplicon coordinate and a map of merged amplicons
bedtools intersect -a "$inAmp_extend_merged_bed" -b <(awk '{OFS="\t"}NR>1{print $1, $2, $3, $5}' "$inAmp" \
| sort -k1,1 -k2,2n) -wa -wb | \
awk '{OFS="\t"}{print $0, $6-$2, $7-$2, $7-$6}' | sort -k1,1 -k2,2n > "$inAmp_sorted_mapped"

# example
# cat "${inAmp/_sorted.bed/_sorted_mapped.bed}" | column -t
# chr1   23882688   23883128   P1_ID3           chr1   23882888   23882928   P1_ID3           200  240  40
# chr1   24469072   24469525   P2_IL22RA1       chr1   24469272   24469325   P2_IL22RA1       200  253  53

### 3. Get sequences
build_name_panel="${build_name}_${panel_name}"
refs="$REF/${build_name_panel}.fa"

mkdir -p "$REF/tmp"
rm -f $refs # Remove any existing reference genome fasta file.

for i in `awk '{OFS=";"}!/^#/{print $1,$2,$3,$4}' "$inAmp_sorted_mapped"`; do 
	chr=$(echo $i| cut -d";" -f1);
	s0=$(bc <<< "$(echo $i| cut -d";" -f2) - $add"); # start
	s1=$(bc <<< "$(echo $i| cut -d";" -f2) - $add + 1"); 
	e=$(bc <<< "$(echo $i| cut -d";" -f3) + $add"); # end
	n=$(echo $i| cut -d";" -f4); # name
	fn="${n}::${chr}:${s0}-${e}"; # filename
	l="${chr}:${s1},${e}"; # location
	wget -O "$REF/tmp/${fn}.fasta" "http://genome.ucsc.edu/cgi-bin/das/${build_name}/dna?segment=$l";
    #echo $n $fn
	sleep 3; # Pauses the script execution for 3 seconds to avoid overloading the server.
done

# append all downloaded sequences to the ref fasta file
for i in `find "$REF/tmp" -name "*.fasta"`; do 
	echo $i; 
	echo -e ">$(basename ${i/.fasta/})" >> "$refs"; 
	awk '!/^</' $i| tr -d "\n" >> "$refs";
	echo -e "\n" >> "$refs" 
done
rm -rf "$REF/tmp"

