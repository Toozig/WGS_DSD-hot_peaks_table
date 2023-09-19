#!/bin/bash


# Script: getSamples.sh
# Description: turns a sample VCF file into a TSV, adds some layer data 
# Usage: bash getSamples.sh <chrom> <sample.vcf.bgz> <peaks.bed>

# Usage function
print_usage() {
    echo "Usage: $(basename "$0") pathFile regionFile"
    echo "    pathFile: Path to the input VCF/BCF file"
    echo "    regionFile: Path to the region BED file"
    echo "    output: output file name"
}

# Check for correct number of arguments
if [ $# -ne 3 ]; then
    echo "Error: Incorrect number of arguments."
    print_usage
    exit 1
fi

# Get command line arguments
pathFile=$1
regionFile=$2
output=$3

# Extract file names for generating output file name
sample_name=$(basename "$pathFile" .txt)
region_name=$(basename "$regionFile" .bed)


# Define headers and formats for different data sources
gnomAD_header='data/read_only/layers_data/headers/gnomAD.header'
gnomAD_format='CHROM,POS,REF,ALT,FILTER,AF,AF_popmax'

region_header='data/read_only/layers_data/headers/interval_ID.header'
region_format='CHROM,FROM,TO,INTERVAL_ID'

RM_bed='data/read_only/layers_data/repeatsMasker.bed'
RM_header='data/read_only/layers_data/headers/repeatsMasker.header'
RM_format='CHROM,FROM,TO,RM_ID,RM_LENGTH,RM_STRAND'

GH_bed='data/read_only/layers_data/GeneHancer_AnnotSV_elements_v5.15.txt'
GH_header='data/read_only/layers_data/headers/geneHancer_AnnotSV_elements.header'
GH_format='CHROM,FROM,TO,GHid,GH_is_elite,GH_type'

# FORMAT="%CHROM\t%POS\t%REF\t%ALT\t%INTERVAL_ID\t%RM_ID\t%GHid\t%GH_is_elite\t%GH_type[\t%GT\t%DP\t%GQ]\n"  with repeatMasker
FORMAT="%CHROM\t%POS\t%REF\t%ALT\t%INTERVAL_ID\t%GHid\t%GH_is_elite\t%GH_type[\t%GT\t%DP\t%GQ]\n"

# Function to check exit code and exit if not 0
check_exit_code() {
  if [ $? -ne 0 ]; then
    echo "Error: $1"
    exit 1
  fi
}


echo $regionFile $pathFile
# Perform various operations on the input data
bcftools view -R $regionFile $pathFile  |   # get only variants in the region files
  bcftools annotate -x FILTER,INFO/AF |   # remove the FILTER and AF field from the vcf
  bcftools annotate -a $regionFile -h $region_header -c $region_format  |   # adds the peaks ID to each variants
  bcftools annotate -a $GH_bed -h $GH_header -c $GH_format | # adds the GeneHancer data
#  bcftools annotate -a $RM_bed -h $RM_header -c $RM_format |  # adds the repeatmasker data
 bcftools query -H -f "$FORMAT"  -o $output   # saves as table

# Check the exit code of the last command
check_exit_code "bcftools merge samples and annotations"

echo "Processing completed."
