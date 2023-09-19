#!/bin/bash

# Script: getGnomad.sh
# Description: Retrive gnomAD FILTER, AF, AF_popmax for specifics locations from gnomAD DB
# Usage: bash getGnomad.sh <chrom> <peaks.bed> <output_file>

# This script extracts specific information from a gnomAD VCF file using bcftools.
# It requires two command line arguments: chromosome and region file.

# Check if the required command line arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <chromosome> <regionFile> <outputPath>"
    exit 1
fi

# Assign input arguments to variables
chrom=$1
regionFile=$2
output=$3

# URL of the gnomAD VCF file for the specified chromosome
#VCF="https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${chrom}.vcf.bgz"
VCF="/cs/prt3/gnomAD_data/chr_vcf/gnomad.genomes.v3.1.2.sites.${chrom}.vcf.bgz"
ln -s "${VCF}.tbi" .
format="%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%AF\t%AF_popmax\n"  # Define desired output format


# Use bcftools query to extract information from the specified region of the VCF file
bcftools query -f "${format}" -R ${regionFile} -o ${output}  ${VCF} 


