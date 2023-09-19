#!/bin/bash

# Description: This script checks if a given VCF sample file header contains references to 'hg38' or related terms.
#              It uses bcftools to perform the check and outputs an error message if the references are absent.
# useage: hg38Checker.sh <sample_file.vcf.bgz>

# Load necessary modules
module load hurcs bcftools

# Filepath to the sample VCF file
sampleFile=$1

# Use bcftools to extract and inspect the headers of the sample VCF file
# Check if 'hg38', 'GRCh38', or 'HG38' appear in any of the headers
header_check=$(bcftools view -h "${sampleFile}" | egrep -c 'hg38|GRCh38|HG38')

# If 'hg38' or related terms are not found in the header, exit with an error message
if [[ $header_check -eq 0 ]]; then
    echo "'hg38', 'GRCh38', or 'HG38' does not appear in the header of the sample file."
    exit 1
fi
