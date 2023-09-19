#!/bin/bash

# Description: This script processes a VCF sample file by removing 'INFO/AF' annotations, normalizing variants,
#              and indexing the processed file. The resulting processed VCF file is stored in an output directory,
#              and the path to the processed file is added to a sample list file.
#
# Usage: bash process_vcf.sh <sampleFile> <outputDir> <sampleListFile>
#    <sampleFile>: Path to the input VCF sample file
#    <outputDir>: Directory to store the processed VCF file
#    <sampleListFile>: File to store the list of processed sample paths

# Input Variables
sampleFile=$1           # Path to the input VCF sample file
outputDir=$2            # Directory to store the processed VCF file
sampleListFile=$3       # File to store the list of processed sample paths

# Load necessary modules
module load hurcs bcftools

# Remove 'INFO/AF' annotations, normalize variants, and save the processed VCF as a compressed file
bcftools annotate -x "INFO/AF" "$sampleFile" | \
    bcftools norm -m +any | \
    bcftools norm -m -any -Oz -o "${outputDir}/processed.${sampleFile}"

# Index the processed VCF using tabix
tabix -p vcf "${outputDir}/processed.${sampleFile}"

# Append the path of the processed file to the sample list file
echo "${outputDir}/processed.${sampleFile}" >> "$sampleListFile"
