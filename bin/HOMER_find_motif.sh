#!/bin/bash

module load hurcs homer
# Set the path to the HOMER program
HOMER_PATH="" # Replace with the actual path to your HOMER binaries

# Check if HOMER is installed and accessible
#if ! command -v "$HOMER_PATH/findMotifsGenome.pl" &>/dev/null; then
#  echo "HOMER not found or the path is incorrect. Please set the HOMER_PATH variable correctly."
#  exit 1
#fi

# Input BED file (hg38 format)
INPUT_BED="$1"  # Replace with the path to your hg38 BED file


bed_name=`basename $INPUT_BED .bed `

# Output directory for HOMER results
OUTPUT_DIR="/cs/prt3/DSD_project/HOMER/$bed_name"


# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

WORKDIR="/cs/prt3/DSD_project/tmp/HOMER/"
 
# Step 0, doing stupid change, dunno why it is needed but...
PEAK_FILE="$WORKDIR/$bed_name.homer"
cp $INPUT_BED $PEAK_FILE

# Step 1: Run HOMER to find motifs in the BED file


echo  "$INPUT_BED" hg38 "$OUTPUT_DIR" -size 200 -mask -preparsedDir "$WORKDIR" 

findMotifsGenome.pl "$PEAK_FILE" hg38 "$OUTPUT_DIR" -size 200 -mask -preparsedDir "$WORKDIR"  -preparse



# Step 2: Extract motif occurrences from HOMER output
MOTIF_FILE="$OUTPUT_DIR/homerResults.html"
PEAKS_WITH_MOTIF="$OUTPUT_DIR/peaks_with_motif.txt"

# Check if the HOMER output file exists
if [ ! -f "$MOTIF_FILE" ]; then

  echo "HOMER did not generate the expected output. Please check the input and try again."
  exit 1
fi
# Parse the HOMER output to extract motif occurrences and corresponding peaks
 
findMotifsGenome.pl "$PEAK_FILE" hg38 "$OUTPUT_DIR" -size 200 -mask -preparsedDir "$WORKDIR" -find "$OUTPUT_DIR/homerMotifs.all.motifs" > $OUTPUT_DIR/peaks_motif.tsv

echo "Motif discovery and peak analysis completed successfully."