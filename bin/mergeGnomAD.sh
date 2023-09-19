#!/bin/bash

output=$1

python3 << EOF
import pandas as pd
import sys

output_name=sys.argv[1]
csv_files = sys.argv[2:]
dfs = []  # List to store individual DataFrames from each CSV file

# Read and append each CSV file to the list
for file in csv_files:
    df = pd.read_csv(file, sep='\t')
    df.columns = [col.split(']')[1] for col in df.columns]
    df = df[df.CHROM == chrom]
    df = df.set_index([col for col in df.columns if ':' not in col])
    dfs.append(df)

# Concatenate DataFrames into a single DataFrame
concatenated_df = pd.concat(dfs, axis=1).sort_index()
index_col = concatenated_df.columns[:4]
concatenated_df = concatenated_df.sort_values(index_col.tolist())
concatenated_df.to_csv(output_name. sep='\t')
EOF


# Compress the output file using bgzip
bgzip ${output}

# Create an index for the compressed output file using tabix
tabix -s1 -b2 -e2 ${output}.gz