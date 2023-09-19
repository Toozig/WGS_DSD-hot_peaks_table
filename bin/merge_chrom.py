#!/usr/bin/env python
import sys
import pandas as pd

# Define constant column names for the DataFrame
INDEX_COLS = ['CHROM', 'POS', 'REF', 'ALT']
GNOMAD_COLS = INDEX_COLS + ['FILTER', 'AF', 'AF_popmax']

# Define indices for command line arguments
GNOMAD_PATH_IDX = 2
OUTPUT_NAME_IDX = 3
INPUT_FILES_IDX = 4

def add_gnomAD(df, gnomAD_path):
    """
    Add gnomAD data to the given DataFrame based on the common INDEX_COLS.
    
    Args:
        df (pd.DataFrame): The input DataFrame.
        gnomAD_path (str): Path to the gnomAD data file.
        
    Returns:
        pd.DataFrame: The DataFrame with added gnomAD data.
    """
    gnomad = pd.read_csv(gnomAD_path, compression='gzip', sep='\t', header=None)
    gnomad.columns = GNOMAD_COLS

    gnomad = gnomad.set_index(INDEX_COLS)
    df = df.reset_index().set_index(INDEX_COLS)
    gnomad = gnomad[gnomad.index.isin(df.index)]
    df = pd.concat([gnomad, df], axis=1)
    return df

def main(csv_files, chrom):
    """
    Concatenate data from multiple CSV files based on the given chromosome.
    
    Args:
        csv_files (list): List of input CSV files.
        chrom (str): Chromosome to filter data on.
        
    Returns:
        pd.DataFrame: Concatenated DataFrame with filtered data.
    """
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
    return concatenated_df

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Error: Insufficient command line arguments.")
        print("Usage: python script_name.py chrom gnomAD_file.csv output_name input_file1.tsv input_file2.tsv ...")
        sys.exit(1)
    
    chrom = sys.argv[1]
    gnomAD_path = sys.argv[GNOMAD_PATH_IDX]
    output_name = sys.argv[OUTPUT_NAME_IDX]
    input_files = sys.argv[INPUT_FILES_IDX:]
    
    df = main(input_files, chrom)
    df = add_gnomAD(df, gnomAD_path)
    df.to_csv(f"{output_name}.tsv", sep='\t')
    print('Done concatenation and gnomAD addition.')
