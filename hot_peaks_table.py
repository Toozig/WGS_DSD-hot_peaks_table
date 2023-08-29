"""
Script Description:
-------------------
This script performs analysis on variant data and pedigree information,
 calculating various statistics for intervals and samples. 
 It then saves the analysis results to an Excel file. 
 The script takes input files for variant data and pedigree information, and optionally an upload path for the output Excel file.

Usage:
------
python hot_peaks_table.py.py VAR_CSV_PATH SAMPLE_METADATA_PATH OUTPUT_PATH [UPLOAD_PATH]

Arguments:
----------
VAR_CSV_PATH (str): Path to the CSV file containing variant data.
SAMPLE_METADATA_PATH (str): Path to the Excel file containing sample metadata (pedigree information).
OUTPUT_PATH (str): Path to the directory where the output Excel file will be saved.
UPLOAD_PATH (str, optional): Path for uploading the output Excel file. Can be omitted.

Example Usage:
--------------
1. Basic usage without specifying upload path:
python hot_peaks_table.py.py data/pipeline_outputs/variants_with_layers/qualityDSD_variants.csv data/read_only/samples/sample_metadata.xlsx \
     /data/pipeline_output/hot_peaks/version1
2. Usage with specifying an upload path:
python hot_peaks_table.py.py data/pipeline_outputs/variants_with_layers/qualityDSD_variants.csv data/read_only/samples/sample_metadata.xlsx \
     /data/pipeline_output/hot_peaks/version1 Nitzan_Gonen_lab/Joint_projects/WGS_on_DSD/Ido/hot_peaks

Note:
-----
- The script requires the Pandarallel library for parallel processing.
- Ensure that the specified paths are correct and accessible.
- If UPLOAD_PATH is omitted, the output Excel file won't be uploaded.
"""

import numpy as np
import pandas as pd
import os
import sys
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=8)
import  bin.Gonen_func as gf

VAR_CSV_PATH = 1
SAMPLE_MATADTA_PATH = 2
OUTPUT_PATH = 3
UPLOAD_PATH = 4
GENES_LOCATIONS_FILE = "hg38_dsd_genes_locations.bed"

def get_interval_stats(df_in, pedg_df):
    """
    Calculate various statistics for intervals based on the provided dataframes.

    Parameters:
        df_in (DataFrame): Input dataframe containing variant information.
        pedg_df (DataFrame): Pedigree dataframe containing sample information.

    Returns:
        pd.Series: Calculated statistics for the intervals.
    """
    df = df_in.drop(columns='INTERVAL_ID')
    
    # Calculate the sum of variants for each sample across all intervals
    sum_per_sample = df.sum(axis=0)
    
    # Determine if any variant exists for each sample across all intervals
    exist_per_sample = df.any(axis=0)
    
    # Create an empty Series to store the calculated statistics
    sum_series = pd.Series(dtype=float)
    
    # Extract proband IDs from the pedigree dataframe
    probands =  pedg_df[pedg_df.fam_relation == 0].ID.astype(str)
    
    # Calculate statistics for all sources
    sum_series['total n probands'] = exist_per_sample[exist_per_sample.index.isin(probands)].sum()
    sum_series['total n non-DSD'] = exist_per_sample.sum() - sum_series['total n probands']
    sum_series['total n proband variants'] = sum_per_sample[sum_per_sample.index.isin(probands)].sum()
    sum_series['total n non-DSD variants'] = sum_per_sample.sum() - sum_series['total n proband variants']
    
    # Calculate statistics for each data source
    for label in pedg_df.source.unique():
        label_id = pedg_df[pedg_df.source == label].ID.astype(str)
        probans_id = label_id[label_id.isin(probands)]
        
        # Statistics for probands and non-DSD samples for each label
        sum_series[f'{label} n probands'] = exist_per_sample[exist_per_sample.index.isin(probans_id)].sum()
        sum_series[f'{label} n non-DSD'] = exist_per_sample[exist_per_sample.index.isin(label_id)].sum() - sum_series[f'{label} n probands']
        sum_series[f'{label} n proband variants'] = sum_per_sample[sum_per_sample.index.isin(probans_id)].sum()
        sum_series[f'{label} n non-DSD variants'] = sum_per_sample.sum() - sum_series[f'{label} n proband variants']
    
    return sum_series

def bool_variant_df(df):
    """
    Create a boolean table of samples versus variants.

    Parameters:
        df (DataFrame): Input dataframe containing variant information.

    Returns:
        DataFrame: Boolean table indicating presence of variants in samples.
    """
    # Extract genotype columns by filtering column names ending with 'GT'
    gt = df[[i for i in df.columns if i.endswith('GT')]].replace(' ',np.nan).replace('',np.nan).notna()
    
    # Rename genotype columns for clarity
    gt.columns = [i.replace(':GT','') for i in gt.columns]
    
    # Concatenate the INTERVAL_ID column with the boolean genotype table
    return pd.concat([df[['INTERVAL_ID']], gt], axis=1)

# ... (similar explanations for other functions)

# some rows have a few genes in 'DSDgenes_1.5mb' and this function splits them to several rows
def handle_multi_gene_rows(df):
    # Create a copy of the DataFrame to work with
    new_df = df.copy()
    df.reset_index(inplace=True)
    # Split values in 'DSDgenes_1.5mb' column and create new rows
    new_rows = []
    for index, row in new_df.iterrows():
        genes = row['DSDgenes_1.5mb'].split(',')
        for gene in genes:
            new_row = row.copy()
            new_row['DSDgenes_1.5mb'] = gene
            new_rows.append(new_row)

    # Create a new DataFrame with the updated rows
    new_df = pd.DataFrame(new_rows, columns=new_df.columns)
    # Reset the index of 'df' to keep 'INTERVAL_ID' as a regular column
    new_df = new_df.join(df.set_index(['CHROM', 'from', 'to'])['INTERVAL_ID'], on=['CHROM', 'from', 'to'])
    return new_df

def calculate_minimal_distance_from_gene(df, locations):
    # Merge locations with df based on 'DSDgenes_1.5mb'
    merged_df = df.merge(locations, left_on='DSDgenes_1.5mb', right_on='gene', how='left')

    # Calculate distances
    merged_df['dist_from_start'] = abs(merged_df['from'] - merged_df['start'])
    merged_df['dist_from_end'] = abs(merged_df['to'] - merged_df['start'])
    # Calculate minimum distance from both start and end of gene
    merged_df['distance_from_nearest_DSD_TSS'] = merged_df[['dist_from_start', 'dist_from_end']].min(axis=1)

    # Drop unnecessary columns
    result_df = merged_df.drop(columns=['start', 'end', 'gene', 'dist_from_start', 'dist_from_end', 'chr'])
    return result_df

def add_dsd_distance(df):
    # Read DSD genes' locations file
    locations = pd.read_table(GENES_LOCATIONS_FILE)
    locations.columns = ['chr', 'start', 'end', 'gene']
    # split df by values in 'DSDgenes_1.5mb'
    new_df = handle_multi_gene_rows(df)

    result_df = calculate_minimal_distance_from_gene(new_df, locations)

    # Keep only rows with minimal distance from DSD gene
    result_df = result_df.sort_values(by=['INTERVAL_ID', 'distance_from_nearest_DSD_TSS']).drop_duplicates(subset='INTERVAL_ID')
    result_df.set_index('INTERVAL_ID', inplace=True)
    
    return result_df

def get_info_table(df):
    relevant_coulmns = ['CHROM','from', 'to', 'length', 'DSDgenes_1.5mb','geneHancer', 'GHid', 'GH_is_elite', 'GH_type']
    df = df[~df.INTERVAL_ID.duplicated()].set_index('INTERVAL_ID')
    return df[relevant_coulmns]

def main(sample_file_path, pedg_path, output_file, upload_path=None):
    """
    Main function to perform the analysis and save results.

    Parameters:
        sample_file_path (str): Path to the sample input file.
        pedg_path (str): Path to the pedigree input file.
        output_file (str): Path for saving the output Excel file.
        upload_path (str, optional): Path for uploading the file. Defaults to None.
    """
    print("Reading files")
    # Read sample data from CSV and pedigree data from Excel
    df = pd.read_csv(sample_file_path, encoding='latin1')
    pedg_df = pd.read_excel(pedg_path)
    
    print("Analyzing peaks")
    # Create a boolean variant dataframe and perform parallel analysis
    peak_df = bool_variant_df(df)
    peak_df = peak_df.groupby('INTERVAL_ID').parallel_apply(lambda x: get_interval_stats(x, pedg_df))
    
    # Combine interval information and analysis results
    result = pd.concat([get_info_table(df), peak_df], axis=1)
    added_result = add_dsd_distance(result)

    print("Saving")
    # Create and save an Excel file with the analysis results
    save_to_excel(added_result, pedg_df, output_file, upload_path)



def get_sample_numbers(pedg_df, source=None):
    df = pedg_df  if source == None else pedg_df[pedg_df.source == source]
    result = {'total': df.shape[0],'probands' :(~df.fam_relation.astype(bool)).sum()}
    return result

def create_sample_dict(pedg_df,columns):
    key_template = '%s\n (n=%s n_prob=%s)'
    
    numbers = get_sample_numbers(pedg_df)
    sample_dict = {key_template % ('total', numbers['total'], numbers['probands']) : [i for i in columns if 'total' in i]} 
    for sample in  pedg_df.source.unique().tolist():
        numbers = get_sample_numbers(pedg_df, sample)
        sample_dict.update({key_template % (sample, numbers['total'], numbers['probands']) : [i for i in columns if sample in i]})
    return sample_dict

def save_to_excel(result, pedg_df, output_file, upload_path):
    interval_dict = {'Peak' : ['CHROM','from','to','length'],
                    'Gene data': ['distance_from_nearest_DSD_TSS','DSDgenes_1.5mb','geneHancer','GHid','GH_is_elite','GH_type']}
    interval_dict.update(create_sample_dict(pedg_df, result.columns))
    gf.create_excel(result, interval_dict, output_file, upload_path=upload_path)


def create_folders_if_not_exist(file_path):
    """
    Create folders in the provided file path if they do not exist.

    Parameters:
    file_path (str): The file path including folders and file name.

    Returns:
    None
    """
    # Extract the directory path from the file path
    directory_path = os.path.dirname(file_path)
    
    # Create the directory path and its subdirectories if they do not exist
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python hot_peaks_table.py.py VAR_CSV_PATH SAMPLE_MATADTA_PATH OUTPUT_PATH [UPLOAD_PATH]")
        sys.exit(1)

    # Get the command line arguments
    sample_file_path = sys.argv[VAR_CSV_PATH]
    pedg_path = sys.argv[SAMPLE_MATADTA_PATH]
    output_path = sys.argv[OUTPUT_PATH]
    upload_path = sys.argv[UPLOAD_PATH] if len(sys.argv) > UPLOAD_PATH else ""

    # Call the main function with the provided arguments
    main(sample_file_path, pedg_path, output_path, upload_path)