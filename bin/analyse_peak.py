"""
This script take a variant DF and output peak analysis.
note that the script does not filterring the table. hence the input should be filtered

"""
import numpy as np
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=8)
import Gonen_func


def get_interval_stats(df_in, pedg_df):
    df = df_in.drop(columns='INTERVAL_ID')
    sum_per_sample = df.sum(axis=0)
    exist_per_sample = df.any(axis=0)
    sum_series = pd.Series(dtype=float)
    probands =  pedg_df[pedg_df.fam_relation == 0].ID.astype(str)
    sum_series['total n probands'] = exist_per_sample[exist_per_sample.index.isin(probands)].sum()
    sum_series['total n samples'] = exist_per_sample.sum()
    sum_series['total n proband variants'] = sum_per_sample[sum_per_sample.index.isin(probands)].sum()
    sum_series['total n variants'] = sum_per_sample.sum()
    for label in pedg_df.source.unique():
        label_id =  pedg_df[pedg_df.source == label].ID.astype(str)
        probans_id = label_id[label_id.isin(probands)]
        sum_series[f'{label} n samples'] = exist_per_sample[exist_per_sample.index.isin( label_id)].sum()
        sum_series[f'{label} n probands'] = exist_per_sample[exist_per_sample.index.isin( probans_id)].sum()
        sum_series[f'{label} n variants'] = sum_per_sample[sum_per_sample.index.isin(label_id)].sum()
        sum_series[f'{label} n proband variants'] = sum_per_sample[sum_per_sample.index.isin(probans_id)].sum()

    return sum_series

def bool_variant_df(df):
    """
    returns a boolean table of samples VS variants.
    True if the variant exist in the sample.
    """
    gt = df[[i for i in df.columns if i.endswith('GT')]].replace(' ',np.nan).replace('',np.nan).notna()
    gt.columns = [i.replace(':GT','') for i in gt.columns ]

    return pd.concat([df[['INTERVAL_ID']],gt], axis=1)

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
                    'Gene data': ['DSDgenes_1.5mb','geneHancer','GHid','GH_is_elite','GH_type']}
    interval_dict.update(create_sample_dict(pedg_df, result.columns))
    Gonen_func.create_excel(result, interval_dict, output_file, upload_path=upload_path)


def get_info_table(df):
    relevant_coulmns = ['CHROM','from', 'to', 'length', 'DSDgenes_1.5mb','geneHancer', 'GHid', 'GH_is_elite', 'GH_type']
    df = df[~df.INTERVAL_ID.duplicated()].set_index('INTERVAL_ID')
    return df[relevant_coulmns]
    

def main(sample_file_path, pedg_path, output_file, upload_path=""):
    print( "reading files")
    df = pd.read_parquet(sample_file_path)
    pedg_df= pd.read_excel(pedg_path)
    print("analysing peaks")
    ## getting a table that tells if a sample have variant or not
    peak_df = bool_variant_df(df)
    ## getting the data on each interval
    peak_df = peak_df.groupby('INTERVAL_ID').parallel_apply(lambda x: get_interval_stats(x,pedg_df))
    ## adding adtional peak info 
    result  = pd.concat([get_info_table(df),peak_df],axis=1)
    print("saving")
    save_to_excel(result, pedg_df, output_file, upload_path)



    

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python analyse_peak.py sample_file_path pedg_path output_file [upload_path]")
        sys.exit(1)

    sample_file_path = sys.argv[SAMPLE_FILE_IDX]
    pedg_path = sys.argv[SAMPLE_FILE_IDX]
    output_file = sys.argv[OUTPUT_PATH_IDX]
    upload_path = sys.argv[UPLOAD_PATH_IDX] if len(sys.argv) >= ARG_MAX_NUM else ""

    main(sample_file_path, pedg_path, output_file, upload_path)