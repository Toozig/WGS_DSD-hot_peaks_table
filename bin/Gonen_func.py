import pandas as pd
import numpy as np
import os
import subprocess
import seaborn as sns
import random
import datetime
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=False, nb_workers=8)

DBXCLI_PATH = "~/dbxcli"


SAMPLE_IDENTIFIER = ":"

def __run_bash_command_with_values(values, bash_command):
    # Check if the number of values matches the number of placeholders
    num_placeholders = bash_command.count('%s')
    if len(values) != num_placeholders:
        raise ValueError(f"Number of values provided ({len(values)}) does not match the number of placeholders ({num_placeholders}) in the bash command.")

    # Substitute placeholders with values in the bash command
    formatted_command = bash_command % tuple(values)

    # Run the formatted command using subprocess
    try:
        subprocess.run(formatted_command, shell=True, check=True)
        print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error executing command:", e)

        
    
## dropbox command ##

def upload_to_dropbox(file_path_local,file_path_dropbox):
    """
    uploads the files to dropbox
    """
    file_name = file_path_local.split('/')[-1]
    target = f"{file_path_dropbox}/{file_name}"
    command = f"{DBXCLI_PATH} put %s %s"
    __run_bash_command_with_values([file_path_local, target], command)
    
    
def download_from_dropbox(file_path_local,file_path_dropbox):
    """
    download the files file from dopbox
    """
    file_name = file_path_dropbox.split('/')[-1]
    target = f"{file_path_local}/{file_name}"
    command = f"{DBXCLI_PATH} get %s %s"
    __run_bash_command_with_values([file_path_dropbox, file_path_local], command)
    

## EXCEL command ##



def __generate_pastel_palette(n):
    pastel_palette = sns.color_palette("pastel", n, as_cmap=True)
    return pastel_palette

def __create_mini_df(df, columns):
    columns =  [col for col in columns if col in df.columns]
    cur_df = df[columns]
    # cur_df.style =  cur_df.style.set_properties(**{'background-color': generate_random_pastel_color()})
    return cur_df

def __create_sample_dict_index(samples) -> dict:
    sam_dict =  {i: [f"{i}:GT",f"{i}:DP",f"{i}:GQ"] for i in samples}
    return sam_dict


def get_date(date_format="%d%m%y"):
    current_date = datetime.datetime.now().strftime(date_format)
    return current_date


def create_excel(df, columns_dict, output_path, color_n=2, upload_path=None):
    # Generate mini DataFrames based on columns_dict
    cur_dict = {i: __create_mini_df(df, columns_dict[i]) for i in columns_dict.keys()}
    # Concatenate mini DataFrames horizontally
    df = pd.concat(cur_dict, axis=1)
    
    # Generate a pastel color palette
    pastel_palette = __generate_pastel_palette(color_n)
    
    # Create a dictionary to assign colors to level 1 columns
    level_1_dict = {list(columns_dict.keys())[i]: pastel_palette[i % color_n] for i in range(len(columns_dict.keys()))}
    
    color_dict = {}
    i = 0
    
    # Assign colors based on the level 1 column values
    for k, v in df.columns:
        color_dict.update({v: f"background-color: {level_1_dict[k]}"})

    # Define a function to apply color styling
    def mycolor(x):
        return pd.Series(color_dict)
    
    # Generate the output filename using the current date
    output = f"{output_path}_{get_date()}.xlsx"
    
    # Apply color styling and save DataFrame to an Excel file
    styled = df.style.apply(mycolor, axis=1)
    create_folders_if_not_exist(output)
    styled.to_excel(output)
    
    # Upload the file to Dropbox if upload_path is provided
    if upload_path is not None:
        upload_to_dropbox(output, upload_path)
    return styled



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


def create_excel_samples(df, columns_dict, samples, output_path, color_n=2, upload_path=None):
    """
    Create Excel when there are samples columns in it.
    
    Args:
        df (pd.DataFrame): The DataFrame to be used.
        columns_dict (dict): A dictionary defining the columns to include in the output.
        samples (list): List of sample identifiers to consider.
        output_path (str): Path where the Excel file will be saved.
        color_n (int): Number of colors for styling (optional). 
        upload_path (str): Path for uploading the file (optional).
    """
    # Get the set of samples present in the DataFrame columns
    df_samples = {i.split(SAMPLE_IDENTIFIER)[0] for i in df.columns}
    
    # Convert the samples list into a set
    samples = set(samples)
    
    # Find the samples that exist both in df and in the samples list
    exist_sample = df_samples & samples
    
    # Find the samples that do not exist in the df
    not_exist = samples - df_samples
    
    # If there are samples that do not exist in df, print a message
    if len(not_exist):
        print(f"The samples {' '.join(not_exist)} do not exist in the DataFrame.")
    
    # Update columns_dict to include only existing samples
    columns_dict.update(__create_sample_dict_index(exist_sample))
    
    # Call the internal __create_excel function to generate the Excel file
    __create_excel(df, columns_dict, output_path, color_n, upload_path)

    
    
    