# Variant Analysis Script

## Description

This script performs analysis on variant data and pedigree information, calculating various statistics for intervals and samples. It then saves the analysis results to an Excel file. The script takes input files for variant data and pedigree information, and optionally an upload path for the output Excel file.

## Features

- Analyzes variant data and pedigree information.
- Calculates statistics for intervals and samples.
- Saves analysis results to an Excel file.
- Supports parallel processing for faster analysis.

## Prerequisites

- Pandarallel library (for parallel processing)

## Usage

1. Clone the repository:

    ```
    git clone https://github.com/yourusername/variant-analysis-script.git
    cd variant-analysis-script
    ```

2. Install Pandarallel (if not installed):

    ```
    pip install pandarallel
    ```

3. Run the script:

    ```
    python script_name.py VAR_CSV_PATH SAMPLE_METADATA_PATH OUTPUT_PATH [UPLOAD_PATH]
    ```

   Replace `script_name.py` with the actual script name.

   Arguments:
   - `VAR_CSV_PATH`: Path to the CSV file containing variant data.
   - `SAMPLE_METADATA_PATH`: Path to the Excel file containing sample metadata.
   - `OUTPUT_PATH`: Path to the directory where the output Excel file will be saved.
   - `UPLOAD_PATH` (optional): Path for uploading the output Excel file.

## Example Usage

1. Basic usage without specifying upload path:

    ```
    python script_name.py data/pipeline_outputs/variants_with_layers/qualityDSD_variants.csv data/read_only/samples/sample_metadata.xlsx ~/bcftools_prog/DSD_project/data/hot_peaks_fixed
    ```

2. Usage with specifying an upload path:

    ```
    python script_name.py data/pipeline_outputs/variants_with_layers/qualityDSD_variants.csv data/read_only/samples/sample_metadata.xlsx ~/bcftools_prog/DSD_project/data/hot_peaks_fixed Nitzan_Gonen_lab/Joint_projects/WGS_on_DSD/Ido/hot_peaks
    ```

I tried to fork!

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
