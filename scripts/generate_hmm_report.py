import pandas as pd
import numpy as np
import argparse
import os
from Bio import SeqIO
from typing import List, Optional, Dict
import logging


# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Script to generate reports from the result file of the script 'analyze_seq.py'.")

    parser.add_argument('-i', '--input', type=str, required=True, help='Input file or regular expression of file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output folder')
    parser.add_argument('-m', '--hmm_info', type=str, required=True, help='Folder with HMM drawings')

    return parser.parse_args()


def check_file_exists(file_path: str, message: str) -> None:
    """Check if a file exists, otherwise log an error message and exit."""
    if not os.path.isfile(file_path):
        logging.error(message)
        exit(1)


def load_csv_file(file_path: str, index_col: int = 0, sep=',') -> pd.DataFrame:
    """Load a CSV file into a DataFrame."""
    return pd.read_csv(file_path, index_col=index_col, sep=sep)

def export_seq_data(seq: str) -> List[str]:
    """Extract and format sequence data."""
    try:
        name, chain, part, length, length_contig, skew, _ = seq.split(';')
        chain = chain.replace('chain_', '')
        part = part.replace('part_', '')
        length = length.replace('len_', '')
        skew = skew.replace('skew_', '')
        length_contig = length_contig.replace('len_contig_', '')
        frame = int(chain.replace('+1', '1').replace('-1', '4')) + int(skew)
        return [name, part, length, frame, length_contig]
    except ValueError as e:
        logging.error(f"Error parsing sequence data: {e}")
        return [None] * 5  # Return a list of None values to maintain DataFrame structure


def get_coordinates(dataframe: pd.DataFrame) -> np.ndarray:
    """Calculate coordinates based on frame and part using vectorized operations."""
    frame_condition = dataframe['Frame'] <= 3
    from_coords = np.where(
        frame_condition,
        dataframe['From'] * 3 + dataframe['Part'],
        (dataframe['Lengh'] - dataframe['From'] * 3) + dataframe['Part']
    )
    to_coords = np.where(
        frame_condition,
        dataframe['To'] * 3 + dataframe['Part'],
        (dataframe['Lengh'] - dataframe['To'] * 3) + dataframe['Part']
    )
    return np.column_stack((from_coords, to_coords))


def create_report(input_df: pd.DataFrame, hmm_df: pd.DataFrame) -> pd.DataFrame:
    """Create a report by merging input data with HMM info."""
    report = input_df.merge(hmm_df, left_on='Query', right_on='Seed_RdRp_HMM Name', how='left').drop('Seed_RdRp_HMM Name', axis=1)
    report = report.sort_values('Score', ascending=False)

    try:
        seq_data = np.array(report.Name.apply(export_seq_data).to_list())
        report[['Name', 'Part', 'Lengh', 'Frame', 'Length_contig']] = seq_data
        report[['Lengh', 'From', 'To', 'Part', 'Frame', 'Length_contig']] = report[
            ['Lengh', 'From', 'To', 'Part', 'Frame', 'Length_contig']
        ].astype('int')
        report[['From', 'To']] = get_coordinates(report)
        report['Alignment_Length'] = np.abs(report['To'] - report['From'])
        report = report.drop(['Part', 'Lengh'], axis=1)
    except Exception as e:
        logging.error(f"An error occurred while processing the report: {e}")

    return report


def main() -> None:
    """Main function to execute the script."""
    args = parse_arguments()

    # Validate input files
    check_file_exists(args.input, "An input file does not exist")
    check_file_exists(args.hmm_info, "HMM info folder does not exist")

    # Load data
    input_file = load_csv_file(args.input)
    input_file['Query'] = input_file['Query'].apply(lambda x: x.replace('.full', ''))
    hmm_info = load_csv_file(args.hmm_info, sep='\t')
    hmm_info = hmm_info[['Seed_RdRp_HMM Name', 'InterProScan_Pfam_InterPro Annotations Description', 'Palmscan_Group']]


    # Process input files
    report = create_report(input_file, hmm_info)

    # Filter and save reports


    report.reset_index(drop=True).to_csv(args.output)

    logging.info("Script executed successfully.")


if __name__ == "__main__":
    main()

