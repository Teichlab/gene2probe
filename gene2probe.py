import pandas as pd
import numpy as np
import pybedtools
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

def read_gtf(path2gtf):
    """
    Read gtf file into pandas
    """
    df = pd.read_csv(path2gtf, 
                     comment='#', 
                     header=None, 
                     sep='\t',
                     names=[
    'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'
    ])
    return df

def extract_feature_from_gtf(attributes, feature='gene_id'):
    """
    Extract feature from gtf's 9th column
    """
    parts = attributes.split(';')
    gene_id = [part.strip().split(' ')[1].strip('"') for part in parts if feature in part]
    return gene_id[0] if gene_id else None


def gtf_2_bed(gtf, name_pref = 'region_'):
    """
    Convert a pandas dataframe in gtf format into one in bed6 format
    """
    gtf_copy = gtf.copy().reset_index(drop=True)
    ## Create name field
    gtf_copy['name'] = name_pref + gtf_copy.index.astype(str)
    ## Subset and rearrange columns as in bed
    df_bed = gtf_copy[['seqname', 'start', 'end', 'name', 'score', 'strand']]
    # Adjust start positions from 1-based GTF to 0-based BED
    df_bed['start'] = df_bed['start'].astype(int) - 1
    return df_bed


def generate_50mers(bed_df):
    """
    Generate all possible 50-mer intervals from a BED-format DataFrame while preserving all original columns
    and appending a suffix to the name for each k-mer.

    Parameters:
    bed_df (pd.DataFrame): DataFrame containing columns 'seqname', 'start', 'end', and optionally 'name', 'score', 'strand'.

    Returns:
    pd.DataFrame: A new DataFrame with the same columns for each 50-mer.
    """
    all_50mers = []

    # Iterate through each row in the DataFrame
    for index, row in bed_df.iterrows():
        seqname = row['seqname']
        start = row['start']
        end = row['end']
        name = row.get('name', 'region')  # Default name if none provided
        score = row.get('score', '.')  # Default score if none provided
        strand = row.get('strand', '.')  # Default strand if none provided

        # Generate 50-mer intervals within the range from start to end
        # Ensure that each interval is exactly 50 bp long
        kmer_index = 0
        for pos in range(start, end - 49):  # Subtract 49 because we want intervals of exact length 50
            kmer_name = f"{name}_{kmer_index}"
            all_50mers.append({
                'seqname': seqname,
                'start': pos,
                'end': pos + 50,
                'name': kmer_name,
                'score': score,
                'strand': strand
            })
            kmer_index += 1

    # Convert the list of dictionaries to a DataFrame
    mers_df = pd.DataFrame(all_50mers)
    return mers_df

# Example usage:
# Assuming you have a DataFrame 'bed_df' loaded with your BED data
# bed_df = pd.DataFrame({
#     'seqname': ['chr1', 'chr1'],
#     'start': [100, 200],
#     'end': [160, 300],
#     'name': ['gene1', 'gene2'],
#     'score': ['.', '.'],
#     'strand': ['+', '-']
# })
# mers_df = generate_50mers(bed_df)
# print(mers_df)

def get_longest_homopolymer(seq):
    """ Returns the length of the longest homopolymer in the sequence. """
    max_count = 1
    current_count = 1
    last_char = ''
    
    for char in seq:
        if char == last_char:
            current_count += 1
            max_count = max(max_count, current_count)
        else:
            current_count = 1
            last_char = char
    
    return max_count


def remove_overlapping_probes(probe_df, probe_id, offset = 100):
    """ Remove all probes in the dataframe that overlap the selected probe"""
    i = np.where(probe_df['name']==probe_id)[0]
    start = int(probe_df['start'].iloc[i])
    end = int(probe_df['end'].iloc[i])
    probe_df_filtered = probe_df[(probe_df['end'] < start - offset) | (probe_df['start'] > end + offset)].copy()
    return(probe_df_filtered)

def write_fasta(names, sequences, output_file):
    if len(names) != len(sequences):
        raise ValueError("The length of names and sequences must be the same.")

    with open(output_file, 'w') as f:
        for name, seq in zip(names, sequences):
            f.write(f">{name}\n")  # Write the header
            f.write(f"{seq}\n")    # Write the sequence