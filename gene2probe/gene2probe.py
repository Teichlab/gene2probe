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
    df_bed.loc[:, 'start'] = df_bed['start'].astype(int) - 1
    return df_bed


def generate_kmers(bed_df, k):
    """
    Generate all possible 50-mer intervals from a BED-format DataFrame while preserving all original columns
    and appending a suffix to the name for each k-mer.

    Parameters:
    bed_df (pd.DataFrame): DataFrame containing columns 'seqname', 'start', 'end', and optionally 'name', 'score', 'strand'.

    Returns:
    pd.DataFrame: A new DataFrame with the same columns for each 50-mer.
    """
    all_kmers = []

    # Iterate through each row in the DataFrame
    for index, row in bed_df.iterrows():
        seqname = row['seqname']
        start = row['start']
        end = row['end']
        name = row.get('name', 'region')  # Default name if none provided
        score = row.get('score', '.')  # Default score if none provided
        strand = row.get('strand', '.')  # Default strand if none provided

        # Generate k-mer intervals within the range from start to end
        # Ensure that each interval is exactly k bp long
        kmer_index = 0
        for pos in range(start, end - (k-1)):  # Subtract k-1 because we want intervals of exact length 50
            kmer_name = f"{name}_{kmer_index}"
            all_kmers.append({
                'seqname': seqname,
                'start': pos,
                'end': pos + k,
                'name': kmer_name,
                'score': score,
                'strand': strand
            })
            kmer_index += 1

    # Convert the list of dictionaries to a DataFrame
    kmers_df = pd.DataFrame(all_kmers)
    return kmers_df

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

def get_sequence_stats(fasta_seq, probe_length, split_nt):
    """
    Read sequences and summarise sequence stats (GC content, longest polymer stretch )for a set of kmers.
    Inputs:
    fasta_seq: Path to fasta file (can be an actual path or obtained through pybedtools with obj.seqfn)
    probe_length (int): length of the probe.
    split_nt (int): 0-based index of where the RHS starts - set to None if probe is not split.

    Returns:
    seq_stats (pd.DataFrame): one row per kmer with the sequences and the stats for each one.
    """
    ## Initialise list of stats to keep track of
    kmer_coord = []
    transcript_seq = []
    probe_seq = []
    gc_content = []
    longest_homopol = []
    if split_nt is not None:
        gc_content_lhs = []
        gc_content_rhs = []

    ## Read fasta file line by line
    with open(fasta_seq) as f:
        for line in f:
            if line.startswith('>'):
                ## Get name of fasta sequence (coordinates)
                kmer_coord.append(line.strip().replace('>', ''))
            else:
                ## Get sequence
                kmer_seq = Seq(line.strip())
                transcript_seq.append(str(kmer_seq))
                ## Reverse complement - this is the sequence of the probe
                probe_seq.append(str(kmer_seq.reverse_complement()))
                ## Estimate GC content (whole probe)
                gc_content.append(gc_fraction(kmer_seq))
                ## If this is a split probe, also estimate the GC content in each half
                if split_nt is not None:
                    gc_content_rhs.append(gc_fraction(kmer_seq[0:split_nt])) ## RHS of the probe is LHS of the target - GC content is the same
                    gc_content_lhs.append(gc_fraction(kmer_seq[split_nt: probe_length])) ## RHS of the probe is LHS of the target - GC content is the same
                ## Estimate longest homopolymer
                longest_homopol.append(get_longest_homopolymer(kmer_seq))

    ## Summarise everything into a dataframe
    seq_stats = pd.DataFrame({
        'kmer_coord': kmer_coord,
        'transcript_seq': transcript_seq,
        'probe_seq': probe_seq,
        'GC_content_full':gc_content,
        'longest_homopolymer': longest_homopol
    })
    if split_nt is not None:
        seq_stats['GC_content_LHS'] = gc_content_lhs
        seq_stats['GC_content_RHS'] = gc_content_rhs
    return seq_stats

def check_for_required_nts(seq_df, required_nts):
    """
    Checks whether kmers have required nts in specified position

    Inputs:
    seq_df (pd.DataFrame): Dataframe of kmers - probe_seq should be a column.
    required_nts (dict): Dictionary with 0-based index (int) as key and nucleotide value (str) as value

    Output:
    has_required_nts (list): List of length = nrow of seq_df - boolean value whether a kmer satisfies nucleotide requirements or not.
    """
    ## Initially mark all probes as matching nucleotide requirement:
    has_required_nts = [True for i in range(seq_df.shape[0])]
    ## For each requirement:
    for nt_idx in required_nts.keys():
        ## For each k-mer:
        for i in range(seq_df.shape[0]):
            ## If it failes the requirement:
            if seq_df['probe_seq'][i][nt_idx] != required_nts[nt_idx]:
                ## Mark as False
                has_required_nts[i] = False
    ## Return
    return has_required_nts


def filter_by_GC_content(seq_df, min_GC, max_GC):
    """
    Filters a set of kmers based on minimum and maximum GC content ranges.

    Inputs:
    seq_df (pd.DataFrame): Dataframe of kmers - should contain at least one column starting with GC_content.
    min_GC (float): Minimum fraction of GC content (should be between 0 and 1).
    max_GC (float): Maximum fraction of GC content (should be between 0 and 1).
    """
    ## First identify all relevant columns (depending on whether we have a split probe or not):
    gc_columns = [col for col in seq_df.columns if col.startswith('GC_content')]

    ## For each column, filter for min/max GC fraction
    for col in gc_columns:
        seq_df = seq_df[(seq_df[col > min_GC]) & (seq_df[col < min_GC])].copy()

    ## Reset index and return df
    seq_df = seq_df.reset_index(drop=True)
    seq_df

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
