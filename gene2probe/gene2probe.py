import pandas as pd
import numpy as np
import pybedtools
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import subprocess

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


def get_region_of_interest(gtf, gene_id, gene_id_type, feature, distance_from_exons=100):
    """
    Subset a gene annotation in gtf format for the gene and feature of interest.
    """
    ## If gene_id_type not in gene_id, gene_name raise error:
    if gene_id_type not in ['gene_name', 'gene_id']:
        raise ValueError("gene_id_type should be one of ['gene_name', 'gene_id'].")
        
    ## Extract regions corresponding to gene of interest (symbol: gene_name, Ensembl ID: gene_ID)
    gene_ids = gtf['attribute'].apply(extract_feature_from_gtf, feature=gene_id_type)

    ## If gene_id not in gtf, raise error:
    if gene_id not in gene_ids.values:
        raise ValueError("The gtf file should contain the provided gene ID.")

    if feature not in ['exon', 'CDS', 'transcript', 'intron']:
        raise ValueError("Feature has to be one of ['exon', 'CDS', 'transcript', 'intron', 'splice_junction']")

    ## Filter for gene_id
    roi = gtf.iloc[np.where(gene_ids == gene_id)[0],:]
        
    ## If feature is exon, CDS or transcript, our job is easy:
    if feature in ['exon', 'CDS', 'transcript']:
        roi = roi[roi['feature']==feature]
        ## Convert to bed-style df
        roi_bed = gtf_2_bed(roi, name_pref = (gene_id + '_'))
    ## If feature is 'intron', we need a bit more work, which is why we use a dedicated function:
    elif feature == 'intron':
        roi_bed = get_introns(roi, dist_from_exon = distance_from_exons)
        roi_bed.loc[:, 'name'] = gene_id + '_' + roi_bed.index.astype(str)

    ## If feature is 'splice-junction', we need a completely different data structure
    elif feature == 'splice_junction':
        ## Split by transcript
        transcript_gtfs = split_by_transcript(roi)
        junction_dfs =  []
        ## For each transcript, get splice junctions
        for tx in transcript_gtfs.keys():
            junction_dfs.append(get_splice_junctions(transcript_gtfs[tx]))
        ## Then concatenate into a single dataframe
        roi_bed = pd.concat(junction_dfs)
        
    ## Return the region of interest
    return roi_bed

def get_introns(gtf, dist_from_exon=100):
    """
    Given a gtf with a single gene, identify the introns by subtracting the exons from the transcript regions (and additionally add some distance from the exons)
    """
    ## Check that only a single gene_id is provided
    gene_ids = gtf['attribute'].apply(extract_feature_from_gtf, feature='gene_name')
    if len(gene_ids.unique()) > 1:
        raise ValueError("The gtf file should contain a single gene name.")
    
    ## Check that we have both 'transcript' and 'exon' features
    if (sum(gtf['feature']=='transcript')==0) or (sum(gtf['feature']=='exon')==0):
        raise ValueError("The gtf file should contain both transcripts and exons to be able to detect introns.")
    
    transcripts = gtf[gtf['feature']=='transcript']
    exons = gtf[gtf['feature']=='exon']
    ## Add flanking region to exons
    exons.loc[:,'start'] = exons.loc[:,'start'] - dist_from_exon
    exons.loc[:,'end'] = exons.loc[:,'end'] + dist_from_exon
    
    ## Convert to bed
    transcripts_bed = pybedtools.BedTool.from_dataframe(gtf_2_bed(transcripts))
    exons_bed = pybedtools.BedTool.from_dataframe(gtf_2_bed(exons))
    
    ## Identify introns as transcript regions not overlapping 
    introns_bed = transcripts_bed.subtract(exons_bed)

    ## Return as a bed-like dataframe (expected)
    return introns_bed.to_dataframe(names=['seqname', 'start', 'end', 'name', 'score', 'strand'])

def split_by_transcript(gtf):
    """
    Given a single gene gtf, it extracts transcript ids, returns a dictionary in the form {'transcript': 'exon_gtf'}
    """
    ## Check that only a single gene_id is provided
    gene_ids = gtf['attribute'].apply(extract_feature_from_gtf, feature='gene_name')
    if len(gene_ids.unique()) > 1:
        raise ValueError("The gtf file should contain a single gene name.")
    ## Identify transcripts
    transcript_ids = gtf['attribute'].apply(extract_feature_from_gtf, feature='transcript_id')
    ## Create a dictionary of gtfs with tx_id as a key and the gtf for each transcript
    transcripts_gtf = {}
    for tx_id in transcript_ids.unique():
        transcripts_gtf[tx_id] = gtf.iloc[np.where(transcript_ids == tx_id)[0],:]
    return transcripts_gtf

def get_splice_junctions(gtf, flank = 35):
    """
    Given a gtf for a single transcript, stitch together consecutive exons
    """
    ## Check that only a single gene_id is provided
    transcript_ids = gtf['attribute'].apply(extract_feature_from_gtf, feature='transcript_id')
    if len(transcript_ids.unique()) > 1:
        raise ValueError("The gtf file should contain a single transcript id.")
    transcript_id = transcript_ids.iloc[0]
    if sum(gtf['feature']=='exon')<2:
        raise ValueError("The gtf file should contain at least two exons.")
    ## Subset to exons only and convert to bed
    exons_bed = gtf_2_bed(gtf.loc[gtf['feature']=='exon',:], name_pref=(transcript_id + '_exon_'))

    junction_dfs = []
    ## For each exon, get the rightmost region of length=flank and concatenate it (row-wise) to the leftmost region of length=flank of the next exon
    ## We also make
    for i in range(exons_bed.shape[0] - 1):
        junction_df = pd.DataFrame({'seqname': exons_bed['seqname'].iloc[0],
                                'start': [max((exons_bed['end'].iloc[i] - flank),exons_bed['end'].iloc[i]), exons_bed['start'].iloc[i+1]],
                                'end': [exons_bed['end'].iloc[i], min((exons_bed['start'].iloc[i+1] + flank) + exons_bed['start'].iloc[i+1])],
                                'name': (transcript_id + '_splice_junction_' + str(i)),
                                'score': '.',
                                'strand': exons_bed['strand'].iloc[0]})
        junction_dfs.append(junction_df)

    ## Concatenate dfs
    junctions_bed = pd.concat(junction_dfs)  
    return junctions_bed 

def define_region_of_interest(coord, strand, name=None):
    """
    Generate a bed file from a string of the form chrX:123-456 and the strand
    """
    chr, ranges = coord.split(':')
    # Split the ranges into start and end
    start, end = ranges.split('-')
    # If no name was provided, use coord:
    if name is None:
        name=coord
    # Convert start and end into integers
    roi_bed = pd.DataFrame({'seqname': [chr], 'start': [int(start)], 'end': [int(end)], 'name': [name], 'score': ['.'], 'strand': [strand]})
    return roi_bed

def define_custom_splice_junction(coord_1, coord_2, strand):
    """
    Given a pair of genomic coordinates, and their strand, stitch them together as a splice junction
    """
    exon_1 = define_region_of_interest(coord_1, strand, name=(coord_1 + '_' + coord_2 + '_junction'))
    exon_2 = define_region_of_interest(coord_2, strand, name=(coord_1 + '_' + coord_2 + '_junction'))
    # Run a couple of sanity checks:
    if exon_1['seqname'].iloc[0] != exon_2['seqname'].iloc[0]:
        raise ValueError("Chromosomes must match between exons")
    if exon_1['strand'].iloc[0] != exon_2['strand'].iloc[0]:
        raise ValueError("Strands must match between exons")
        
    junction_bed = pd.concat([exon_1, exon_2]).reset_index(drop=True)
    return(junction_bed)

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
        seq_df = seq_df[(seq_df[col] > min_GC) & (seq_df[col] < max_GC)].copy()

    ## Reset index and return df
    seq_df = seq_df.reset_index(drop=True)
    return seq_df

def run_blast(fasta, blastdb, path2blastn, outfile):
    """
    Given a set of sequences in fasta format, and a blast database, run a very lenient version of blastn (in strand-aware mode) and read in the result.
    """
    ### First, specify the BLAST command
    command = [
        path2blastn,
        '-query', fasta, ## Sequences of targeted regions
        '-db', blastdb,                     # BLAST database
        '-out', outfile,  # Output file
        '-outfmt', '6',
        '-task', 'blastn',
        '-strand', 'plus',  # Make sure that we are BLASTing against transcripts in the same strand
        '-evalue', '1',  # Very lenient e-value to detect even distant potential off-targets
        '-dust', 'no'       # Turn off low-complexity filter
    ]

    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)

    ## And read in the blast output
    blast_columns = [
        "name",   # Query Seq-id
        "sseqid",   # Subject Seq-id
        "pident",   # Percentage of identical matches
        "length",   # Alignment length
        "mismatch", # Number of mismatches
        "gapopen",  # Number of gap openings
        "qstart",   # Start of alignment in query
        "qend",     # End of alignment in query
        "sstart",   # Start of alignment in subject
        "send",     # End of alignment in subject
        "evalue",   # Expect value
        "bitscore"  # Bit score
    ]
    
    blast_res = pd.read_csv(
        outfile,
        sep='\t', 
        header=None, 
        names=blast_columns
        )

    ## Extract the gene ID of the target:
    blast_res['sgeneid'] = blast_res['sseqid'].str.split('::').str[0]

    ## And return the result
    return blast_res

def detect_offtargets(blast_res, gene_id, min_mismatches = 5):
    """
    Given a dataframe of blast results, which includes a 'sgeneid' column, detect offtargets as alignments of length >= max_alignment with genes different than the target
    """
    offtargets = blast_res['name'][(blast_res['sgeneid']!=gene_id) & (blast_res['mismatch'] < min_mismatches)].unique().tolist()
    return offtargets

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
