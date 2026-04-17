'''
Function definitions necessary for the peptide annotation.  
'''
import pandas as pd
from Bio import SeqIO
import re


def get_peptide_lengths(peptides_df: pd.DataFrame) -> list[int]:
    '''Return all peptide lengths.

    Args:
        peptides_df: Pandas dataframe containing all peptides.

    Returns:
        A list containing all peptide lengths in the input.
    '''
    return list(set(peptides_df['sequence'].apply(lambda seq: len(seq)).to_list()))


def proteome_to_dict(proteome: str) -> dict[str,str]:
    '''Read reference proteome into dictionary.
    
    Args: 
        proteome: Path to the reference proteome.

    Returns:
        The reference proteome as a dictionary.
    '''
    proteome_dict = {}
    proteome = SeqIO.parse(open(proteome),'fasta')
    for protein in proteome:
        proteome_dict[protein.id] = str(protein.seq)
    return proteome_dict


def get_seq_length(seq: str, length:int) -> list[str]:
    '''Return all sequences of a certain length that occur in the protein sequence.

    Args:
        seq: Protein sequence.
        length: Length of peptides.

    Returns:
        List of all subsequences of a specified length.
    '''
    seq_length = len(seq)
    seqs = []
    for i in range(seq_length-length+1):
        seqs.append(seq[i:i+length])
    return seqs


def create_fasta_dict(proteome_dict: dict, lengths: list[int]) -> dict:
    '''Create dictionary containing the peptides as keys and the accessions as values.

    Args:
        proteome_dict: Dictionary containing the proteome.
        lengths: Lengths of peptides that should be in the dictionary.

    Returns:
        A dictionary containing all peptides as keys and the accessions as values.
    '''
    lookup_dict = {}
    for length in lengths:
        for acc, seq in proteome_dict.items():
            seqs = get_seq_length(seq, length)
            for pep in seqs:
                if pep in lookup_dict:
                    lookup_dict[pep] += f';{acc}'
                else:
                    lookup_dict[pep] = acc
    return lookup_dict


def annotate_peptides(fasta_dict: dict, peptide_df: pd.DataFrame) -> pd.DataFrame:
    '''Annotates peptides in peptide_csv.

    Args:
        fasta_dict: Dictionary containing peptides as keys and accessions as values.
        peptide_csv: Path to file containing the peptides to be annotated
    
    Returns:
        The input dataframe with protein annotations.
    '''
    peptide_df['accessions'] = peptide_df['sequence'].apply(lambda seq: fasta_dict[seq] if seq in fasta_dict.keys() else 'unmapped')
    return peptide_df


def get_pos(accessions: list[str], peptide: str, fasta_dict: dict) -> list[list]:
    '''Compute positions of peptide in sequence.

    Args:
        accessions: List of protein accessions to which peptide is mapped.
        peptide: Peptide sequence of interest.
        fasta_dict: Dictionary containing reference proteome.

    Returns:
        Lists containing the accessions and start and end position of a peptide in the proteome.
    '''
    peptide_len = len(peptide)
    peptide = peptide.replace('I','L').replace('L','(I|L)')
    starts = ''
    ends = ''
    accessions_str = ''
    for accession in accessions.split(';'):
        if accession != 'unmapped':
            seq = fasta_dict[accession]
            groups_pos = re.finditer(f'(?=({peptide}))', seq)
            for group_pos in groups_pos:
                pep_start, pep_end = group_pos.span()
                accessions_str += f';{accession}'
                starts += f';{pep_start}'
                ends += f';{pep_end+peptide_len-1}'
        else:
            accessions_str = f';unmapped'
    return accessions_str[1:], starts[1:], ends[1:]


def add_positions(fasta_dict: dict, annotated_df: pd.DataFrame, accession_column: str, seq_column: str) -> pd.DataFrame:
    '''Adds positions of peptides in the proteins.

    Args:
        fasta_dict: Dictionary containing peptides as keys and accessions as values.
        annotated_csv: Path to csv containing the peptide sequence and the proteins they map to.
        accession_column: Header of the column containing the protein accessions.
        seq_column: Header of the column containing the peptide sequences.

    Returns:
        The input dataframe with annotations.
    '''

    annotated_df[['accessions', 'start', 'end']] = annotated_df.apply(lambda row: get_pos(row[accession_column], row[seq_column], fasta_dict), result_type='expand', axis=1)
    return annotated_df


