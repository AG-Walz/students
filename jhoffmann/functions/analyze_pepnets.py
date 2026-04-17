'''
Functions used in analyze_cluster.ipynb to analyze the performance of pepnets.
'''
import numpy as np
import matplotlib.pyplot as plt
import re

def group_repetitive(start: list[int], end: list[int]) -> tuple[list[int], list[int]]:
    '''Groups repetitive peptides together.

    Args:
        start: The start positions of the peptide appearances.
        end: The end positions of the peptide appearances.

    Returns:
        The start and end positions of the repetitive group.
    '''

    updated_start = []
    updated_end = []

    lists = list(zip(start, end))
    lists = sorted(lists, key=lambda x: int(x[0]))
    start, end = zip(*lists)

    updated_start.append(start[0])

    for pep_pos in range(len(start)-1):

        # two start positions are not part of one repetitive region if the next start position is higher than the current end position 
        if int(start[pep_pos + 1]) > int(end[pep_pos]):
            updated_end.append(end[pep_pos])
            updated_start.append(start[pep_pos + 1])

    # add the last occurrences end position to the end positions
    updated_end.append(end[-1])

    return list(updated_start), list(updated_end)


def position(entire: str, sequences: str, start: int, end: int) -> tuple[list[int], list[int]]:
    '''Calculates the position of sequences in the entire sequence.

    Args:
        entire: Protein sequence.
        sequences: List of peptide sequences.
        start: Start of the peptide cluster.
        end: End of the peptide cluster.
    
    Returns:
        A tuple containing the start and end positions of the sequences in the 
        entire protein sequence. If a sequence appears multiple times in the 
        protein the position closest to the cluster start and end position is
        returned.
    '''

    starts_f = []
    ends_f = []

    for seq in sequences:

        # allow for imperfect matches with more or less general amino acids (e.g. B)
        seq_r = re.sub(r'(Z|Q|E)','(Z|Q|E)',re.sub(r'(N|B|D)','(N|B|D)', seq.replace('X','*')))

        # all appearances
        seq_pos = [pos.span() for pos in re.finditer(f'(?=({seq_r}))', entire)]
        starts = [int(pos[0]) for pos in seq_pos]
        ends = [int(pos[1])+len(seq)-1 for pos in seq_pos]
        starts_u, ends_u = group_repetitive(starts , ends)

        # take occurrence closest to cluster position
        if len(starts_u) > 1:
            diff = []
            for s, e in zip(starts_u, ends_u):
                diff.append(np.sqrt((s-start[0])**2+(e-end[0])**2))
            o = diff.index(min(diff))
            starts_u = [starts_u[o]]
            ends_u = [ends_u[o]]

        for group_pos in zip(starts_u, ends_u):

            pep_start, pep_end = group_pos
            starts_f.append(pep_start)
            ends_f.append(pep_end)

    return starts_f, ends_f


def update_landscape(landscape: list[int], starts: list[int], ends: list[int]) -> list[int]:
    '''Update protein landscape.

    Args: 
        landscape: Current landscape.
        starts: All start positions of the peptides mapped to a peptide cluster.
        ends: All end positions of the peptides mapped to a peptide cluster.
    
    Returns:
        The updated peptide landscape.
    '''
    for start, end in zip(starts, ends):
        start = start - min(starts)
        end = end - min(starts)
        for i in range(start, end+1):
            landscape[i] += 1
    return landscape


def id_consensus(landscape: list[int], min_epi_len: int, entire: str) -> tuple[str, int, int]:
    '''Calculate consensus sequences as in epicore.

    Args:
        landscape: The landscape of a peptide cluster.
        min_epi_length: The minimal epitope length.
        entire: The entire cluster sequence. # TODO: make sure not the entire protein is returned here

    Returns:
        A tuple containing the consensus sequence plus it's start and end position.
    '''
    
    # get all landscape values
    total_counts = np.unique(landscape)
    total_counts[::-1].sort()
        
    # find total coverage for which consensus epitope is at least min_epi_len long
    for total_count in total_counts:

        Z = landscape < total_count

        # get lengths of peptide sequences with coverage above the current threshold
        seqs_idx = np.where(np.diff(np.hstack(([False],~Z,[False]))))[0].reshape(-1,2)
                
        # get length of longest peptide subsequences with current count
        ce_start_pos = seqs_idx[np.diff(seqs_idx, axis=1).argmax(),0]
        current_pep_length = np.diff(seqs_idx, axis=1).max()
                
        # check if min_epi_length is fulfilled for that sequence
        if current_pep_length >= min_epi_len:

            # get position of epitope in protein sequences
            pep_in_prot_start = ce_start_pos
            pep_in_prot_end = pep_in_prot_start + current_pep_length

            # get consensus epitopes
            consensus_seq = entire[pep_in_prot_start:pep_in_prot_end]
            return consensus_seq, pep_in_prot_start, pep_in_prot_end-1

    return entire, 0, len(entire)-1


def compute_coverage(pep_start: int, pep_end: int, consensus_start: int, consensus_end: int) -> int:
    '''Compute the coverage of the consensus sequence by a peptide sequence.

    Args:
        pep_start: Peptide start position.
        pep_end: Peptide end position.

    Returns:
        The coverage of the consensus sequence by the peptide sequence.
    '''

    if pep_start < consensus_start:
        overlap = max(min(consensus_end-consensus_start+1, pep_end-consensus_start+1),0)
    else:
        overlap = max(min(pep_end-pep_start+1, consensus_end-pep_start+1),0)
    return overlap/(consensus_end-consensus_start+1)


def plot_coverage(coverage_hist_all: list[float]):
    '''Plot the coverages of the consensus sequences by the peptides.

    Args:
        coverage_hist_all: List containing all coverages.
        out: Path were the plot gets saved.
    '''
    figure, axis = plt.subplots(1, 1, dpi=150)
    axis.hist(coverage_hist_all, 100, color='red', label='all consensus sequences')
    axis.set_yscale('log')
    axis.set_ylabel('Number of peptides')
    axis.set_xlabel(f'Consensus sequence coverage')
    axis.grid()
    plt.tight_layout()
    
    
def plot_coverages(pepnets, epicore):
    '''Plot the coverages of the consensus sequences by the peptides.

    Args:
        pepnets: List containing all coverages for pepnets.
        epicore: List containing all coverages for epicore.
    '''
    figure, axis = plt.subplots(1, 1)
    axis.hist(pepnets, [i*0.01 for i in range(0,101,1)], color='blue', label='pepnets', alpha=0.5)
    axis.hist(epicore, [i*0.01 for i in range(0,101,1)], color='red', label='epicore', alpha=0.5)
    axis.set_yscale('log')
    axis.set_ylabel('Number of peptides')
    axis.set_xlabel(f'Consensus sequence coverage')
    axis.grid()
    plt.tight_layout()
    plt.legend()