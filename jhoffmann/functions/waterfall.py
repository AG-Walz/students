import pandas as pd
import ast
import itertools
import matplotlib.pyplot as plt
from matplotlib import colors
import plotly.graph_objects as go
import ast
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as mpatches
import numpy as np


def sanity_check(merged_df: pd.DataFrame):
    '''Check if dataframe containing sequence and consensus frequencies is possible.

    Args:
        merged_df: A pandas dataframe containing sequence and consensus frequencies.
    '''
    if len(merged_df[merged_df['benign_frequency_epicore']<merged_df['benign_frequency_single']]) > 0:
        print('Something went wrong!')
    if len(merged_df[merged_df['malignant_frequency_epicore']<merged_df['malignant_frequency_single']]) > 0:
        print('Something went wrong!')
    if len(merged_df[(merged_df['benign_frequency_epicore']>1)&(merged_df['benign_frequency_epicore']<0)]) > 0:
        print('Something went wrong!')
    if len(merged_df[(merged_df['malignant_frequency_epicore']>1)&(merged_df['malignant_frequency_epicore']<0)]) > 0:
        print('Something went wrong!')
    if len(merged_df[(merged_df['benign_frequency_single']>1)&(merged_df['benign_frequency_single']<0)]) > 0:
        print('Something went wrong!')
    if len(merged_df[(merged_df['malignant_frequency_single']>1)&(merged_df['malignant_frequency_single']<0)]) > 0:
        print('Something went wrong!')


def consensus_frequency(epitopes_csv: str, allotype='') -> pd.DataFrame:
    '''Calculate the frequencies of the peptide groups.

    Args:
        epitopes_csv: Path to the csv containing the epicore result with one 
            peptide group per row.
        allotype: Specifies if allotypes should be computed for a specific allotype. 

    Returns:
        A dataframe containing the frequencies of the peptide groups.
    '''
    # read in consensus epitopes
    epitopes_df = pd.read_csv(epitopes_csv, index_col=[0], usecols=['grouped_peptides_sequence','grouped_peptides_sample','grouped_peptides_condition','whole_epitopes', 'consensus_epitopes'])
    epitopes_df['grouped_peptides_sequence'] = epitopes_df['grouped_peptides_sequence'].apply(lambda x: ast.literal_eval(x))
    epitopes_df['grouped_peptides_sample'] = epitopes_df['grouped_peptides_sample'].apply(lambda x: ast.literal_eval(x))
    epitopes_df['grouped_peptides_condition'] = epitopes_df['grouped_peptides_condition'].apply(lambda x: ast.literal_eval(x))

    # summarize all peptide groups with the same whole sequence
    epitopes_df = epitopes_df.groupby(['whole_epitopes']).agg(list).reset_index()
    epitopes_df['grouped_peptides_sample'] = epitopes_df['grouped_peptides_sample'].apply(lambda x: list(itertools.chain(*x)))
    epitopes_df['grouped_peptides_condition'] = epitopes_df['grouped_peptides_condition'].apply(lambda x: list(itertools.chain(*x)))

    # compute counts
    epitopes_df['sample_condition'] = epitopes_df.apply(lambda row: list(set([sample+condition for sample, condition in zip(row['grouped_peptides_sample'],row['grouped_peptides_condition'])])), axis=1)
    if allotype == '':
        epitopes_df['benign_count'] = epitopes_df['sample_condition'].apply(lambda cell: sum('benign' in com for com in cell))
        epitopes_df['malignant_count'] = epitopes_df['sample_condition'].apply(lambda cell: len(set(com.split('~')[0] for com in cell if 'malignant' in com)))
    else:
        epitopes_df['benign_count'] = epitopes_df['sample_condition'].apply(lambda cell: sum('benign' in com for com in cell))
        epitopes_df['malignant_count'] = epitopes_df['sample_condition'].apply(lambda cell: len(set(com.split('~')[0] for com in cell if allotype in com)))
    
    # get number of total malignant and benign samples
    if allotype == '':
        malignant_total = len(set(com.split('~')[0] for com in epitopes_df['sample_condition'].explode('sample_condition').unique()  if 'malignant' in com))
        benign_total = sum('benign' in com for com in epitopes_df['sample_condition'].explode('sample_condition').unique())
    else:
        malignant_total = len(set(com.split('~')[0] for com in epitopes_df['sample_condition'].explode('sample_condition').unique() if allotype in com))
        benign_total = sum('benign' in com for com in epitopes_df['sample_condition'].explode('sample_condition').unique())
    
    # compute frequencies
    epitopes_df['benign_frequency'] = epitopes_df['benign_count'] / benign_total
    epitopes_df['malignant_frequency'] = epitopes_df['malignant_count'] / malignant_total

    epitopes_df = epitopes_df.drop(columns='sample_condition')
    return epitopes_df


def peptide_frequency(peptide_csv, allotype=''):
    ''' Calculate the frequencies of peptides.

    Args:
        peptide_csv: The Path to the peptide csv.
        allotype: Specifies if allotypes should be computed for a specific allotype. 

    Returns:
        A pandas dataframe containing the frequencies of the peptides.
    '''

    # read in data
    peptide_df = pd.read_csv(peptide_csv, index_col=[0], usecols=['sequence', 'sample', 'condition', 'accessions'])

    # compute counts
    peptide_df = peptide_df.groupby(['sequence']).agg(list).reset_index()
    peptide_df['sample_condition'] = peptide_df.apply(lambda row: list(set([sample+cond for sample, cond in zip(row['sample'],row['condition'])])), axis=1)
    peptide_df['benign_count'] = peptide_df['sample_condition'].apply(lambda cell: sum('benign' in com for com in cell))
    benign_total = sum('benign' in com for com in peptide_df['sample_condition'].explode('sample_condition').unique())

    # get number of total malignant and benign samples
    if allotype == '':
        peptide_df['malignant_count'] = peptide_df['sample_condition'].apply(lambda cell: len(set(com.split('~')[0] for com in cell if 'malignant' in com)))
        malignant_total = len(set(com.split('~')[0] for com in peptide_df['sample_condition'].explode('sample_condition').unique() if 'malignant' in com))
    else:
        peptide_df['malignant_count'] = peptide_df['sample_condition'].apply(lambda cell: len(set(com.split('~')[0] for com in cell if allotype in com)))
        malignant_total = len(set(com.split('~')[0] for com in peptide_df['sample_condition'].explode('sample_condition').unique() if allotype in com))

    # compute frequencies
    peptide_df['benign_frequency'] = peptide_df['benign_count'] / benign_total
    peptide_df['malignant_frequency'] = peptide_df['malignant_count'] / malignant_total

    peptide_df = peptide_df.drop(columns='sample_condition')
    
    return peptide_df


def plot_waterfall(df: pd.DataFrame, column: str, xlabel:str, axis):
    '''Plot the waterfall plot to the data stored in df.

    Args: 
        df: A dataframe containing the benign and malignant frequencies of peptide
            groups and/or peptide sequences.
        column: The header of the column containing the peptide or consensus sequences.
        xlabel: The string that is put as a xlabel on the plot.
    '''
    # group peptide (groups) that have the same frequencies
    binned_df = df.drop_duplicates([column]).groupby(['benign_frequency','malignant_frequency']).agg(list).reset_index()
    binned_df['count'] = binned_df[column].str.len()

    # compute waterfall plot order
    binned_df['ratio'] = (binned_df['benign_frequency'] + 0.000001) / (binned_df['malignant_frequency'] + 0.000001)
    binned_df = binned_df.sort_values('ratio')

    # select only peptide (groups) that occur in malignant samples
    binned_df = binned_df[binned_df['malignant_frequency']!=0]
    axis.bar(binned_df['count'].cumsum()-binned_df['count'], binned_df['malignant_frequency'],binned_df['count'], align='edge', color='red', label='malignant')
    axis.bar(binned_df['count'].cumsum()-binned_df['count'], -binned_df['benign_frequency'],binned_df['count'], align='edge', color='blue', label='benign')
    axis.legend()
    n_groups = max(binned_df['count'].cumsum())
    axis.set_xticks([i for i in range(n_groups) if i%6000==0])
    axis.set_yticks([-0.6,-0.3,0,0.3,0.6,0.9],[0.6,0.3,0,0.3,0.6,0.9])
    axis.set_xlabel(xlabel)
    axis.set_ylabel('Frequency')


def candidate_region_epicore(df, axis, malignant_threshold, benign_threshold, original=False):
    '''Plot the candidate region of a dataset on consensus sequence level.

    Args:
        df: A dataframe containing the peptide level and consensus level frequencies.
        axis: The axis of a matplotlib plot.
        malignant_threshold: The malignant cutoff value for the candidate region.
        benign_threshold: The benign cutoff value for the candidate region.
        original: Boolean indicating if the order of the candidates on peptide level should be kept
    '''

    df = df.groupby(['whole_epitopes']).agg({'malignant_frequency_epicore':'first', 'benign_frequency_epicore':'first', 'malignant_frequency_single':'max', 'benign_frequency_single':'max'})
    df['ratio'] = (df['benign_frequency_epicore'] + 0.000001) / (df['malignant_frequency_epicore'] + 0.000001)
    df = df.sort_values('ratio')
    df = df[(df['benign_frequency_epicore']<=benign_threshold)&(df['malignant_frequency_epicore']>=malignant_threshold)]
    if original:
        df = df.sort_values(['malignant_frequency_single', 'benign_frequency_single'], ascending=[False, True])
    axis.bar([i for i in range(len(df))], df['malignant_frequency_epicore'], color='red', label='epicore malignant')
    axis.bar([i for i in range(len(df))], df['malignant_frequency_single'], color='white')
    axis.bar([i for i in range(len(df))], df['malignant_frequency_single'], color='red', alpha=0.4, label='peptide malignant')
    if max(df['benign_frequency_epicore']) > 0:
        axis.bar([i for i in range(len(df))], -df['benign_frequency_epicore'], color='black', label='epicore benign')
    if max(df['benign_frequency_single']) > 0:
        axis.bar([i for i in range(len(df))], -df['benign_frequency_single'], color='white')
        axis.bar([i for i in range(len(df))], -df['benign_frequency_single'], color='black', alpha=0.4, label='peptide benign')
    axis.legend(fontsize=10)
    axis.set_xticks([i for i in range(0,len(df), 5)])
    axis.set_xlabel('Peptide group', fontsize=10)
    axis.set_ylabel('Frequency', fontsize=10)


def candidate_region_epicore_highlighted(df, axis, malignant_threshold, benign_threshold, candidate_groups):
    '''Plot the candidate region of a dataset on consensus sequence level with some candidates being highlighted.

    Args:
        df: A dataframe containing the peptide level and consensus level frequencies.
        axis: The axis of a matplotlib plot.
        malignant_threshold: The malignant cutoff value for the candidate region.
        benign_threshold: The benign cutoff value for the candidate region.
        original: Boolean indicating if the order of the candidates on peptide level should be kept
    '''

    df = df.groupby(['whole_epitopes']).agg({'malignant_frequency_epicore':'first', 'benign_frequency_epicore':'first', 'malignant_frequency_single':'max', 'benign_frequency_single':'max'}).reset_index()
    df['ratio'] = (df['benign_frequency_epicore'] + 0.000001) / (df['malignant_frequency_epicore'] + 0.000001)
    df = df.sort_values('ratio')
    df = df[(df['benign_frequency_epicore']<=benign_threshold)&(df['malignant_frequency_epicore']>=malignant_threshold)]
    axis.bar([i for i in range(len(df))], df['malignant_frequency_epicore'], color='red', label='epicore malignant')
    axis.bar([i for i in range(len(df))], df['malignant_frequency_single'], color='white')
    axis.bar([i for i in range(len(df))], df['malignant_frequency_single'], color='red', alpha=0.4, label='peptide malignant')
    if max(df['benign_frequency_epicore']) > 0:
        axis.bar([i for i in range(len(df))], -df['benign_frequency_epicore'], color='black', label='epicore benign')
    if max(df['benign_frequency_single']) > 0:
        axis.bar([i for i in range(len(df))], -df['benign_frequency_single'], color='white')
        axis.bar([i for i in range(len(df))], -df['benign_frequency_single'], color='black', alpha=0.4, label='peptide benign')

    axis.bar([i for i in range(len(df))], [row['malignant_frequency_epicore'] if row['whole_epitopes'] in candidate_groups else 0 for i, row in df.iterrows()], color='green')
    axis.bar([i for i in range(len(df))], [row['malignant_frequency_single'] if row['whole_epitopes'] in candidate_groups else 0 for i, row in df.iterrows()], color='white')
    axis.bar([i for i in range(len(df))], [row['malignant_frequency_single'] if row['whole_epitopes'] in candidate_groups else 0 for i, row in df.iterrows()], color='green', alpha=0.4, label='Warehouse peptide')
    
    axis.legend(fontsize=10)
    axis.set_xticks([i for i in range(0,len(df), 5)])
    axis.set_xlabel('Peptide group', fontsize=10)
    axis.set_ylabel('Frequency', fontsize=10)


def candidate_region_peptide(df, axis, malignant_threshold, benign_threshold, lv_threshold):
    '''Plot the candidate region of a dataset on peptide level.

    Args: 
        df: A dataframe containing the peptide level and consensus level frequencies.
        axis: The axis of a matplotlib plot.
        malignant_threshold: The malignant cutoff value for the candidate region.
        benign_threshold: The benign cutoff value for the candidate region.
        lv_threshold: The threshold value for peptide groups to be highlighted as 
            peptide groups including benign length variants. 
    '''

    # get candidate region
    df = df.groupby(['sequence']).agg({'malignant_frequency_epicore':'first', 'benign_frequency_epicore':'first', 'malignant_frequency_single':'first','benign_frequency_single':'first'})
    df['ratio'] = (df['benign_frequency_single'] + 0.000001) / (df['malignant_frequency_single'] + 0.000001)
    df = df.sort_values('ratio')
    df = df[(df['benign_frequency_single']<=benign_threshold)&(df['malignant_frequency_single']>=malignant_threshold)]

    # plot candidate region
    barcolors = [colors.to_rgb('red')[0:3]+(1,) if b <= lv_threshold else colors.to_rgb('red')[0:3]+(0.2,) for b in df['benign_frequency_epicore']]
    labels = ['_no benign lv' if b == 0 else '_benign lv' for b in df['benign_frequency_epicore']]
    labels[labels.index('_no benign lv')] = 'No benign lv'
    labels[labels.index('_benign lv')] = 'Benign lv'
    axis.bar([i for i in range(len(df))], df['malignant_frequency_single'], color=barcolors, label=labels)
    axis.legend(fontsize=15)
    axis.set_xlabel('Peptide', fontsize=15)


def sequence_flow(frequency_df, benign_frequency_sequence, benign_frequency_consensus, malignant_frequency_sequence, malignant_frequency_consensus):
    '''Calculate the sequence flows between categories.

    Args:
        frequency_df: A pandas DataFrame containing in each row one sequence.
        benign_frequency_sequence: The header of the column containing the 
            frequency of a sequence in benign samples.
        benign_frequency_consensus:The header of the column containing the 
            frequency of a consensus sequence the peptide sequence of the row 
            contributes to in benign samples.
        malignant_frequency_sequence: The header of the column containing the 
            frequency of a sequence in malignant samples.
        malignant_frequency_consensus: The header of the column containing the 
            frequency of a consensus sequence the peptide sequence of the row 
            contributes to in malignant samples.
    
    Returns:
        Returns a list containing the flows between categories on sequence and consensus level.
    '''
    candidate_candidate = len(frequency_df[(frequency_df[malignant_frequency_sequence] >= 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] >= 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))
    candidate_benign = len(frequency_df[(frequency_df[malignant_frequency_sequence] >= 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] == 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    candidate_both = len(frequency_df[(frequency_df[malignant_frequency_sequence] >= 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] > 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    candidate_malignant = len(frequency_df[(frequency_df[malignant_frequency_sequence] >= 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] < 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))

    malignant_candidate = len(frequency_df[(frequency_df[malignant_frequency_sequence] < 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] >= 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))
    malignant_benign = len(frequency_df[(frequency_df[malignant_frequency_sequence] < 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] == 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    malignant_both = len(frequency_df[(frequency_df[malignant_frequency_sequence] < 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] > 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    malignant_malignant = len(frequency_df[(frequency_df[malignant_frequency_sequence] < 0.2) & (frequency_df[benign_frequency_sequence] == 0) & (frequency_df[malignant_frequency_consensus] < 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))

    benign_candidate = len(frequency_df[(frequency_df[malignant_frequency_sequence] == 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] >= 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))
    benign_benign = len(frequency_df[(frequency_df[malignant_frequency_sequence] == 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] == 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    benign_both = len(frequency_df[(frequency_df[malignant_frequency_sequence] == 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] > 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    benign_malignant = len(frequency_df[(frequency_df[malignant_frequency_sequence] == 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] < 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))

    both_candidate = len(frequency_df[(frequency_df[malignant_frequency_sequence] > 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] >= 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))
    both_benign = len(frequency_df[(frequency_df[malignant_frequency_sequence] > 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] == 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    both_both = len(frequency_df[(frequency_df[malignant_frequency_sequence] > 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] > 0) & (frequency_df[benign_frequency_consensus] > 0)].drop_duplicates('grouped_peptides_sequence'))
    both_malignant = len(frequency_df[(frequency_df[malignant_frequency_sequence] > 0) & (frequency_df[benign_frequency_sequence] > 0) & (frequency_df[malignant_frequency_consensus] < 0.2) & (frequency_df[benign_frequency_consensus] == 0)].drop_duplicates('grouped_peptides_sequence'))

    return candidate_candidate, candidate_malignant, candidate_benign, candidate_both, malignant_candidate, malignant_malignant, malignant_benign, malignant_both, benign_candidate, benign_malignant, benign_benign, benign_both, both_candidate, both_malignant, both_benign, both_both



def plot_flow_all(flows):
    '''Plot the flow of sequences between categories from sequence to consensus sequence level.

    Args:
        flow: The sequence flows.
    '''
    fig = go.Figure(data=[go.Sankey(
        valueformat = '.0f',
        valuesuffix = 'TWh',
        node = dict(
            pad = 15,
            thickness = 20,
            line = dict(color = "black", width = 0.5),
            label = ['Candidate', 'Malignant <0.2', 'Benign', 'Both','Candidate_epicore', 'Malignant_epicore<0.2', 'Benign_epicore', 'Both_epicore'],
            color = 'blue'
        ),
        link = dict(
            source = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3],
            target = [4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7],
            value =  flows,
            color = ['red']*5+['grey']*3+['red'] + ['grey']*3+['red'] + ['grey']*3
        )
    )])
    fig.show()


def mask_peptide(sequence, core_start, core_end, label, core_pos):
    '''Mask peptides.

    Args: 
        sequence: The peptide sequence to be masked.
        core_start: The start position of the subsequence being masked.
        core_end: The end position of the subsequence being masked.
        label: The label inserted in the masked sequence
    
    Returns: 
        The input sequence with a masked core and the specified label.
    '''
    masked_sequence = sequence[0:core_start]
    core_length = core_end-core_start+1
    if core_length != 0:
        masked_sequence += '-'*core_pos
        masked_sequence += label
        masked_sequence += '-'*(core_length-len(label)-core_pos)
        masked_sequence += sequence[core_end+1:]
    else:
        masked_sequence = sequence
    return masked_sequence



def candidate_heatmap(epitopes_csv: str, candidate: str, ax: plt.axes, core_pos, label, legend_anchor, allotype='', allotype_samples=None, n_benign=None, vmax=0.2):
    '''Create a heatmap indicating which length variant is present in which sample.

    Args:
        epitopes_csv: The path to the file containing by epicore computed peptide groups. 
        candidate: The sequence of the peptide group.
        ax: Matplotlib axis for the heatmap.    
    '''

    # load peptide groups
    df = pd.read_csv(epitopes_csv, index_col=[0])
    group = df[(df['whole_epitopes']==candidate)].iloc[0]
    group_df = pd.DataFrame({'samples':ast.literal_eval(group['grouped_peptides_sample'])
                            ,'sequences':ast.literal_eval(group['grouped_peptides_sequence'])
                            ,'conditions':ast.literal_eval(group['grouped_peptides_condition'])
                            ,'start':ast.literal_eval(group['grouped_peptides_start'])
                            ,'end':ast.literal_eval(group['grouped_peptides_end'])
                            ,'consensus_start': group['core_epitopes_start']
                            ,'consensus_end': group['core_epitopes_end']})


    # mask sequence
    group_df['consensus_start'] = group_df.apply(lambda cell: max(0,cell['consensus_start']-cell['start']), axis=1)
    group_df['consensus_end'] = group_df.apply(lambda cell: min(cell['consensus_end']-cell['start'], len(cell['sequences'])-1), axis=1)
    group_df['end'] = group_df['end'].apply(lambda cell: cell - group_df['start'].min())
    group_df['start'] = group_df['start'].apply(lambda cell: cell - group_df['start'].min())
    group_df = group_df.groupby('sequences').agg({'samples':list,'conditions':list, 'start': 'first', 'end': 'first', 'consensus_start':'first', 'consensus_end':'first'}).reset_index()
    group_df = group_df.sort_values(['start', 'end'], ascending=[True, False])
    group_df['core_pos'] = core_pos
    group_df = group_df.explode(['samples', 'conditions'])
    group_df['consensus_end'] = group_df.apply(lambda row: max(min(row['consensus_end'], row['end']),-1), axis=1)
    group_df['sequences'] = group_df.apply(lambda row: mask_peptide(row['sequences'], row['consensus_start'], row['consensus_end'], label, row['core_pos']), axis=1)

    malignant_df = group_df[group_df['conditions'].str.contains('malignant')]
    malignant_df = malignant_df[malignant_df['conditions'].str.contains(allotype, regex=False)]
    benign_df = group_df[group_df['conditions'].str.contains('benign')].drop_duplicates(['sequences', 'samples'])

    # build heatmap matrix
    samples = malignant_df['samples'].to_list()
    sequences = malignant_df['sequences'].to_list()
    c = [f'{sample}_{sequence}' for sample, sequence in zip(samples, sequences)]
    unique_sequences = group_df['sequences'].unique()
    matrix = np.zeros((len(unique_sequences), len(allotype_samples)+1))


    allotype_samples = sorted(allotype_samples, key=lambda sample: int(sample[-2:]))
    # add malignant occurrences
    for s_i, sample in enumerate(allotype_samples):
        for se_i, sequence in enumerate(unique_sequences):
            if f'{sample}_{sequence}' in c:
                matrix[se_i, s_i] = 1
    # compute benign frequencies
    for i, sample in enumerate(allotype_samples):
        for j, sequence in enumerate(unique_sequences):
            if matrix[j,i] == 0:
                matrix[j,i] = - len(benign_df[benign_df['sequences']==sequence])/n_benign

    # identify peptide covering the most benign samples and highlight it in the heatmap
    predominant_lv = matrix.sum(axis=1).argmax()
    matrix[predominant_lv][np.where(matrix[predominant_lv]==0)] = 0.5

    # create heatmap
    cmap=ListedColormap(["#BF1701",'none',"#b1d1fc7f",'#0a437a'])

    # ensure correct mappings for colors
    norm = [-1.1,-0.2,0.2,0.7,1.1]
    norm = BoundaryNorm(norm, 4)
    # split data for two color scales
    data1 = np.ma.masked_array(matrix, (matrix<=0))
    data2 = np.ma.masked_array(matrix, (matrix>0))
    im1 = ax.imshow(data1, cmap=cmap, norm=norm)
    cmap = colors.LinearSegmentedColormap.from_list('',['white', "#810404"])
    im2 = ax.imshow(np.abs(data2), cmap=cmap, vmin=0, vmax=vmax)

    # align peptide sequences for better overview
    sequences_aligned = []
    for sequence in unique_sequences:
        sequences_aligned.append(sequence+' '*(group_df['end'].max()-group_df[group_df['sequences']==sequence].iloc[0]['end']))

    # add legend and ticks
    ax.set_xticks([i for i in range(len(allotype_samples))], allotype_samples, rotation=90)
    ax.set_yticks([i for i in range(len(unique_sequences))])
    ax.set_yticklabels(sequences_aligned, fontfamily='monospace')
    legend_elements = []
    if 1 in matrix: 
        legend_elements.append(mpatches.Patch(color="#0455a1",label='Malignant occurrence'))
    if -1 in matrix: 
        legend_elements.append(mpatches.Patch(color="#BF1701",label='Benign occurrence'))
    if 0.5 in matrix: 
        legend_elements.append(mpatches.Patch(color="#b1d1fc7f",label='Predominant length variant'))

    return im2, legend_elements




def response_heatmap(epitopes_csv: str, candidate: str, ax: plt.axes, core_pos, label, allotype='', allotype_samples=None, n_benign=None, vmax=0.2, responses=None):
    '''Create a heatmap indicating T-cell responses to specific length variants.

    Args:
        epitopes_csv: The path to the file containing by epicore computed peptide groups. 
        candidate: The sequence of the peptide group.
        ax: Matplotlib axis for the heatmap.    
    '''

    # load peptide groups
    df = pd.read_csv(epitopes_csv, index_col=[0])
    group = df[(df['grouped_peptides_sequence'].str.contains(candidate))].iloc[0]
    group_df = pd.DataFrame({'samples':ast.literal_eval(group['grouped_peptides_sample'])
                            ,'sequences':ast.literal_eval(group['grouped_peptides_sequence'])
                            ,'conditions':ast.literal_eval(group['grouped_peptides_condition'])
                            ,'start':ast.literal_eval(group['grouped_peptides_start'])
                            ,'end':ast.literal_eval(group['grouped_peptides_end'])
                            ,'consensus_start': group['core_epitopes_start']
                            ,'consensus_end': group['core_epitopes_end']})


    # mask sequence
    group_df['consensus_start'] = group_df.apply(lambda cell: max(0,cell['consensus_start']-cell['start']), axis=1)
    group_df['consensus_end'] = group_df.apply(lambda cell: min(cell['consensus_end']-cell['start'], len(cell['sequences'])-1), axis=1)
    group_df['end'] = group_df['end'].apply(lambda cell: cell - group_df['start'].min())
    group_df['start'] = group_df['start'].apply(lambda cell: cell - group_df['start'].min())
    group_df = group_df.groupby('sequences').agg({'samples':list,'conditions':list, 'start': 'first', 'end': 'first', 'consensus_start':'first', 'consensus_end':'first'}).reset_index()
    group_df = group_df.sort_values(['start', 'end'], ascending=[True, False])
    group_df['core_pos'] = core_pos
    group_df = group_df.explode(['samples', 'conditions'])
    group_df['consensus_end'] = group_df.apply(lambda row: max(min(row['consensus_end'], row['end']),-1), axis=1)
    group_df['sequences'] = group_df.apply(lambda row: mask_peptide(row['sequences'], row['consensus_start'], row['consensus_end'], label, row['core_pos']), axis=1)

    malignant_df = group_df[group_df['conditions'].str.contains('malignant')]
    malignant_df = malignant_df[malignant_df['conditions'].str.contains(allotype, regex=False)]
    benign_df = group_df[group_df['conditions'].str.contains('benign')].drop_duplicates(['sequences', 'samples'])

    # build heatmap matrix
    samples = malignant_df['samples'].to_list()
    sequences = malignant_df['sequences'].to_list()
    c = [f'{sample}_{sequence}' for sample, sequence in zip(samples, sequences)]
    unique_sequences = group_df['sequences'].unique()
    matrix = np.zeros((len(unique_sequences), len(allotype_samples)+1))


    allotype_samples = sorted(allotype_samples, key=lambda sample: int(sample[-2:]))
    for s_i, sample in enumerate(allotype_samples):
        for se_i, sequence in enumerate(unique_sequences):
            if f'{sample}_{sequence}' in c:
                if sample in responses:
                    matrix[se_i, s_i] = 2
                else:
                    matrix[se_i, s_i] = 1


    predominant_lv = matrix.sum(axis=1).argmax()
    matrix[predominant_lv][np.where(matrix[predominant_lv]==0)] = 0.5

    # plot heatmap
    cmap=ListedColormap(["#BF1701",'none',"#b1d1fc7f",'#0a437a', "#053105"])
    norm = [-1.1,-0.2,0.2,0.7,1.1,2.1]
    norm = BoundaryNorm(norm, 5)
    data1 = np.ma.masked_array(matrix, (matrix<=0))
    data2 = np.ma.masked_array(matrix, (matrix>0))
    im1 = ax.imshow(data1, cmap=cmap, norm=norm)#, aspect ='auto')
    cmap = colors.LinearSegmentedColormap.from_list('',['white', "#810404"])
    im2 = ax.imshow(np.abs(data2), cmap=cmap, vmin=0, vmax=vmax)

    # align peptide sequences for better overview
    sequences_aligned = []
    for sequence in unique_sequences:
        sequences_aligned.append(sequence+' '*(group_df['end'].max()-group_df[group_df['sequences']==sequence].iloc[0]['end']))

    # add legend and ticks
    ax.set_xticks([i for i in range(len(allotype_samples))], allotype_samples, rotation=90)
    ax.set_yticks([i for i in range(len(unique_sequences))])
    ax.set_yticklabels(sequences_aligned, fontfamily='monospace')
    legend_elements = []
    if 1 in matrix: 
        legend_elements.append(mpatches.Patch(color="#0455a1",label='Peptide present in patient without response'))
    if -1 in matrix: 
        legend_elements.append(mpatches.Patch(color="#BF1701",label='Benign occurrence'))
    if 0.5 in matrix: 
        legend_elements.append(mpatches.Patch(color="#b1d1fc7f",label='Predominant length variant'))
    if 2 in matrix: 
        legend_elements.append(mpatches.Patch(color="#053105",label='Peptide present in patient with response'))

    return im2, legend_elements