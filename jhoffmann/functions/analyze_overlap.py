''''
This script contains the overlap computations performed by epicore when --QC is set.
'''

import pandas as pd



def included(
    start: int, end: int, overlap_df: pd.DataFrame, included_df: pd.DataFrame
) -> bool:
    """Checks of a peptide is included in all peptides of the peptide group.

    Args:
        start: The start position of the peptide.
        end: The end position of the peptide.
        overlap_df: Dataframe containing information about peptide overlap.
        included_df: Dataframe containing information about peptide inclusions.

    Returns:
        Boolean that indicates if a peptide overlaps that overlaps with all peptides is completely included in any of the peptides.
    """

    # get peptides that overlap with all peptides in the group
    included_df = included_df[overlap_df[overlap_df.all(axis=1)].index]

    # get row of interest and check if any value is true
    return included_df.loc[f"{start}-{end}"].any()


def get_largest_overlap(
    starts: list[int], ends: list[int], pep_start: int, pep_end: int, intern: bool
) -> int:
    """Computes the maximal overlap of a peptide with the peptides of a group.

    Args:
        starts: The start positions of the peptides in a group.
        ends: The end positions of the peptides in a group.
        pep_start: The start position of a peptide.
        pep_end: The end position of a peptide.
        intern: Boolean indicating if overlap is computed within a group.

    Returns:
        The maximal overlap between the peptide and any peptide of the group.
    """
    max_overlap = 0

    # iterate over all peptides in the peptide group.
    for start, end in zip(starts, ends):
        if (start == pep_start) & (end == pep_end):  # peptide occurs in the group
            if intern:
                continue
            else:
                return 0

        if start <= pep_start:  # peptide occurs before the current
            max_overlap = max(
                min(pep_end - pep_start + 1, end - pep_start + 1), max_overlap
            )

        else:  # peptide occurs after the current
            max_overlap = max(min(pep_end - start + 1, end - start + 1), max_overlap)

    return min(max_overlap, pep_end - pep_start + 1)


def get_minimal_overlap(
    starts: list[int], ends: list[int], pep_start: int, pep_end: int
) -> int:
    """Computes the minimal overlap of a peptide with the peptides of a group.

    Args:
        starts: The start positions of the peptide in one group.
        ends: The end positions of the peptide in one group.
        pep_start: The start position of a peptide.
        pep_end: The end position of a peptide.

    Returns:
        The minimal overlap between the peptide and any peptide of the group.
    """

    if len(starts) == 1:
        return 0

    min_overlap = 100
    # iterate over all peptides in the peptide group.
    for start, end in zip(starts, ends):

        if (start == pep_start) & (end == pep_end):  # peptide occurs in the group
            continue

        if (start <= pep_start) & (end <= pep_end):  # peptide occurs before the current
            min_overlap = max(0, min(end - pep_start + 1, min_overlap))

        elif (start > pep_start) & (end > pep_end):  # peptide occurs after the current
            min_overlap = max(0, min(pep_end - start + 1, min_overlap))

        else:
            min_overlap = min(end - start + 1, min_overlap)

    return min(min_overlap, pep_end - pep_start + 1)


def consensus_coverage(
    start: int, end: int, consensus_start: int, consensus_end: int
) -> int:
    """Computes the consensus sequence coverage by a peptide.

    Args:
        start: The start position of a peptide.
        end: The end position of a peptide.
        consensus_start: The start of the consensus sequence.
        consensus_end: The end of the consensus sequence.

    Returns:
        The consensus sequence coverage.
    """
    consensus_len = consensus_end - consensus_start + 1
    if start <= consensus_start:
        return (
            max(0, min(end - consensus_start + 1, consensus_end - consensus_start + 1))
            / consensus_len
        )
    else:
        return max(0, min(consensus_end - start + 1, end - start + 1)) / consensus_len


def build_overlap_df(
    starts: list[int], ends: list[int]
) -> tuple([pd.DataFrame, pd.DataFrame]):
    """Creates two dataframes tracking peptides overlap and inclusion.

    Args:
        starts: The start positions of the peptides in a peptide group.
        ends: The end positions of the peptides in a peptide group.

    Returns:
        Returns a tuple containing two dataframes. The first dataframe contains
        information about if two peptides overlap with each other. The second
        dataframe contains information about if one peptide is completely
        included in the other.
    """
    # create df tracking if peptide is included or overlaps
    column_names = [f"{start}-{end}" for start, end in zip(starts, ends)]
    overlap_df = pd.DataFrame(columns=column_names)
    included_df = pd.DataFrame(columns=column_names)

    # iterate all peptide pairs of the peptide group
    for i, (start_1, end_1) in enumerate(zip(starts, ends)):
        for start_2, end_2 in zip(starts[i + 1 :], ends[i + 1 :]):

            # check if peptides overlap
            overlap = min(end_2, end_1) - start_2 + 1
            if overlap > 0:
                overlap_bool = True
            else:
                overlap_bool = False
            overlap_df.loc[f"{start_2}-{end_2}", f"{start_1}-{end_1}"] = overlap_bool
            overlap_df.loc[f"{start_1}-{end_1}", f"{start_2}-{end_2}"] = overlap_bool

            # check if peptide is included in other peptide
            included_df.loc[f"{start_2}-{end_2}", f"{start_1}-{end_1}"] = end_2 <= end_1

    included_df.loc[f"{starts[0]}-{ends[0]}", f"{starts[0]}-{ends[0]}"] = None
    return overlap_df, included_df


def all_overlap(row: pd.Series) -> pd.DataFrame:
    """Compute sequence overlap between peptide and peptide of previous and next group.

    Args:
        row: Pandas series containing the peptide groups of one protein.

    Returns:
        A pandas dataframe containing information about the overlap of peptides
        within and outside of a peptide group.
    """

    # define dataframe containing the overlap to the previous, current and next peptide group
    pep = pd.DataFrame(
        columns=[
            "previous",
            "intern_max",
            "intern_min",
            "next",
            "len_pep",
            "sequence",
            "included",
        ]
    )
    k = 0

    # iterate all peptide groups
    for i, (group_starts,group_ends,group_sequences) in enumerate(
        zip(row["grouped_peptides_start"],row["grouped_peptides_end"],row["grouped_peptides_sequence"])
    ):
        pos_df = pd.DataFrame({"starts": group_starts, "ends": group_ends, "sequences": group_sequences})
        pos_df = pos_df.drop_duplicates()
        group_starts = pos_df["starts"].to_list()
        group_ends = pos_df["ends"].to_list()
        group_sequences = pos_df["sequences"].to_list()
        overlap_df, included_df = build_overlap_df(group_starts, group_ends)

        # iterate all peptide of the current group
        for _, (start, end, sequence) in enumerate(
            zip(group_starts, group_ends, group_sequences)
        ):

            # calculate maximal overlap to previous group
            if i > 0:
                pep.loc[k, "previous"] = get_largest_overlap(
                    row["grouped_peptides_start"][i - 1],
                    row["grouped_peptides_end"][i - 1],
                    start,
                    end,
                    False,
                )
            else:
                pep.loc[k, "previous"] = 0

            # calculate maximal overlap to next group
            if i < len(row["grouped_peptides_start"]) - 1:
                pep.loc[k, "next"] = get_largest_overlap(
                    row["grouped_peptides_start"][i + 1],
                    row["grouped_peptides_end"][i + 1],
                    start,
                    end,
                    False,
                )
            else:
                pep.loc[k, "next"] = 0

            # calculate the maximal and minimal overlap within a group
            pep.loc[k, "intern_max"] = get_largest_overlap(
                group_starts, group_ends, start, end, True
            )
            pep.loc[k, "intern_min"] = get_minimal_overlap(
                group_starts, group_ends, start, end
            )
            pep.loc[k, "len_pep"] = end - start + 1
            pep.loc[k, "sequence"] = sequence
            pep.loc[k, "included"] = included(start, end, overlap_df, included_df)
            k += 1

    return pep