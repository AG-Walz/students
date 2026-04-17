import matplotlib.pyplot as plt


def length_distribution(class1_df, class2_df):
    '''Plot the length distribution of the dataset.

    Args: 
        class1_df: Pandas dataframe containing the class1 peptides.
        class2_df: Pandas dataframe containing the class2 peptides.
    '''

    figure, axis = plt.subplots(1,2)

    # classI
    class1_l = class1_df.apply(lambda row: int(row['end'].split(';')[0])+1-int(row['start'].split(';')[0]), axis=1).to_list()
    class1_shortest = min(class1_l)
    class1_longest = max(class1_l)
    axis[0].hist(class1_l, bins=[i for i in range(class1_shortest,class1_longest +2,1)])
    axis[0].set_xlabel('peptide_length')
    axis[0].set_ylabel('number of peptides')
    axis[0].set_title('class1 peptides')


    # classII
    class2_l = class2_df.apply(lambda row: int(row['end'].split(';')[0])+1-int(row['start'].split(';')[0]), axis=1).to_list()
    class2_shortest = min(class2_l)
    class2_longest = max(class2_l)
    axis[1].hist(class2_l, bins=[i for i in range(class2_shortest,class2_longest +2,1)])
    axis[1].set_xlabel('peptide_length')
    axis[1].set_ylabel('number of peptides')
    axis[1].set_title('class2 peptides')

    plt.tight_layout()
    plt.show()
