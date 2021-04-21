import pandas as pd
from tqdm import tqdm


def merge_columns(df, mapping_dict):
    '''Function that merges columns by taking the mean of them.
    Merges based on if they have the same element in the meta_samples meta_column.
    Returns:
        new_df: Dataframe but with the columns of interest merged. Also column names are now based
                on the unique meta_samples[meta_column].
    '''

    vals = mapping_dict.keys()#np.unique(meta_samples[meta_column].values)
    new_df = pd.DataFrame(index=df.index, columns=vals)
    for i in vals:
        if not mapping_dict[i] == []:#(meta_samples[meta_column] == i).any():
            new_col = (df.loc[:, mapping_dict[i]])
            new_col = new_col.mean(axis=1)
            new_df.loc[:, i] = new_col
    new_df = new_df.loc[:, ~(new_df.isnull().all())]
    return new_df


def load_merged_peaks():
    return


def check_strand(x, anno): #anno_df_otherrna
    if x["Parent"] in anno.index:
        return x["Strand"] == anno.loc[x["Parent"],"Strand"]
    else:
        return False


def load_homer_peaks():
    return


def collect_stats():
    stats = pd.DataFrame(index=["ncRNA"],
                         columns=["NCBI annotation", "Total",
                                  "Mean Distance", "Median Distance",
                                  "Mode Distance",
                                  "Standard Deviation Distance",
                                  "Mean absolute Distance",
                                  "Median absolute Distance",
                                  "Mode absolute Distance",
                                  "Standard Deviation absolute Distance",
                                  "Number of TSS in concordance",
                                  "Number of TSS <= 10nt",
                                  "Number of TSS <= 150nt",
                                  "Number of TSS > 150nt",
                                  "Number of TSS > 10nt"])
    return


def expand_anno_id(df, break_char="=",colname=8,cpus=32):
    #from pandarallel import pandarallel
    #pandarallel.initialize(nb_workers=cpus)
    df = df.copy()

    for ind, val in tqdm(df.iterrows()):
        curr = val[colname].split(';')
        for i in curr:
            v = i.strip().replace('"', "")
            if len(v) == 0:
                continue
            curr_split = v.split(break_char)
            df.at[ind, curr_split[0]] = curr_split[1]

    return df
