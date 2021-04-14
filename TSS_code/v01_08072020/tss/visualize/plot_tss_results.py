import pandas as pd
import numpy as np
from cycler import cycler
from collections import Counter
from matplotlib import pyplot as plt
import matplotlib as mpl
import pickle
import os
from matplotlib import rcParams
import pandas as pd
from os.path import dirname
import seaborn as sns
from matplotlib_venn import venn2
from tss.visualize.fig_utils import helper_save
#rcParams['figure.figsize'] = 8, 6
#mpl.style.use('fivethirtyeight')
#mpl.style.use('ggplot')
#mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
from mplh.color_utils import get_colors
from mplh.fig_utils import legend_from_color, helper_save
from collections import defaultdict
from os.path import join

############################################################
### Plot the peak results
### 1. Peaks per each annotated tss
### 2. For each sample, how many annotated tss were covered
### 3. For each tissue, how many annotated tss were covered
### 4. For each peak, what fraction of tissues/samples covered it
############################################################
def plot_peaks_per_landmark(f_in,landmark_name,f_save=None):
    f = plt.figure(dpi=300)
    ax = f.add_subplot(1, 1, 1)
    # ax.hist(gene_df['Number of SS']);

    df = pickle.load(open(f_in, 'rb'))
    df_small = df[df['Number of SS'] <= 10]
    df_small = df_small.groupby(['Number of SS']).count()['hasGene']
    print(df_small.head())
    df_small['>10'] = np.sum(df['Number of SS'] > 10)
    df_small.plot.bar(color='b')
    ax.set_xlabel('Number of peaks per ' + landmark_name)
    ax.set_ylabel('Counts')
    # ax.set_xticks(range(0,max(gene_df['Number of SS'])+1,2))
    print('Number of genes with greater than 10 peaks:', np.sum(df['Number of SS'] > 10))
    print('Percent identified: ', 1.0 * np.sum(df['hasGene']) / df.shape[0])
    plt.title('Peaks per ' + landmark_name)
    plt.tight_layout()
    helper_save(f_save)
    return


def fraction_of_samples_with_peak(df,sample_expression, f_save, landmark_name, meta_samples=None, sample_name = None):
    f = plt.figure(dpi=300)
    ax = (np.sum((sample_expression > 0)) / df.shape[0]).plot.bar(color='b')
    ax.set_title('Fraction of ' + landmark_name + ' with a peak')
    names = map(lambda x: os.path.basename(x),sample_expression)
    ax.set_xticklabels(names)
    plt.tight_layout()
    helper_save(f_save)
    return


###
# Looks at the landmark-matrix and plots fraction of landmarks with a peak
###
def plot_sample_by_genes(f_mat_in,f_in,landmark_name, f_save=None,to_close = False):
    gene_tissue_matrix = pd.read_csv(f_mat_in,sep='\t',index_col=0)
    df = pickle.load(open(f_in,'rb'))
    f = plt.figure()
    ax = (np.sum((gene_tissue_matrix > 0)) / df.shape[0]).plot.bar()
    ax.set_title('Fraction of ' + landmark_name + ' with a peak')
    ax.set_xticklabels(list(map(lambda x: os.path.basename(x),gene_tissue_matrix.columns)))
    plt.tight_layout()
    helper_save(f_save)
    if to_close:
        plt.close()
    return


def collapse_on_experimental_type(f_in,landmark_name, f_save=None,to_close=False):
    """ Determines how many tissues per landmark. To do this, needs to collapse values onto the actual tissue.
        Assumes that name in samples columns is 'tissue_...'. Also does the same for samples without CHO
    """
    df = pickle.load(open(f_in,'rb'))
    not_in_cho = []
    in_cho = []
    all_t = []
    for g in df['samples']:
        curr = [i.split('_')[0] for i in g]
        curr = np.unique(np.array(curr))
        if 'CHO' in curr:
            in_cho.append(len(curr))
        else:
            not_in_cho.append(len(curr))
        all_t.append(len(curr))

    f = plt.figure(dpi=300)
    pd.Series(Counter(all_t)).plot.bar(color='y')
    plt.xlabel('Tissues per ' + landmark_name + ' start site')
    plt.ylabel('Number of ' + landmark_name)
    #if f_save is not None:
        #plt.savefig('Results/Figures/num_tissues.pdf', bbox_inches='tight')
    helper_save(f_save)
    if to_close:
        plt.close()
    f = plt.figure(dpi=300)
    pd.Series(Counter(not_in_cho)).plot.bar(color='y')
    plt.xlabel('Tissues per ' + landmark_name + ' start site')
    plt.title('Not in CHO')
    plt.ylabel('Number of ' + landmark_name)
    helper_save(f_save+'_not_in_cho')
    if to_close:
        plt.close()
    #plt.savefig('Results/Figures/num_tissues_not_in_cho.pdf', bbox_inches='tight')
    return


###
# Get the TSS across tissues
###
def plot_tss_across_tissues(f_in, tissues, landmark_name,f_save=None,
                            tissue_list='maxSamples'):
    df = pickle.load(open(f_in,'rb'))
    tissues_genes = dict()
    tissues_genes_unique = dict() # The number of unique genes- i.e
                                  # found in only that tissue
    for t in tissues:
        tissues_genes[t] = 0

    for ind, val in df.iterrows():
        curr_ts = set()
        for t in val[tissue_list]:
            curr_ts.add(t.split('_')[0])
        for t in curr_ts:
            tissues_genes[t] += 1
            # if '1h' in t or 'KLA' in t:
            #     print(t)

    no_peak = []
    for t in tissues_genes:
        if tissues_genes[t] == 0:
            no_peak.append(t)

    for t in no_peak:
        tissues_genes.pop(t, None)

    f, ax = plt.subplots(dpi=300)
    f.patch.set_facecolor("w")
    ax.set_facecolor("w")
    x = range(len(tissues_genes)+1)
    names = list(tissues_genes.keys()) #Add the total number of genes
    names.append('Cumulative fraction')
    y = 1.0*np.array(tissues_genes.values())/(df.shape[0])
    y = np.append(y, [1.0 * np.sum(df['hasGene']) / (df.shape[0])])
    barlist = plt.bar(x, y, align='center')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    barlist[-1].set_color('purple')
    plt.xticks(range(len(tissues_genes)+1), names, rotation=90)
    plt.ylabel('Fraction of ' + landmark_name + ' covered by tissue',{'fontsize': 22})
    plt.title('TSS across tissues',{'fontsize': 22})
    ax.grid("off")
    ax.yaxis.grid(color="grey")
    helper_save(f_save)


def plot_tss_across_tissues_plus_unique(f_in, tissues, landmark_name,
                                f_save=None, tissue_list='Tissues',
                                        is_unique=True):
    df = pickle.load(open(f_in,'rb'))
    tissues_genes = dict()
    tissues_genes_unique = dict() # The number of unique genes- i.e
                                  # found in only that tissue
    for t in tissues:
        tissues_genes[t] = 0
        tissues_genes_unique[t] = 0
    for ind, val in df.iterrows():
        curr_ts = val[tissue_list].split(",")
        if curr_ts[0] == "":
            continue
        for t in curr_ts:
            tissues_genes[t] += 1
        if len(curr_ts) == 1:
            tissues_genes_unique[curr_ts[0]] += 1

    no_peak = []
    for t in tissues_genes:
        if tissues_genes[t] == 0:
            no_peak.append(t)

    for t in no_peak:
        tissues_genes.pop(t, None)

    f,ax = plt.subplots(dpi=300)
    f.patch.set_facecolor("w")
    ax.set_facecolor("w")

    ## Plot
    x = range(len(tissues_genes)+1)
    names = list(tissues_genes.keys()) #Add the total number of genes
    names.append('Cumulative fraction')
    y = 1.0*np.array(tissues_genes.values())/(df.shape[0])
    y = np.append(y, [1.0 * np.sum(df['hasGene']) / (df.shape[0])])
    barlist = plt.bar(x, y, align='center', color="b")


    ## Unique genes
    if is_unique:
        y2 = np.array([1.0*tissues_genes_unique[k] for k in names[
                                                            :-1]])\
                / \
            df.shape[0]
        #y = 1.0 * np.array(tissues_genes_unique.values()) / (
        # df.shape[0])
        y2 = np.append(y2, [1.0 * np.sum(y2)])
        barlist2 = plt.bar(x, y2, align='center',color='g')

        np.savetxt(f_save + ".vals", np.array((y, y2)),
                   delimiter=",")
    else:
        np.savetxt(f_save + ".vals", y, delimiter=",")
    ## Plot settings
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    barlist[-1].set_color('purple')
    plt.xticks(range(len(tissues_genes)+1), names, rotation=90)
    plt.ylabel('Fraction of ' + landmark_name + ' covered by tissue',{'fontsize': 22})
    plt.title('TSS across tissues',{'fontsize': 22})
    ax.grid("off")
    ax.yaxis.grid(color="grey")
    helper_save(f_save)
    return


## Creates a Venn Diagram of genes found in CHO vs other Tissues
def find_unique_cho(f_in, tissue_list='Tissues', f_save=""):
    df = pickle.load(open(f_in, 'rb'))
    tissue_genes = set()
    cho_genes = set()
    for ind, val in df.iterrows():
        curr_ts = val[tissue_list].split(",")
        if curr_ts[0] == "":
            continue
        if "CHO" in curr_ts:
                cho_genes.add(ind)
                curr_ts.remove("CHO")
        for t in curr_ts:
            tissue_genes.add(ind)
            break
    venn2((cho_genes,tissue_genes),set_labels=("CHO", "Other"))
    helper_save(f_save)
    return tissue_genes,cho_genes


def merge_tissues(tissues, tissues_merge:dict):
    """Creates a collapsed list of tissues, where tissues_merge contains keys
    that are found in tissues and the value is the new value to return. This will be a
    many-to-one dictionary"""
    tissues_collapse = set()
    if tissues_merge is None:
        return tissues
    for t in tissues:
        if t not in tissues_merge:
            tissues_collapse.add(t)
        else:
            tissues_collapse.add(tissues_merge[t])
    #print("tissues_collapse", tissues_collapse)
    return list(tissues_collapse)

######
def plot_tss_across_tissues_plus_seen_in_others(df, tissues, landmark_name,
                                                f_save=None,
                                                tissue_list='Tissues',
                                                not_in_cho=False,
                                                verbose=True,
                                                to_plot=False,
                                                tis_n=(), use_markers=False, merge_type=None):
    """Plots barplot for each tissue and how many genes were seen. For each tissue,
    Will also have subbars in there that represents the fraction of the TSS
    that were seen in n tissues or less, where each comes from the tis_n list.
    If it is empty, will just plot the barplots, unless unique is also True, then
    will add which ones are unique (i.e. <= 1 tissue).
    A tissue_genes dict with each tissue as a key and the value is the overall number of genes in that tissue.
    There will be tissue_genes_unique with tissues as keys but each value will be a dictionary,
    where the key will be a value in tis_n, which basically is number of genes seen in
    less than or equal to n tissues that is seen in this tissue. There will also be an overall one as well."""





    tissues_genes = dict()
    tissues_genes_unique = dict()  # The number of unique genes- i.e # found in only that tissue
    tissues_genes_unique_vals = defaultdict(list)
    # Merge samples dictionary
    if merge_type == "bmdm":
        print("Merging bmdm")
        merge_dict = {"BMDMwt": "BMDM", "BMDM1hKLA": "BMDM"}
        tissues = merge_tissues(tissues, merge_dict)
    else:
        merge_dict = {}
    tis_n = list(tis_n)
    tis_n.sort()
    #if to_merge is not None:
    #     for t in tissues:
    #         if t in to_merge:
    #             tissues_genes[to_merge[t]] = 0
    #             tissues_genes_unique[to_merge[t]] = dict()
    #             for n in tis_n:
    #                 tissues_genes_unique[to_merge[t]][n] = 0
    #         else:
    #             tissues_genes[t] = 0
    #             tissues_genes_unique[t] = dict()
    #             for n in tis_n:
    #                 tissues_genes_unique[t][n] = 0
    # else:
    #     for t in tissues:
    #         tissues_genes[t] = 0
    #         tissues_genes_unique[t] = dict()
    #         for n in tis_n:
    #             tissues_genes_unique[t][n] = 0


    for t in tissues:
        tissues_genes[t] = 0
        tissues_genes_unique[t] = dict()
        for n in tis_n:
            tissues_genes_unique[t][n] = 0


    cumulative = dict()
    for n in tis_n:
        cumulative[n] = 0
    cumulative['Overall'] = 0

    # Loop through each gene
    for ind, val in df.iterrows():
        curr_ts = val[tissue_list].split(",")
        if curr_ts[0] == "":
            continue
        curr_ts = merge_tissues(curr_ts, merge_dict)
        curr_n = len(curr_ts)
        for t in curr_ts:
            tissues_genes[t] += 1
            for n in tis_n:
                if curr_n <= n:
                    tissues_genes_unique[t][n] += 1
        for n in tis_n:
            if curr_n <= n:
                cumulative[n] += 1
        if curr_n > 0:
            cumulative['Overall'] += 1
        if curr_n == 1:
            tissues_genes_unique_vals[curr_ts[0]].append(ind)


    ## Double check if there was a tissue with no peak and remove
    no_peak = []
    for t in tissues_genes:
        if tissues_genes[t] == 0:
            no_peak.append(t)
    for t in no_peak:
        tissues_genes.pop(t, None)

    # Calculate the peaks not in CHO
    if not_in_cho:
        tissue_genes = set()
        cho_genes = set()
        for ind, val in df.iterrows():
            curr_ts = val[tissue_list].split(",")
            if curr_ts[0] == "":
                continue
            curr_ts = merge_tissues(curr_ts, merge_dict)
            if "CHO" in curr_ts:
                cho_genes.add(ind)
                curr_ts.remove("CHO")
            for t in curr_ts:
                tissue_genes.add(ind)
                break
        not_in_cho_val = len(tissue_genes - cho_genes)

    #Calculate cumulative fraction:
    #y = 1.0*np.array(tissues_genes.values())/(df.shape[0])
    #tissues_genes["Cumulative"] = np.sum(df['hasGene'])

    # turn into dataframe:
    results_df = pd.DataFrame(tissues_genes_unique)
    results_df = results_df.append(pd.Series(tissues_genes, name="Overall"))
    results_df = pd.concat((results_df, pd.Series(cumulative, name="Cumulative")), axis=1)
    results_df = results_df.iloc[::-1]

    #results_df = results_df[~(results_df.sum()== 0)]

    if to_plot:
        plot_bar(results_df,df.shape[0], landmark_name, f_save=f_save,include_frac=True,
                 use_markers=use_markers)
    if f_save is not None:
        results_df.to_csv(f_save+".csv")
        #outdir = dirname(f_save)
        for t in tissues_genes_unique_vals:
            with open(f_save.replace('.png','') +  t + ".uniqueGenes.txt", 'w') as f:
                f.write("\n".join(tissues_genes_unique_vals[t]))

    return results_df, df.shape[0]


def plot_bar(results_df,num_genes,landmark_name,f_save=None,include_frac=True,use_markers=False):
    color_map, name_map = get_colors('categorical', n_colors=len(results_df), names=results_df.index)
    patterns = ["+" , "-", ".", "*","x", "o", "O","/" ]
    f, ax = plt.subplots(figsize=(12,12))
    count = 0
    for ind, val in results_df.iterrows():
        print(ind)
        x = range(len(val))
        names = list(val.keys()) #Add the total number of genes
        #names.append('Cumulative fraction')
        y = 1.0*np.array(val.values)
        #y = np.append(y, [1.0 * np.sum(val)])

        if ind != 'Overall':
            if use_markers:
                barlist = plt.bar(x, y, align='center', color=color_map[ind],hatch=patterns[count])#,ax=ax)
                count += 1
            else:
                barlist = plt.bar(x, y, align='center',color=color_map[ind])
        else:
            barlist = plt.bar(x, y, align='center', color=color_map[ind])

    #plt.xticks(range(len(tissues_genes)+1), names, rotation=90)
    plt.xticks(range(results_df.shape[1]), names, rotation=90)
    plt.ylabel('Fraction of ' + landmark_name + ' covered by tissue',{'fontsize': 22})
    plt.title('TSS across tissues',{'fontsize': 22})

    #twin-axis for fractin values
    ax2 = ax.twinx()
    ax2.yaxis.grid([])
    ax2.set_ylim(ax.get_ylim())
    ax2.tick_params(labelleft=False,labelright=True)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    vals = ax.get_yticks()
    new = []
    for v in vals:
        new.append(f'{v/num_genes:.2f}')
    ax2.set_yticklabels(new)

    # Change color map names
    color_map_proper_names = dict()
    for c in color_map:
        if type(c) == int:
            color_map_proper_names[f'<= {c} tissues'] = color_map[c]
        else:
            color_map_proper_names[c] = color_map[c]
    legend_from_color(color_map_proper_names, ax, f=None)
    helper_save(f_save)
    return

#######################################################################
## Take a series and the desired width and use the index of the
# series as the location on the x-axis and sum them up over w of them
# e.g. series index is -75,-74...75 and w =10 then -75:-66 will be
# summed and -65:-56 will be summed...
def create_bins(x,w=10):
    instances = x
    inds = x.index.values
    grouping = []
    group_inds = []
    for i in range(0, len(instances), w):
        grouping.append(instances[i:min(len(instances), i + w)].sum())
        group_inds.append(inds[i])

    group_inds.append(inds[-1])

    group_inds_mid = []
    for ind, val in enumerate(group_inds[:-1]):
        group_inds_mid.append(
            (group_inds[ind] + group_inds[ind + 1]) / 2)

    plt.plot(group_inds_mid, grouping)

    return grouping, group_inds_mid


############################################################
### A. Plot the distance to annotated tss
### B. Plot nucleotide frequency
### C. Plot max Stat value
############################################################


###
# C. Plots the value from the bed file
###
def plot_stat(f_bed_in,f_save=None):
    df = pd.read_csv(f_bed_in,sep='\t')
    sns.distplot(np.log2(df[4]))
    plt.title('Experimental genetic transcription start sites')
    plt.xlabel('Log2 TSS count per million')
    plt.ylabel('Frequency')
    helper_save(f_save)