import pandas as pd
import numpy as np
from cycler import cycler
from collections import Counter
from matplotlib import pyplot as plt
import matplotlib as mpl
import pickle
import os
from matplotlib import rcParams
import seaborn as sns
from matplotlib_venn import venn2

rcParams['figure.figsize'] = 8, 6
#mpl.style.use('fivethirtyeight')
mpl.style.use('ggplot')
#mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')


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

    f,ax = plt.subplots(dpi=300)
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
    barlist = plt.bar(x, y, align='center')

    ## Unique genes
    if is_unique:
        y2 = np.array([1.0*tissues_genes_unique[k] for k in names[
                                                            :-1]])\
                / \
            df.shape[0]
        #y = 1.0 * np.array(tissues_genes_unique.values()) / (
        # df.shape[0])
        y2 = np.append(y2, [1.0 * np.sum(y2) / (df.shape[0])])
        barlist2 = plt.bar(x, y2, align='center',color='g')

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


def helper_save(f_save):
    """ Function to save as png and svg if f_save is not None."""
    if f_save is not None:
        if '.png' in f_save:
            plt.savefig(f_save)
            plt.savefig(f_save.split('.png')[0] + '.svg')
        else:
            plt.savefig(f_save + '.png',bbox_inches='tight')
            plt.savefig(f_save + '.svg')


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