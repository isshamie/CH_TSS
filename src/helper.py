import yaml
import sys
import os

## First go back up a folder
import os
import sys
import pandas as pd
import matplotlib
import seaborn as sns
import pickle
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from itertools import product
import glob
import re
from matplotlib import rcParams
import inspect
rcParams['figure.figsize'] = 8, 6
import tqdm
### notebook specific configuration ###
from os.path import basename
#mpl.style.use('ggplot')
mpl.style.use('fivethirtyeight')
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')



def output_of_findPeaks(in_f, type_f , out_f = ''):
    """Function that takes the output of Homer (tested with v4.10, and 4.10.3, 5-16-2018)
    from findPeaks or mergePeaks and creates a more standard output and saves it as .tsv file

    Input:
    1) in_f: Name of the input file
    2) type_f: 'merge' or 'orig', where merge means the file comes from merged peaks,
                and orig means it comes from mergePeaks
    3) out_f: file name to save to. If empty, out_f will be in_f + '_STD', where STD is standard.

    """

    if type_f == 'orig':
        df = pd.read_csv(in_f,sep='\t',skiprows = 36)

        df = df.rename(columns={'start': 'Start', 'end': 'End', 'strand': 'Strand', 'chr': 'Chr', '#PeakID': 'ID',
                           'Normalized Tag Count': 'Stat'})
    return 42

def merge_files():
    return 42


#################################################
def get_tss_files(data_folder, t, s, RNA_dir=''):
    '''
    Function to get list of tss,input, and rna tag files. Deals with the annoying naming
    '''
    if s == 'GRO':
        tss_tag = glob.glob(
            os.path.join(data_folder, t, s, 'f03_tags/') + '/*5GRO*')
        inp_tag = glob.glob(
            os.path.join(data_folder, t, s, 'f03_tags') + '/*_GRO*')
    else:
        tss_tag = glob.glob(
            os.path.join(data_folder, t, s, 'f03_tags/') + '/*_mSTART_*')
        inp_tag = glob.glob(
            os.path.join(data_folder, t, s, 'f03_tags') + '/*mSTARTinput*')

    # Naming is different for CHO
    if t == 'CHO':
        if s == 'GRO':
            tss_tag = glob.glob(
                os.path.join(data_folder, t, s, 'f03_tags/') + '/*5GRO*')
            inp_tag = glob.glob(
                os.path.join(data_folder, t, s, 'f03_tags') + '/*-GRO*')
        else:
            tss_tag = glob.glob(
                os.path.join(data_folder, t, s, 'f03_tags/') + '/*-mSTART-*')
            inp_tag = glob.glob(
                os.path.join(data_folder, t, s, 'f03_tags') + '/*mSTART_input*')

    rna_tag = ''
    if not RNA_dir == '' and t in tissues_with_RNA:
        comb = glob.glob(os.path.join(RNA_dir, 'Combined*' + tissues_with_RNA[t]))
        if len(comb) == 1:  # glob.glob('Combined*Brain')
            rna_tag = comb[0]
        else:
            rna_tag = glob.glob(os.path.join(RNA_dir, '*' + tissues_with_RNA[t] + '*'))
            rna_tag = (np.array(rna_tag)[map(lambda x: os.path.isdir(x), rna_tag)])[0]  # Not the bam but directory

    return tss_tag, inp_tag, rna_tag
