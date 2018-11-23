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


#################################################
def get_tss_files(data_folder, t, s, tissues_with_RNA, RNA_dir=''):
    """
    Function to get list of tss,input, and rna tag files. Deals with the annoying naming
    """
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


#################################################
### Plot helpers
#################################################
def determine_rows_cols(num_samples):
    nrows = 2
    ncols = 2
    while nrows * ncols <= num_samples:
        ncols = ncols + 1
        if nrows * ncols <= num_samples:
            nrows = nrows + 1

    return nrows,ncols