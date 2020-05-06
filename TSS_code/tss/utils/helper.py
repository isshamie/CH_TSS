import yaml
import sys
import os
import time
import sarge


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
### Plot helpers
#################################################
def determine_rows_cols(num_samples):
    nrows = 2
    ncols = 2
    while nrows * ncols < num_samples:
        ncols = ncols + 1
        if nrows * ncols <= num_samples:
            nrows = nrows + 1

    return nrows,ncols
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



def mdir(d, replace=True, save_old=True):
    """ Function to make a directory. If the directory exists,
    it will either do nothing or replace the directory and save the
    old in a subfolder depending on the parameters
    Parameters:
        @d: type: str
            directory name.
        @replace: type=bool
                  default=True
                  To replace the old file or not
        @save_old: type=bool
                   default=True
                   Only on if replace=True. If true, will save the old
                   directory inside
                   'd/old_{index}' where the index will be the old
                   numbered version.

    Returns:
        True if the folder was not already present OR was
        replaced.
    """
    if os.path.exists(d):
        if replace:
            if save_old:
                olds = glob.glob(os.path.join(d, "old_*"))
                curr_index = len(olds)
                count = 1
                new_old = "old_%d" % (curr_index + count)
                new_old = os.path.join(d, new_old)

                # Keep incrementing until a new directory is found.
                while os.path.exists(new_old):
                    count += 1
                    new_old = new_old[:new_old.rfind("_") +1] + str(
                        count)

                os.mkdir(new_old)
                print("Directory already there. Moving contents into "
                      "%s" %
                      new_old)

                #This command moves everything but the old folder
                cmd = "mv `ls -1 %s/* | grep -v old_` %s" % (d,
                                                              new_old)
                print(cmd)
                os.system(cmd)
                time_file(new_old, name="DATE_moved.txt")

            else:
                cmd = "rm -rf %s/* " % d
                os.system(cmd)
        else:
            print("Directory already there. Keeping as is")
            return False
    else:
        os.makedirs(d)
    return True


def time_file(outdir, name="DATE"):
    """ Function that generates a simple text file with a time stamp
    up to hours and minutes.
        @outdir: directory to save file in
        @name: file name to save to
        """
    t = time.localtime()
    timestamp = time.strftime('%b-%d-%Y_%H%M', t)
    with open(os.path.join(outdir, name), 'w') as f:
        f.write(timestamp)
    return


def append_name(f, prefix="", suffix="", to_rm=(".csv",".tsv")):
    """ Appends a prefix and suffix to string to the beginning of a
    filename
    that may be a full path """

    for i in to_rm:
        f = f.replace(i, "")

    if not prefix == "":
        prefix = prefix + "_"
    if not suffix == "":
        suffix = "_" + suffix

    f = os.path.join(os.path.dirname(f), prefix + os.path.basename(f)
                     + suffix)
    return f