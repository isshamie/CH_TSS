import glob
import os
import numpy as np
import pandas as pd
import snakemake


#################################################
def read_peak_file(peak_f):
    df = pd.read_csv(peak_f,sep='\t',skiprows=36, index_col=0)
    cols = df.columns.values
    cols[:4] = ['Chr','Start','End','Strand']
    df.columns = cols
    return df


# def read_bed_file(bed_f):
#     df = pd.read_csv(bed_f,sep='\t', header=None)
#     df.columns = ["Chr", "Start","End","ID","Stat","Strand"]
#     df = df.set_index("ID")
#     df = df[['Chr', 'Start', 'End', 'Strand', 'Stat']]
#     return df


def read_bed_file(bed_f):
    try:
        df = pd.read_csv(bed_f,sep='\t', header=None)
    except pd.errors.ParserError:
        df = pd.read_csv(bed_f,sep='\t', header=None, skiprows=1)
    columns = df.columns.values.astype(str)
    columns[:6] = ["Chr", "Start","End","ID","Stat","Strand"]
    df.columns = columns
    df =  df.set_index("ID")
    #df = df[['Chr', 'Start', 'End', 'Strand', 'Stat']]
    return df


def write_bed_file(df, bed_f, use_index=True):
    if use_index:
        df["ID"] = df.index
    df = df[['Chr', 'Start', 'End', "ID", 'Stat', 'Strand']]

    df.to_csv(bed_f, sep="\t", header=None, index=None)


def bed_subset_region(bed_f, f_save=None, region=(-2, +3)):
    bed = read_bed_file(bed_f)
    center = np.ceil((bed["Start"] + bed["End"]) / 2)
    bed["Start"] = (center + region[0]).astype(int)
    bed["End"] = (center + region[1]).astype(int)
    print(bed.columns)
    print(bed.index.name)
    if f_save is not None:
        write_bed_file(bed, f_save)

    return bed


def get_tss_files(data_folder, t, s, tissues_with_RNA, RNA_dir=''):
    """
    Function to get a list of tss,input, and rna tag files. Deals
    with the annoying naming
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
            rna_tag = (np.array(rna_tag)[[os.path.isdir(x) for x in rna_tag]])[0]  # Not the bam but directory

    return tss_tag, inp_tag, rna_tag
