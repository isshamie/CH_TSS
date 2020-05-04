from tss.utils import Homer

from tss.visualize import plot_tss_results
import yaml
import sys
import os
import sys
import pandas as pd
import seaborn as sns
import pickle
import numpy as np
import matplotlib as mpl
from os.path import basename


def run(p):
    p_stage = p["histogram"]
    p_stage_f = p_stage["filenames"]
    p_reference = p["reference"]
    p_global = p["global"]

    tissues = p_global["TISSUES"]

    # Filenames:
    gene_f_in = p["gene_centric"]["gene centric peaks from samples plus tissues"]
    txn_f_in = p["gene_centric"]["txn centric peaks from samples plus tissues"]

    gene_f_save = p_stage_f["gene fraction unique"]
    venn_gene_f_save = p_stage_f["gene venn"]
    txn_f_save = p_stage_f["txn fraction unique"]
    venn_txn_f_save = p_stage_f["txn venn"]

    plot_tss_results.plot_tss_across_tissues_plus_unique(gene_f_in, tissues,
                                                         landmark_name='genes',f_save=gene_f_save)

    plot_tss_results.find_unique_cho(gene_f_in,f_save=venn_gene_f_save)

    plot_tss_results.plot_tss_across_tissues_plus_unique(txn_f_in, tissues, landmark_name='Transcript',f_save=txn_f_save)
    plot_tss_results.find_unique_cho(txn_f_in,f_save=venn_txn_f_save)
    return
