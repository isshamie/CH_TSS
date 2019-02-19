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
from os.path import basename
mpl.style.use('fivethirtyeight')
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')


def plot_cdf_rna(exon_peaks, tss_peaks, peaks, name, is_save=True):
    f = plt.figure()
    num_bins = 20

    counts, bin_edges = np.histogram(exon_peaks['Log2Ratio vs. RNA'], bins=num_bins, normed=True)
    cdf = np.cumsum(counts)
    plt.plot(bin_edges[1:], cdf / cdf[-1], label='Exon')

    counts, bin_edges = np.histogram(tss_peaks['Log2Ratio vs. RNA'], bins=num_bins, normed=True)
    cdf = np.cumsum(counts)
    plt.plot(bin_edges[1:], cdf / cdf[-1], color='r', label='TSS')

    counts, bin_edges = np.histogram(peaks['Log2Ratio vs. RNA'], bins=num_bins, normed=True)
    cdf = np.cumsum(counts)
    plt.plot(bin_edges[1:], cdf / cdf[-1], color='g', label='TSS')

    plt.legend()
    plt.title(name + ' Log2Ratio Start-seq over Total RNA CDF')
    if is_save:
        plt.savefig('Results/csRNATSS/Figures/rna_' + name + '.png')
        plt.close()
    return


def plot_cdf_input(exon_peaks, tss_peaks, peaks, name, is_save=True):
    f = plt.figure()
    num_bins = 20
    counts, bin_edges = np.histogram(exon_peaks['Log2Ratio vs. Input'], bins=num_bins, normed=True)
    cdf = np.cumsum(counts)
    plt.plot(bin_edges[1:], cdf / cdf[-1], label='Exon')

    counts, bin_edges = np.histogram(tss_peaks['Log2Ratio vs. Input'], bins=num_bins, normed=True)
    cdf = np.cumsum(counts)
    plt.plot(bin_edges[1:], cdf / cdf[-1], color='r', label='TSS')

    counts, bin_edges = np.histogram(peaks['Log2Ratio vs. Input'], bins=num_bins, normed=True)
    cdf = np.cumsum(counts)
    plt.plot(bin_edges[1:], cdf / cdf[-1], color='g', label='TSS')

    plt.legend()
    plt.title(name + ' Log2Ratio Start-seq over Input CDF')
    if is_save:
        plt.savefig('Results/csRNATSS/Figures/input_' + name + '.png')
        plt.close()
    return


def plot_sdf_rna(exon_peaks, tss_peaks, peaks, name, is_save=True, num_bins=100):
    f = plt.figure()

    ## Setup bin edges
    num_bins = 100
    concat = np.concatenate((np.array(exon_peaks['Log2Ratio vs. RNA']), np.array(peaks['Log2Ratio vs. RNA'])))
    step = np.ceil(np.ceil(concat.max()) - np.floor(concat.min())) / num_bins
    bins = np.arange(np.floor(concat.min()), np.ceil(concat.max()), step)

    ## Exon peaks
    counts_exon, bin_edges = np.histogram(exon_peaks['Log2Ratio vs. RNA'], bins=bins, density=True)
    cdf_exon = np.cumsum(counts_exon)
    cdf_exon = cdf_exon / cdf_exon[-1]
    sdf_exon = 1 - cdf_exon
    plt.plot(bin_edges[1:], sdf_exon, label='Exon')

    ## TSS peaks
    counts_tss, bin_edges = np.histogram(peaks['Log2Ratio vs. RNA'], bins=bins, density=True)
    cdf_tss = np.cumsum(counts_tss)
    cdf_tss = cdf_tss / cdf_tss[-1]
    sdf_tss = 1 - cdf_tss
    plt.plot(bin_edges[1:], sdf_tss, color='g', label='TSS')

    plt.plot(bin_edges[1:], sdf_tss - sdf_exon, color='gray', label='Diff')

    plt.legend()
    plt.title(name + ' Log2Ratio Start-seq over Total RNA SDF')

    if is_save:
        plt.savefig('Results/csRNATSS/Figures/SDF_rna_' + name + '.png')
        plt.close()
    return


def plot_sdf_input(exon_peaks, tss_peaks, peaks, name, is_save=True, num_bins=100):
    f = plt.figure()

    ## Setup bin edges
    concat = np.concatenate((np.array(exon_peaks['Log2Ratio vs. Input']), np.array(peaks['Log2Ratio vs. Input'])))
    step = np.ceil(np.ceil(concat.max()) - np.floor(concat.min())) / num_bins
    bins = np.arange(np.floor(concat.min()), np.ceil(concat.max()), step)

    ## Exon peaks
    counts_exon, bin_edges = np.histogram(exon_peaks['Log2Ratio vs. Input'], bins=bins, density=True)
    cdf_exon = np.cumsum(counts_exon)
    cdf_exon = cdf_exon / cdf_exon[-1]
    sdf_exon = 1 - cdf_exon
    plt.plot(bin_edges[1:], sdf_exon, label='Exon')

    ## TSS peaks
    counts_tss, bin_edges = np.histogram(peaks['Log2Ratio vs. Input'], bins=bins, density=True)
    cdf_tss = np.cumsum(counts_tss)
    cdf_tss = cdf_tss / cdf_tss[-1]
    sdf_tss = 1 - cdf_tss
    plt.plot(bin_edges[1:], sdf_tss, color='g', label='TSS')

    plt.plot(bin_edges[1:], sdf_tss - sdf_exon, color='gray', label='Diff')

    plt.legend()
    plt.title(name + ' Log2Ratio Start-seq over Total Input SDF')

    if is_save:
        plt.savefig('Results/csRNATSS/Figures/SDF_input_' + name + '.png')
        plt.close()
    return