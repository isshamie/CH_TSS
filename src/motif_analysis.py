import numpy as np
import pickle
import tqdm
import sys
import os
import pandas as pd
## Bio import
from Bio import motifs as mt
from Bio.Seq import Seq
from Bio import SeqIO

sys.path.append("/home/isshamie/software/homebrew/parallel_functions/")
import parallel_functions as pf
import time

def read_pssm(motif_file ,background={'A' :0.25 ,'C' :0.25 ,'G' :0.25 ,'T' :0.25}):
    ''' Loads pfm motif file and coverts into pssm '''
    with open(motif_file) as handle:
        m = mt.read(handle, "pfm")
    pwm = m.counts.normalize(pseudocounts=0.5)
    ## Get pssm with backgroud
    pssm = pwm.log_odds(background)
    return pssm


def compute_score(pssm ,seq ,norm_len=-1):
    ''' Computes the pssm score for each position.
    Params:
    pssm: Bio.motifs.matrix.PositionSpecificScoringMatrix
    seq: Bio.Seq.Seq or string.
    norm_len: If not -1, it is length of the sequene to be used.
              If the length of seq is greater than norm_len, then the center plus (norm_len/2) on the left and norm_len/2 on the right of center.
    '''
    if type(seq) == str:
        seq = Seq(seq ,pssm.alphabet)

    if not norm_len == -1 and len(seq) > norm_len:
        center = len(seq)
    if len(pssm[0]) > len(seq):
        return np.zeros(len(seq))
    return pssm.calculate(seq)


def search_for_instance(pssm,seq):
    ''' Searches for all the instances of the pssm in the sequence.
        Returns a list of positions and sequences, where positions is where it happened and
        sequences are the seqs that passed the threshold'''
    return list(zip(*list(pssm.search(seq, threshold=3.0))))

    # for pos,score in pssm.search(seq, threshold=3.0):
    #     print('pos',pos)
    #     #print(score)
    #     positions.append(pos)
    #     scores.append(scores)
    #     # print pos, seq
    # return positions, score


def wrap_fa_motif(fa_file ,motif_file,f_save=None):
    '''Reads in fasta file and motif file and computes the score of the pssm starting at each position in each sequence.
    Saves to f_save as pickle file
    Returns:
    motifs: A dict where each id from the fa_file is a key and the value is a vector of scores where each element
    is the score of the motif starting at that position in the sequence.
    '''
    pssm = read_pssm(motif_file)
    motifs = dict()
    t0 = time.time()
    for record in tqdm.tqdm(SeqIO.parse(fa_file, "fasta" ,alphabet = pssm.alphabet)):
        # print(record.seq)
        # seq = Seq(record.seq,pssm.alphabet)
        score = compute_score(pssm,record.seq)
        motifs[record.id] = score
    if f_save is not None:
        pickle.dump(motifs ,open(f_save ,'wb'))
    print("Time taken: %d" % (time.time() - t0))
    return motifs


def wrap_fa_motif_instance(fa_file ,motif_file ,f_save=None):
    '''Reads in fasta file and motif file and finds all instances of the motif for each sequence.
    Saves to f_save as pickle file
    Returns:
    motifs: A dict where each id from the fa_file is a key and the value is the positions list of
    where an instance of the motif occurred.
    '''
    t0 = time.time()
    pssm = read_pssm(motif_file)
    motifs = dict()
    for record in tqdm.tqdm(SeqIO.parse(fa_file, "fasta" ,alphabet = pssm.alphabet)):
        positions = search_for_instance(pssm ,record.seq)
        motifs[record.id] = positions
    if f_save is not None:
        pickle.dump(motifs ,open(f_save ,'wb'))
    print("Time taken: %d" % (time.time() - t0))
    return motifs


def wrap_both_instance_and_score(fa_file,motif_file,f_save=None):
    t0 = time.time()
    pssm = read_pssm(motif_file)
    motifs = dict()
    motifs_instance = dict()
    for record in tqdm.tqdm(SeqIO.parse(fa_file, "fasta" ,alphabet = pssm.alphabet)):
        positions ,scores = search_for_instance(pssm ,record.seq)
        motifs_instance[record.id] = positions
        score = compute_score(pssm ,record.seq)
        motifs[record.id] = score

    if f_save is not None:
        pickle.dump(motifs,open(f_save + '_score','wb'))
        pickle.dump(motifs, open(f_save + '_instance', 'wb'))
    print("Time taken: %d" % (time.time() - t0))
    return motifs


def all_motifs_of_interest_instance(motifs_list,fa_file, motifs_folder, f_save):
    for i in motifs_list:
        print(i)
        f_out = '%s_%s_instance.p' % (f_save,i)
        motif_f = motifs_folder + i + '.motif.pfm'
        wrap_fa_motif_instance(fa_file,motif_f,f_save=f_out)
    return


def all_motifs_of_interest_score(motifs_list, fa_file, motifs_folder, f_save):
    for i in motifs_list:
        print(i)
        f_out = '%s_%s.p' % (f_save,i)
        motif_f = motifs_folder + i + '.motif.pfm'
        wrap_fa_motif(fa_file,motif_f,f_save=f_out)
    return



################
# Motif Results analysis
##############
def motif_counts(seq_len, motif_instances, bins=(-75, 75), compl_strand=False):
    """ Returns counts of motifs within the a certain range.
        If compl_strand=True, it returns a tuple with number of positive hits as the first element
        And negative elements as the second.
    """
    mid = np.ceil(seq_len / 2)
    curr_motif_instance = np.array(motif_instances)  # From list to array

    if compl_strand:
        pos = curr_motif_instance[np.apply_along_axis(lambda x: x >= 0, 0, curr_motif_instance)]
        neg = curr_motif_instance[np.apply_along_axis(lambda x: x < 0, 0, curr_motif_instance)]
        pos_count = np.sum(np.apply_along_axis(lambda x: (bins[0] < mid - x) & (mid - x < bins[1]), 0, pos))
        neg = np.abs(neg)
        neg_count = np.sum(np.apply_along_axis(lambda x: (bins[0] < mid - x) & (mid - x < bins[1]), 0, neg))
        return (pos_count, neg_count)
        # return list(map(lambda x: bins[0] < mid - x and mid-x < bins[1] ,motif_instances)).count(True)
    else:
        return np.sum(np.apply_along_axis(lambda x: (bins[0] < mid - x) & (mid - x < bins[1]), 0, curr_motif_instance))


def wrap_motifs_counts(df,anno_df,motifs,new_col):
    not_in_peaks = 0
    for peak, val in df.iterrows():
        seq_len = anno_df.loc[peak,'Length']
        if peak in motifs:
            motif_instances = motifs[peak]
            #motif_instances = df.loc[peak,'instances']
            df.at[peak,new_col] = motif_counts(seq_len,motif_instances,compl_strand=True)[0]
        else:
            df.at[peak,new_col] = 0
            not_in_peaks += 1
    #print('Not found in peaks', not_in_peaks)
    return df


def create_peak_by_motif_df(all_motifs, f_anno, f_save=None):
    #all_motifs = glob.glob('Results/motifs/all_peaks_merged_anno*instance.p')
    anno_peaks = pd.read_csv(f_anno, sep='\t', index_col=0)
    anno_peaks['Length'] = anno_peaks['End'] - anno_peaks['Start']
    motif_names = map(lambda x: os.path.basename(x).split('_')[-2],all_motifs)
    peak_motif_counts = pd.DataFrame(index=anno_peaks.index,columns = motif_names)
    for ind,i in enumerate(all_motifs):
        print(ind)
        print(i)
        f_out = i
        motifs = pickle.load(open(f_out,'rb'))

        print(motif_names[ind])
        #peak_motif_counts = wrap_motifs(anno_peaks,motifs,motif_names[ind])
        peak_motif_counts = pf.parallel_df(peak_motif_counts,func=wrap_motifs_counts,func_args=[anno_peaks[['Length']],motifs,motif_names[ind]],num_processes=2)
    if f_save is not None:
        peak_motif_counts.to_csv(f_save)
    return peak_motif_counts