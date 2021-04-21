import pandas as pd
from pandarallel import pandarallel
import click
import numpy as np
from tss.manuscript.utils import expand_anno_id


# def get_distance_anno(peaks_df, in_anno_df, nb_workers=32):
#     """ Get nearest distance and if same strand
#
#     :param peaks_df:
#     :return:
#     """
#     pandarallel.initialize(nb_workers=nb_workers)
#     anno_chr_str = in_anno_df.set_index(['Chr', "Strand"])
#
#     def get_nearest_info(x, isSame, col_id="transcript_id"):
#         # peaks of one gene
#         print('x')
#         print(x)
#         curr_str = x['Strand']
#         assert (curr_str == "+" or curr_str == "-")
#         if not isSame:
#             if curr_str == "-":
#                 curr_str = "+"
#             else:
#                 curr_str = "-"
#             name='TSS_Rev'
#         else:
#             name="TSS"
#
#         print('key', (x['Chr'], curr_str))
#         print((x['Chr'], curr_str))
#         if (x['Chr'], curr_str) not in anno_chr_str.index:
#             return pd.Series([np.nan, ""],
#                              index=['Distance to Nearest '+name,
#                                     f'Nearest ' + name])
#         curr_anno = anno_chr_str.loc[(x['Chr'], curr_str)].reset_index()
#         dist_to_peak = x-curr_anno['actual_start']
#         min_idx = np.abs(dist_to_peak).idxmin()
#         min_dist = dist_to_peak.loc[min_idx]
#         min_id = curr_anno.loc[min_idx, col_id]
#         return pd.Series([min_dist, min_id], index=[f'Distance to Nearest {name}', f'Nearest {name}'])
#
#     #nearest_peaksFw  = peaks_df.parallel_apply(get_nearest_info, axis=1, args=(True,))
#     print('peaks_df')
#     print(peaks_df.head())
#     nearest_peaksFw = peaks_df.apply(get_nearest_info, axis=1,
#                                               args=(True,))
#     nearest_peaksRev = peaks_df.apply(get_nearest_info, axis=1,
#                                                args=(False,))
#     #nearest_peaksRev = peaks_df.parallel_apply(get_nearest_info, axis=1, args=(False,))
#     peaks_anno_dist = pd.merge(nearest_peaksFw, nearest_peaksRev, left_index=True, right_index=True)
#     return peaks_anno_dist


def get_distance_anno(peaks_df, in_anno_df, nb_workers=32):
    """ Get nearest distance and if same strand

    :param peaks_df:
    :return:
    """

    anno_chr_str = in_anno_df.sort_values(['Chr', 'Strand']).set_index(
        ['Chr', "Strand"])
    assert ((peaks_df['Start'] - peaks_df['End'] == 0).all())

    def get_nearest_info(x, isSame, col_id="transcript_id"):
        # peaks of one gene
        curr_str = x['Strand']
        assert (curr_str == "+" or curr_str == "-")
        if not isSame:
            if curr_str == "-":
                curr_str = "+"
            else:
                curr_str = "-"
            name = 'TSS_Rev'
        else:
            name = "TSS"

        if (x['Chr'], curr_str) not in anno_chr_str.index:
            return pd.Series([np.nan, ""],
                             index=['Distance to Nearest ' + name,
                                    f'Nearest ' + name])
        curr_anno = anno_chr_str.loc[(x['Chr'], curr_str)].reset_index()
        dist_to_peak = x['Start'] - curr_anno['actual_start']

        min_idx = np.abs(dist_to_peak).idxmin()
        min_dist = dist_to_peak.loc[min_idx]
        min_id = curr_anno.loc[min_idx, col_id]
        return pd.Series([min_dist, min_id],
                         index=[f'Distance to Nearest {name}',
                                f'Nearest {name}'])
    nearest_peaksFw = peaks_df.parallel_apply(get_nearest_info, axis=1,
                                              args=(True,))
    nearest_peaksRev = peaks_df.parallel_apply(get_nearest_info, axis=1,
                                               args=(False,))
    peaks_anno_dist = pd.merge(nearest_peaksFw, nearest_peaksRev,
                               left_index=True, right_index=True)
    return peaks_anno_dist


def create_anno_tss(refTSS_f, gff_f, pc, ncrna=(), break_char="="):

    # Read in annotation file
    if '.gtf' in gff_f or '.gtf.gz' in gff_f:
        print(
            "Only works for gff files, not gtf. Convert using gffread from the command line.")
    genome_ann = pd.read_csv(gff_f, comment='#', sep='\t', header=None)
    col_names = list(genome_ann.columns)
    print('col_names', col_names)
    col_names[:9] = ['Chr', 'How', 'Annotation', 'Start', 'End', '.1',
                     'Strand', '.2', 'ID']
    genome_ann.columns = col_names
    genome_ann.sort_values(['Chr', 'Start', 'End'], inplace=True)

    genome_ann = genome_ann[
        genome_ann['Annotation'].isin(list(pc) + list(ncrna))]
    # genome_ann = genome_ann[['ID', 'Chr', 'Start', 'End', 'Strand', 'Annotation']]
    # Make the center of the peak the start site with flanking base pairs
    genome_ann['actual_start'] = genome_ann.apply(
        lambda x: x['Start'] if x['Strand'] == '+' else x['End'],
        axis=1)
    ## If the strand is positive, then the end moves to the start, but if its negative, the start moves to the end.
    genome_ann['End'] = genome_ann['actual_start']
    genome_ann['Start'] = genome_ann['actual_start']

    genome_ann = expand_anno_id(genome_ann, colname="ID")
    genome_ann = genome_ann.set_index('transcript_id', drop=False)
    genome_ann.index.name = None
    genome_ann.to_csv(refTSS_f, sep="\t")
    return genome_ann


def homeranno(peak_f, ref_fa, ref_anno, anno_f):
    from tss.utils import Homer
    Homer.annotate_peaks(peak_file=peak_f, output_file=anno_f,
                         ref_fa=ref_fa, ref_anno=ref_anno)



@click.group()
def cli1():
    print('i')
    pass


@cli1.command('genome')
@click.argument('out_f', type=click.Path(exists=False))
@click.argument('gff', type=click.Path(exists=True))
@click.argument('pc', type=click.STRING)
@click.argument('ncrna', type=click.STRING, default="")
@click.option('--breakchar', type=click.STRING, default="=")
def create_genome(out_f, gff, pc, ncrna, breakchar):
    print('generating genome')
    print('gff', gff)
    create_anno_tss(out_f, gff, pc.split(','),
                    ncrna=ncrna.split(','), break_char=breakchar)
    return

@click.group()
def cli2():
    pass

@cli2.command('distance')
@click.argument('peaks_f', type=click.Path(exists=True))
@click.argument('reftss_f', type=click.Path(exists=True))
@click.argument('out_peaks_distance', type=click.Path(exists=False))
@click.option('--cpu', default=32)
def create_distance(peaks_f, reftss_f, out_peaks_distance, cpu):
    pandarallel.initialize(nb_workers=cpu, verbose=2, progress_bar=True)
    reftss = pd.read_csv(reftss_f, sep='\t', index_col=0)
    reftss['actual_start'] = reftss['actual_start'].astype(int)
    peaks_df = pd.read_csv(peaks_f, sep="\t", index_col=0)
    # peaks_df['actual_start'] = peaks_df['actual_start'].astype(int)
    peaks_anno_dist = get_distance_anno(peaks_df, reftss,
                                        nb_workers=cpu)
    peaks_anno_dist.to_csv(out_peaks_distance, sep='\t')
    return





cli = click.CommandCollection(sources=[cli1, cli2])
if __name__ == '__main__':
    cli()
