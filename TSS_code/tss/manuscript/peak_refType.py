import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import click


def create_refType(peaks_distance,refTSS,homeranno, refType,
                   intergenic_dist,bins, pc, ncrna):
    df = pd.read_csv(peaks_distance, sep='\t', index_col=0)
    ref_df = pd.read_csv(refTSS, sep='\t', index_col=0)
    homeranno_df = pd.read_csv(homeranno, sep='\t', index_col=0)

    print('df before dropping no nearest tss', df.shape)
    df = df.dropna(subset=["Nearest TSS"])
    print('df after dropping no nearest tss', df.shape)
    df['Annotation'] = df["Nearest TSS"].map(
        ref_df["Annotation"].to_dict())

    df['homeranno'] = homeranno_df.loc[df.index, 'Annotation']

    # get intergenic
    intergenic = df[(df['homeranno'].str.upper() == "INTERGENIC") & (
                df['Distance to Nearest TSS'] < intergenic_dist)]


    promoter = df[(df['Distance to Nearest TSS'] > bins[0]) & (
                df['Distance to Nearest TSS'] < bins[1])]
    pc_promoters = promoter[promoter['Annotation'] == pc]
    ncrna_promoters = promoter[promoter['Annotation'].isin(ncrna)]
    print(pc_promoters.shape)
    pc_promoters.head()

    print('Total', pc_promoters.shape[0])
    print('annotated as promoter by homer',
          pc_promoters['homeranno'].str.contains("promoter").sum())

    print(ncrna_promoters.shape)
    ncrna_promoters.head()

    nonpromoters = df.drop(promoter.index)
    divergent = nonpromoters[
        (nonpromoters['Distance to Nearest TSS_Rev'] > bins[0]) & (
                    nonpromoters['Distance to Nearest TSS_Rev'] < bins[
                1])]
    print('divergent', divergent.shape)

    df['refType'] = None
    df.loc[intergenic.index, 'refType'] = "Intergenic"
    df.loc[pc_promoters.index, 'refType'] = "mRNA"
    df.loc[ncrna_promoters.index, 'refType'] = "ncRNA"
    df.loc[divergent.index, 'refType'] = "divergent"

    df.to_csv(refType, sep='\t')
    sns.countplot(df['refType'])
    plt.savefig(refType + ".png")
    return


def add_gene_body(df, cds_f, intron_f):
    pd.read_csv(cds_f)

    return

@click.command()
@click.argument('peaks_distance', type=click.Path(exists=True))
@click.argument('reftss', type=click.Path(exists=True))
@click.argument('homeranno', type=click.Path(exists=True))
@click.argument('reftype', type=click.Path())
@click.argument('intergenic_dist', type=click.INT)
@click.argument('bins', type=click.STRING)
@click.argument('pc', type=click.STRING)
@click.argument('ncrna', type=click.STRING)
def main(peaks_distance, reftss, homeranno, reftype,
                   intergenic_dist, bins, pc, ncrna):
    ncrna = ncrna.split(',')
    bins = [int(x) for x in bins.split(",")]
    if bins[0]>0:
        bins[0] = -1*bins[0]
    print('bins', bins)
    #intergenic_dist = int(intergenic_dist)
    create_refType(peaks_distance, reftss, homeranno, reftype,
                   intergenic_dist, bins, pc, ncrna)
    return

if __name__== "__main__":
    main()
