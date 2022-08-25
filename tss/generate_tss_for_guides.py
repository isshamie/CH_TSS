import sys
import glob
import yaml
import os
import pandas as pd
import helper
import time
from create_output import *


#######################################################################
#######################################################################
def info_on_gene(goi, peaks_df, tss_df, start_anno, thresh=150):
    # Chr, Start, End, ID, Evidence, Strand, Distance to TSS
    tss_df = tss_df[
        ["Chr", "Start", "End", "ID", "Strand", "Gene", "Transcript",
         "Is Experimental"]].copy()
    tss_df["Evidence"] = ""

    gene_struct = pd.DataFrame(columns=tss_df.columns)

    curr_g = tss_df[tss_df["Gene"] == goi]
    isos = curr_g.groupby("Transcript")
    #print("Number of Isoforms", len(isos))

    for ind, val in isos:
        is_seen = False
        #If theres at least one exp TSS nearby then dont need to use it
        #print(val)

        if (peaks_df.loc[val["ID"].dropna(), "Distance to TSS"] <
            thresh).sum() \
                > 0:
            #print(val)
            is_seen = True

        for ind2, val2 in val.iterrows():
            curr = val2.copy()
            if "p0@" in ind2:
                # Add to dataframe
                curr["Evidence"] = "No experimental peak."
                gene_struct = pd.concat(
                    (gene_struct, pd.DataFrame(curr).transpose()),
                    sort=True, axis=0)
            else:
                if np.abs(peaks_df.loc[
                              curr["ID"], "Distance to TSS"]) < thresh:
                    curr[
                        "Evidence"] = "Experimental peak; Reference less than %d nt away" % thresh
                    gene_struct = pd.concat(
                        (gene_struct, pd.DataFrame(curr).transpose()),
                        sort=True, axis=0)
                else:
                    # Create two peaks- one for reference and one for experimental
                    curr[
                        "Evidence"] = "Experimental peak; Reference greater than %d nt away " % thresh
                    gene_struct = pd.concat(
                        (gene_struct, pd.DataFrame(curr).transpose()),
                        sort=True, axis=0)

                    # Create the RefSeq one and the name
                    if not is_seen:  # This makes sure that the isoform is not repeated
                        new = curr.copy()
                        new[
                            "Evidence"] = "Reference TSS further than %d nt away of an experimental TSS" % thresh
                        # new["Start"] = curr[]
                        curr_ref = start_anno[
                            (start_anno["gene"] == curr["Gene"]) & (
                                        start_anno["transcript_id"] ==
                                        curr["Transcript"])]
                        new["Start"] = curr_ref["Start"].values[0]
                        new["End"] = curr_ref["End"].values[0]
                        new["Strand"] = curr_ref["Strand"].values[0]
                        new["ID"] = curr_ref.index.values[0]
                        new["is Experimental"] = False
                        new = new.rename("Ref_" + new.name)
                        gene_struct = pd.concat((gene_struct,
                                                 pd.DataFrame(
                                                     new).transpose()),
                                                sort=True, axis=0)
                    is_seen = True
    gene_struct = gene_struct[
        ["Gene", "Transcript", "Evidence", "Chr", "Start", "End", "ID",
         "Strand", "Is Experimental"]]
    return gene_struct


def wrap_create_targets(peaks_f, tss_f, meta_f, bed_f, f_save=None,
                        gois=None, flank_seq=(1000)):

    # meta_f = "Results/output/TSS1.meta"
    # bed_f = "Results/output/TSS1.bed"

    start_anno = pd.read_csv(tss_f, sep="\t", index_col=0)
    start_anno["gene"] = start_anno["gene"].apply(lambda x: x.upper())
    peaks_df = pd.read_csv(peaks_f,sep="\t",index_col=0)
    tss_df = pd.concat((pd.read_csv(meta_f, sep="\t", index_col=0),
                        read_bed_file(bed_f)), axis=1)
    tss_df["Gene"] = tss_df["Gene"].apply(lambda x: x.upper())

    tss_df = tss_df[
        ["Chr", "Start", "End", "ID", "Strand", "Gene", "Transcript",
         "Is Experimental"]].copy()
    tss_df["Evidence"] = ""


    if gois is None:
        gois = set(start_anno["gene"].values)

    gene_struct = pd.DataFrame(columns=tss_df.columns)
    for goi in gois:
        #print(goi)
        gene_struct = pd.concat((gene_struct,info_on_gene(goi,
                                                          peaks_df, tss_df,
                                start_anno, thresh=150)), sort=False,
                                axis=0)
    gene_struct.sort_values(["Gene", "Start","End"], inplace=True)
    gene_struct["ID"] = gene_struct.index

    descriptives = "There are %s silenced secMs" % len(gene_struct['Gene'].unique())
    descriptives = descriptives + "\n" + \
            "There are %s number of expTSSs further than 150 nts away from the reference" % (
            gene_struct.index.str.contains('Ref_')).sum()

    # Extract fasta with 1000bp


    tss_df["Middle"] = int((tss_df["Start"] + tss_df["End"])/2)
    tss_df["Start"] -= flank_seq
    tss_df["End"] += flank_seq

    # Save variables
    if f_save is not None:
        if not f_save.endswith(".tsv"):
            f_save = f_save + ".tsv"
        gene_struct.to_csv(f_save, sep="\t", index=None)

        #Supplemental files
        stats_f = helper.append_name(f_save, prefix="stats",
                                     suffix="txt")
        with open(stats_f,'w') as f:
            f.write(descriptives)
        helper.time_file(os.path.dirname(f_save), name="DATE.txt")
        return
    return gene_struct


def plot_tss_per_isoform(f_save):
    regions = pd.read_csv(f_save, sep="\t")
    ((regions.groupby(
        ["Gene", "Transcript"]).size()).value_counts()).plot.bar()
    plt.title("Number of TSSs per alternative start site (isoform)")
    helper_save(f_save, remove_bg=True)
    return


def wrapper_yaml(yaml_dir):
    for param_f in glob.glob(yaml_dir):
        with open(param_f, 'r') as f:
            doc = yaml.load(f)
        #function = doc["function"]
        #cmd = "python {fnc} {f}"
    return


def main(d):
    # for arg in sys.argv[1:]:
    #     print arg

    files = glob.glob(os.path.join(d, "*.yaml"))
    print(files)

    for param_f in files:
        with open(param_f, 'r') as f:
            doc = yaml.load(f)

            # silenced_f = os.path.join(supplemental,
            #                           "silenced_genes_CHO/DEseq2_HEK_vs_CHO_SecM.tsv")
            #

            goi_f = doc["goi_f"]
            name = param_f.replace(".yaml", "")
            if doc["name"] is not None:
                name = doc["name"]

            meta_f = doc["meta_f"] #""Results/output/TSS1.meta"
            bed_f = doc["bed_f"]  # "Results/output/TSS1.bed"
            peaks_f = doc["peaks_f"]
            ref_tss = doc["ref_tss"]

            outdir = "Results/target_genes"
            if doc["outdir"] is not None:
                outdir = doc["outdir"]
            outdir = os.path.join(outdir, name) #Add name
            gois = pd.read_csv(goi_f, header=None)[0].values
            gois = list(map(lambda x: x.upper(), gois))

            #Make the directory
            helper.mdir(outdir, replace=True, save_old=True)

            #Copy the yaml file into it
            cmd = "cp {param_f} {outdir}/" .format(param_f=param_f,
                                                   outdir=outdir)
            os.system(cmd)

            f_save = os.path.join(outdir, "secM_targets_" + name+".tsv")
            wrap_create_targets(peaks_f, ref_tss, meta_f, bed_f,
                                f_save=f_save, gois=gois)

            plot_tss_per_isoform(f_save)
    return

####


if __name__ == "__main__":
    main(sys.argv[1])

