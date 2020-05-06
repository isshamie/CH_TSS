import os


def run_findMotifs(bed_f, out_dir, ref_fa, args=None,  bg=None, len_mo='small'):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if len_mo == 'small':
            cmd = f"nohup findMotifsGenome.pl {bed_f} {ref_fa} {out_dir} -size 200 -len 3,4,5 "
    else:
        cmd = f"nohup findMotifsGenome.pl {bed_f} {ref_fa} {out_dir} -size 200 -len 6,8,10,12     "
    if bg is not None:
        cmd = f"{cmd} -bg {bg} "
    if args is not None:
        cmd = cmd + " ".join(args) + " "
    cmd = f"{cmd} > {out_dir}.log"
    print(cmd)
    os.system(cmd)
    return


# curr_dir = os.path.join(save_dir, "eTSS_bg_rTSS_nocpg")
# run_findMotifs(exp_f_save, out_dir=curr_dir,ref_fa=ref_fa, bg=mrna_f_save)
#
#
# exp_dir = os.path.join(save_dir, "eTSS_motifs_nocpg")
# run_findMotifs(exp_f_save, out_dir=exp_dir,ref_fa=ref_fa)
#
# ref_dir = os.path.join(save_dir, "rTSS_motifs_nocpg")
# run_findMotifs(mrna_f_save, out_dir=ref_dir,ref_fa=ref_fa)
#
# hg38_out = os.path.join(save_dir, "hg38_nocpg")
# hg38_genome = "/data/isshamie/genome/hg38/GCF_000001405.38_GRCh38.p12_genomic.fna"
# run_findMotifs(hg38_f_save, out_dir=hg38_out,ref_fa=hg38_genome)
#