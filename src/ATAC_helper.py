import os
from IPython.display import Image
import textwrap
import inspect

def run_heatmap(tss_bed,bigwig,outfile,dist=2000,numProc=6,
                dispImage = False,computeMat_args = '',plotHeat_args = '',
               is_tss_f = '', to_compute=True, label=None):    
    if to_compute:
        cmd = 'computeMatrix reference-point -S {bw} -R {bd} -a {d} -b {d} -out {out} --numberOfProcessors {p} {args}'.format(
            bw=bigwig, bd = tss_bed,d=dist,out=outfile,p=numProc,args=computeMat_args)
        print(cmd)
        os.system(cmd)
    #!{cmd}   
    if label is None:
        label=outfile
    
    cmd = """plotHeatmap\
     -m {out}\
     -out {out}.png \
     --heatmapHeight 15 \
     --refPointLabel "TSS Distance (bp)" \
     --legendLocation none \
     --xAxisLabel TSS Distance (bp) \
     --samplesLabel "{label}" \
     --plotTitle 'ATAC signal' {args}""".format(out=outfile,args = plotHeat_args, label=label)

    cmd = """plotHeatmap\
     -m {out} \
     -out {out}.png \
     --heatmapHeight 15  \
     --legendLocation none \
     --xAxisLabel "TSS Distance (bp)" \
     --regionsLabel TSS \
     --samplesLabel "{label}" \
     {args}""".format(out=outfile,args = plotHeat_args, label=label)

    cmd = inspect.cleandoc(cmd)
    print(cmd)
    os.system(cmd)
    #!{cmd}
    
    
    if dispImage:
        return Image(filename=outfile + '.png')
    return
