import os
from IPython.display import Image


def run_heatmap(tss_bed,bigwig,outfile,dist=2000,numProc=6,
                dispImage = False,computeMat_args = '',plotHeat_args = '',
               is_tss_f = ''):
    
    
    cmd = 'computeMatrix reference-point -S {bw} -R {bd} -a {d} -b {d} -out {out} --numberOfProcessors {p} {args}'.format(
        bw=bigwig, bd = tss_bed,d=dist,out=outfile,p=numProc,args=computeMat_args)
    print(cmd)
    os.system(cmd)
    #!{cmd}   
    
    cmd = """plotHeatmap\
     -m {out}\
     -out {out}.png \
     --heatmapHeight 15  \
     --refPointLabel TSS.center \
     --regionsLabel TSS \
     --plotTitle 'ATAC signal' {args}""".format(out=outfile,args = plotHeat_args)

    print(cmd)
    os.system(cmd)
    #!{cmd}
    
    
    if dispImage:
        return Image(filename=outfile + '.png')
    return
