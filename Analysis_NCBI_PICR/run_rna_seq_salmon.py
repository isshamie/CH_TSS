import sarge
import os
import sys
import glob
from os.path import basename


rna_in_dir = "/data/isshamie/dropbox/RNAseq/2017_12_15_hamster_Seq/"
rna_out_dir = "Results/RNA_align/"
salmon_index = "Results/RNA_align/transcripts_index/"

def salmon(input_files,out_path,salmon_index,thread,lib=''):
    '''run salmon'''
    if lib == '':
        library = '-l A'
    else:
        library = '-l ' + lib
    if len(input_files) == 2:
        cmd=('salmon quant -i {index} {lib} -p {t} -1 {f1} -2 {f2} -o {out_path} --validateMappings').format(
        index=salmon_index,lib=library,t=str(thread),f1=input_files[0],f2=input_files[1],out_path=out_path)
    else:
        cmd=('salmon quant -i {index} {lib} -p {t} -r {f1} -o {out_path}').format(
        index=salmon_index,lib=library,t=str(thread),f1=input_files[0],out_path=out_path)
    print(cmd)
    sarge.run(cmd)

in_f = glob.glob(os.path.join(rna_in_dir,"*_R1_001.fastq.gz"))
for i1 in in_f:
    print(glob.glob(i1))
    i2 = i1.replace("_R1_","_R2_")
    print(glob.glob(i2))
    out_f = os.path.join(rna_out_dir, basename(i1).split("_")[1])
    print(out_f)
    
    salmon([i1,i2],out_f,salmon_index,thread=16,lib='')
