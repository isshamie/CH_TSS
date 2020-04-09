import sys
import yaml
from os.path import join, basename
import os
import glob

params_f = "parameters/params.yaml"
with open(params_f,'r') as f:
    params = yaml.load(f)
print("here", params)

QC_out = join("Results","QC")
if not os.path.exists(QC_out):
	os.mkdir(QC_out)

for i in glob.glob(params['RNAseq_folder']+"/*"):
	print(i)
	cmd = f"fastqc {i} -o {QC_out}" #{join(QC_out, basename(i.replace('.gz','')))}"
	print(cmd)
	os.system(cmd)

outdir = join("Results", "multiqc")
if not os.path.exists(outdir):
    os.mkdir(outdir)

cmd = f"multiqc {QC_out} --outdir {outdir}"
print(cmd)
os.system(cmd)

