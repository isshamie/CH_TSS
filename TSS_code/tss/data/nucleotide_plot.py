from tss.visualize.fig_utils import wrap_plots
from tss.utils.Homer import homer_nucleotide
import matplotlib.pyplot as plt


f1 = "Results/output/TSS1.exp.bed"
f2 = "Results/Figures/Figure3/alt/A.TSS1_mrna"

f = plt.figure()
wrap_plots([f1,f2],homer_nucleotide,"Results/Figures/Figure3/C.combine",["Experimental","RefSeq"],ref_fa,1000,(0.1,0.65),False)

f1 = "Results/output/TSS1.exp.bed"
f2 = "Results/Figures/Figure3/alt/A.TSS1_mrna"

f = plt.figure()
wrap_plots([f1,f2],homer_nucleotide,"Results/Figures/Figure3/C_200.combine",["Experimental","RefSeq"],ref_fa,200,(0.1,0.65),False)