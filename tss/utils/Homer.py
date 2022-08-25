# coding: utf8
import sarge, sys
import os
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')
from tss.visualize.plot_tss_results import *
from tss.data.data_io import *
import pyupset as pyu
from tss.config import HOMER_PATH
print(HOMER_PATH)

#Append HOMERS path variable
#os.environ["PATH"] = HOMER_PATH + ":" + os.environ["PATH"]

#####
#See Homer.ucsd.edu for more info on all these functions
#####


########################################
def find_peaks(tag_dir, out_file, peak_style, control_dir,
               otherParams=('')):
    '''
    find peaks with homer from the tag directory.
    * tag_dir: Created from makeTagDirectory
    * out_file: Output tsv file
    * peak_style: 'tss';'histone'; 'dnase' ..
    * control_dir: The background tag directory
    * otherParams: Additional parameters such as how much fold increase over the background to be called a peak 
    '''
    cmd = 'findPeaks {tag} -style {style} -o {out} -i {control} -norm 1000000 {other}'.format(
                tag=tag_dir,style=peak_style,out=out_file,control=control_dir,
                other=' '.join(otherParams))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    return


########################################   
def merge_peaks(input_files, output_file, dist='given', type_merge=''):
    """
    * input_files: a list of peak files, name format is 5gro_and_gro
    * output_file: final merged peak file
    # dist: distance between peaks to merge. Default is given, which means overlap
    """
    # Output file name add merge
    if not output_file.endswith('.merge'):
        output_file = output_file + '.merge'
    out_venn = output_file + '_venn'
    out_matrix = output_file + '_mat'
    
    # If only one file, do not want to run bc wont have proper layout, so just duplicate the file and merge
    is_one = False
    if len(input_files) == 1:
        is_one = True
        tmp = '%s_copy' % (input_files[0])
        cmd = 'cp {input} {copy}'.format(input=input_files[0], copy=tmp)
        input_files.append(tmp)
        print('Only one file to merge. Creating copy of original.')
        print(cmd);sys.stdout.flush()
        sarge.run(cmd)

    #subprocess.call("source /home/isshamie/.bashrc")
    # Run Homer command
    cmd = ('mergePeaks -d {dist} -venn {out_venn} '
           '-matrix {out_matrix} -strand {in_files} > {out}').format(
           dist=str(dist), in_files=' '.join(input_files), out_venn=out_venn, out_matrix=out_matrix, out=output_file)

    print(cmd);sys.stdout.flush()
    #sarge.run(cmd)
    os.system(cmd)

    # remove the copied peak file
    if is_one:
        cmd = 'rm {copy}'.format(copy=tmp)
        print(cmd);sys.stdout.flush()
        sarge.run(cmd)

    # Load in mergePeaks and change the column names.
    peaks = pd.read_csv(output_file, sep='\t')
    col_names = list(peaks.columns.values)
    meta_info = col_names[0] #the parameters are written as the first column name. Want to save it as a separate file
    col_names[0] = 'ID'
    col_names[1] = 'Chr'
    col_names[2] = 'Start'
    col_names[3] = 'End'
    col_names[4] = 'Strand'
    peaks.columns = col_names

    # Remove the duplicate column if its one
    if is_one:
        peaks = peaks.loc[:, ~(peaks.columns.str.contains('_copy'))]
        peaks['Parent files'] = peaks.columns.values[8]  # Just the file name

    # Save meta file
    with open(output_file + '.meta', 'w') as f:
        f.write(meta_info)
    f.close()

    # Change ID to the index, which is unique
    peaks['ID'] = peaks.index.astype(int)

    # If all, reduce the set to only peaks that are present in all the files.
    # Save that to the file, and save the non-exclusive ones to output_file_raw
    if type_merge == 'all' and 'Parent files' in peaks.columns.values:
        peaks.to_csv(output_file+'_raw', sep='\t', index=False)
        num_breaks = max(peaks['Parent files'].apply(lambda x: x.count('|')))
        filt = '^([^\|]*\|){%d}[^\|]*$' % (num_breaks)  # regex to have the number of files
        peaks = peaks[peaks['Parent files'].str.contains(filt)]
    # Save file
    peaks.to_csv(output_file, sep='\t', index=False)  # Keep the index, which is unique
    return


########################################
def annotate_peaks(peak_file, output_file, ref_fa, ref_anno):
    """
    This function annotates peaks in terms of the annotation (e.g. promoter, exon..).
    Also marks closest TSS to each peak as well. 
    """
    if ref_anno.endswith('gtf'):
        anno = '-gtf ' + ref_anno
    else:
        anno = '-gff ' + ref_anno
    cmd = 'annotatePeaks.pl {peaks} {genome} {annotation} > {out}'.format(
        peaks=peak_file, genome=ref_fa,
        annotation=anno, out=output_file)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    return


########################################
def annotate_filter(anno_file,marks=('promoter','exon','Intergenic')):
    '''Separates annotated peaks into different files given the mark it is annotated by.'''
    anno = pd.read_csv(anno_file,sep='\t')
    anno.dropna(axis=1,inplace=True,how='all')
    anno.fillna('',inplace=True)

    for mark in marks:
        tmp = anno[anno['Annotation'].astype('str').str.match(mark)]
        tmp.to_csv(anno_file+'_' + mark,index=None,sep='\t')

    return

########################################
# Not needed
# def rm_5GRO_ctrl(GRO5_tag,ctrl_tag):
#     '''this function remove ctrl_tag coverage from GRO-Cap tag
#     '''
#     before_sub = GRO5_tag+'/genome.tages_before_sub.tsv'
#     after_sub = GRO5_tag+'/genome.tags.tsv'
#     os.rename(after_sub,before_sub)
#     df_raw = pd.read_csv(before_sub,sep='\t',header=None)
#     df_ctr = pd.read_csv(ctrl_tag+'/genome.tags.tsv',sep='\t',header=None)
    

######################################## 
def hist(tag_dir,hist_out,ref_fa,anno,mode='peak',peak='',region=2000,res=10,pc=0,hist_norm = 100):
    '''this function gets tag coverage around tss. Creates histogram, normalized histogram, and also stores the raw coverage for each read in a matrix form
    * tag_dir: tag directory
    * hist_out: Base filename to save to. Will save additional files with different extensions. All are tsv files
    * ref_fa: Reference assembly
    * anno: gff file or gtf file
    * mode: 'peak', 'tss'. Peak uses peak centers for histogram, tss uses the annotated start sites as center.
    * pc: number of tags to consider at each position
    * region: length to conisder in the x axis. e.g. 4000 means -2000 to 2000 around tss
    * res: resolution of the histogram in basepairs.
    * pc: Max number of tags for each bp. 0 means no max
    * hist_norm: In normalized histogram, this number is min tag count to normalize by for each read
    '''
    if anno.endswith('gtf'):
        anno = '-gtf ' + anno
    else:
        anno = '-gff ' + anno
    
    if mode == 'tss':
        peak = 'tss'
    #elif mode == 'peak':
    if peak == '':
        raise ValueError('input is empty')
    
    cmd = ('annotatePeaks.pl {peak} {ref_fa} {anno} -fragLength 1 -size {size} -hist {bin} -d {dir} -pc {pc} > {out}').format(
            peak=peak,ref_fa=ref_fa,anno=anno, size=str(region),
            bin=str(res),dir=tag_dir,pc=str(pc), out=hist_out)
    print(cmd)
    sarge.run(cmd)
    
    cmd1 = ('annotatePeaks.pl {peak} {ref_fa} {anno} -fragLength 1 -size {size} -hist {bin} -histNorm {hist_norm} -d {dir} -pc {pc} > {out}').format(
            peak=peak,ref_fa=ref_fa,anno=anno,size=str(region),bin=str(res),dir=tag_dir,pc=str(pc),hist_norm=str(hist_norm),out=hist_out+'Norm')
    print(cmd1)
    sarge.run(cmd1)

    cmd2 = ('annotatePeaks.pl {peak} {ref_fa} {anno} -fragLength 1 -size {size} -hist {bin} -ghist -d {dir} -pc {pc} -strand + > {out}').format(
                    peak=peak,ref_fa=ref_fa,anno=anno,size=str(region),bin=str(res),dir=tag_dir,pc=str(pc),out=hist_out + 'MatS')
    print(cmd2)
    sarge.run(cmd2) 

    cmd3 = ('annotatePeaks.pl {peak} {ref_fa} {anno} -fragLength 1 -size {size} -hist {bin} -ghist -d {dir} -pc {pc} -strand - > {out}').format(
            peak=peak,ref_fa=ref_fa,anno=anno,size=str(region),bin=str(res),dir=tag_dir,pc=str(pc),out=hist_out + 'MatAS')
    print(cmd3)
    sarge.run(cmd3)
    return


# def annotate_peaks_shift_center(peak_file,output_file,ref_fa,annotation,col_center='Distance to TSS',no_col=False,center_id='NCBI'):
#     '-center <motif file> (This will re-center peaks on the specified motif, or remove peak'


def filter_anno_peak(in_peak_file,filter_peak_file):
    '''this function extracts the reliable TSS from the peak file
    The rule is: for each 5GRO, get overlap of peaks against different GROseq.
    Then get union set from pevious peaks.
    '''
    df = pd.read_csv(in_peak_file,sep='\t',header=0)
    gro5 = []
    gro = []
    for f in df['Focus Ratio/Region Size']:
        files = f.split('|')
        for sub_f in files: # sub_f is peak file result
            peaks = sub_f.split('_and_')  # peaks has 5gro file and gro file name
            for p in peaks:
                if p in gro5 or p in gro:
                    continue
                else:
                    if '5GRO' in p:
                        gro5.append(p)
                    else:
                        gro.append(p)
    print((gro5,gro))

    def extract_peak(gro5,gro,fns):
        '''fns is the splited filename in 6th column of annotated peak file'''
        keep = []
        res = True
        for g5 in gro5:
            keep.append([g5+'_and_'+g in fns for g in gro])
        for k in keep:
            if False in k:
                res=False
        return res
        
    # filter out the annopeak
    cri = df['Focus Ratio/Region Size'].map(lambda x: extract_peak(gro5,gro,x))
    df = df[cri]
    df.to_csv(filter_peak_file,sep='\t',index=False)
    return


########################################
def createPeakFileFromGFF(annotation_file, output_file,
                          anno_of_interest='mRNA', is_start=True,
                          shift=0):
    """
    Function to create a peak file based on GFF and specific annotation of interest
    
    Input
    * annotation_file: gff3 file. Not gtf!
    * output_file: File to save to
    * anno_of_interest: a text term that has to be in the 'Annotation' column of the annotation file
    * is_start: The start site is the center of the peak

    Creates a peak file of the Homer format and saves it to output_file

    """
    # Read in annotation file
    if '.gtf' in annotation_file or '.gtf.gz' in annotation_file:
        print("Only works for gff files, not gtf. Convert using gffread from the command line.")
    genome_ann = pd.read_csv(annotation_file, comment='#', sep='\t',
                             header=None)
    col_names = list(genome_ann.columns.values)
    col_names[0] = 'Chr'
    col_names[1] = 'How'
    col_names[2] = 'Annotation'
    col_names[3] = 'Start'
    col_names[4] = 'End'
    col_names[6] = 'Strand'
    col_names[8] = 'ID'

    genome_ann.columns = col_names
    genome_ann.sort_values(['Chr', 'Start', 'End'], inplace=True)

    curr = genome_ann[genome_ann['Annotation'] == anno_of_interest]
    curr = curr[['ID', 'Chr', 'Start', 'End', 'Strand']]
    # Make the center of the peak the start site with flanking base pairs
    if is_start:
        # curr['End'] = curr['Start'] + 1
        # curr['Start'] -= 1
        curr2 = curr.copy()


        curr['actual_start'] = curr.apply(
            lambda x: x['Start'] if x['Strand'] == '+' else x['End'],
            axis=1)
        ## If the strand is positive, then the end moves to the start, but if its negative, the start moves to the end.
        curr['End'] = curr2.apply(
            lambda x: x['Start'] + shift if x['Strand'] == '+' else x[
                                                'End'] + shift, axis=1)
        curr['Start'] = curr2.apply(
            lambda x: x['Start'] - shift if x['Strand'] == '+' else x[
                                                'End'] - shift, axis=1)
    curr.to_csv(output_file, index=None, sep='\t')
    return


########################################
#parseGTF.pl with Homer
def run_parsegtf(annotation, output):
    cmd = f"parseGTF.pl {annotation} tss > {output}"
    print(cmd)
    os.system(cmd)
    return


########################################
def peakFileToPeakFile(desired_peak_file, tag_peak_file,
                       distance=1000, is_save=1, f_save=None,
                       is_peak=True, is_bed=False):
    """
    Function to filter the peaks in a file to only the ones that are within a certain distance
    of at least one peak in another file. 

    Input:
    * desired_peak_file: File which contains the peaks to be filtered
    * tag_peak_file: File that contains peaks to see if nearby the desired peaks
    * distance: Maximum distance to consider between peaks
    * is_save: To save the file or not. If on, will save peaks to file desired_peak_file+'filt'

    Output:
    updated_desired_peaks: The filtered peaks.
    """

    desired_peaks = pd.read_csv(desired_peak_file, sep='\t')
    if is_peak: #Either peak or merge file
        tag_peaks = read_peak_file(tag_peak_file)
    elif is_bed:
        tag_peaks = read_bed_file(tag_peak_file)
    else:
        tag_peaks = pd.read_csv(tag_peak_file, sep='\t')
    if not 'actual_start' in tag_peaks.columns.values:
        tag_peaks['actual_start'] = tag_peaks.apply(lambda x: x['Start'] if x['Strand'] == '+' else x['End'], axis=1)
    if not 'actual_start' in desired_peaks.columns.values:
        desired_peaks['actual_start'] = desired_peaks.apply(lambda x: x['Start'] if x['Strand'] == '+' else x['End'],
                                                            axis=1)

    inds_to_keep = []
    for ind, val in desired_peaks.iterrows():
        curr = tag_peaks[tag_peaks['Chr'] == val['Chr']]
        diff = val['actual_start'] - curr['actual_start']
        if np.sum(np.abs(diff) < distance) > 0:
            inds_to_keep.append(ind)

    updated_desired_peaks = desired_peaks.loc[inds_to_keep]
    
    if is_save:
        if f_save is None:
            updated_desired_peaks.to_csv(desired_peak_file + 'filt', index=None, sep='\t')
        else:
            updated_desired_peaks.to_csv(f_save + 'filt', index=None, sep='\t')
    return updated_desired_peaks


########################################
def convert_peak_to_bed_file(peak_file,out_name):
    ''' Convert homer peak tsv file to a bed file ''' 
    cmd = 'pos2bed.pl %s > %s' % (peak_file,out_name)
    print(cmd)
    sarge.run(cmd)
    return

########################################
def convert_intersection_HOMER_to_pyupset(df,f_name=''):
    ''' Takes in merge tsv file and creates a pyupset figure that shows the overlaps of the different files'''

    vals = df.iloc[-1,-1].split('|')
    vals = df.columns.values[:-2] #-2 so Dont include Name and Total
    df_out = dict()
    for i in vals:
        df_out[i] = [] 

    count = 0
    for i in range(1,len(vals)+1):
        for j,val in df.iterrows():
            #print(j)
            curr_count = val['Total']
            #curr_count = int(df[df['Name']  == '|'.join(j)]['Total'])

            samples = df.columns.values[~val.isnull()][:-2] #-2 so Dont include Name and Total
            for n in samples:
                #print(n)
                #print(count)
                df_out[n].extend(list(range(count,count+curr_count)))

            count += curr_count
    for i in df_out:
        df_out[i] = pd.DataFrame(df_out[i],columns=['Name'])#pd.DataFrame(df_out[i])
   
    if not f_name == '':
        pyu.plot(df_out)
        plt.savefig(f_name)
    return df_out


########################################
def homer_trim_cmd(input_file,seq='AGATCGGAAGAGCACACGTCT'):
    ''' Trims the fastq.gz file based on adaptor sequence '''
    cmd = ('{HOMER_PATH}/homerTools trim -3 {seq} -mis 2 '
           '-minMatchLength 4 -min 20  {in_file}').format(
        HOMER_PATH=HOMER_PATH,seq=seq, in_file=input_file)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    # gzip and save as similar to trimmomatic, trim_{input}
    cmd = 'gzip %s.trimmed' % (input_file) #CHBrainnegmaybe7neg1_mSTARTinput_JHS1082_SD_CACGAT_S104_L005_R1_001.fastq.gz.trimmed'
    print(cmd);sys.stdout.flush(); 
    sarge.run(cmd)
    cmd = 'mv %s.trimmed.gz trim_%s' % (input_file,input_file)
    print(cmd);sys.stdout.flush(); 
    sarge.run(cmd)
    return


#######################################
def homer_nucleotide(input_file, output_file, ref_fa, size=1000,
                     desired_lims=None, only_plot=False, f=None):
    """
    Input
    * input_file: peak file
    * output_file where to save the nucleotide frequency plot
    * ref_fa: Reference annotation

    Will also plot the nuc freqs and save to {output_file}.png
    """

    # cmd = 'annotatePeaks.pl f06_annoPeaks/merge_bg_2.anno_promoter /data/genome/hamster/picr/picr.fa -size 1000 -hist 1 -di > nuc_freq.txt'
    if not only_plot:
        cmd = 'annotatePeaks.pl %s %s -size %s -hist 1 -di > %s' % (
                input_file, ref_fa, size, output_file)
        print(cmd)
        sarge.run(cmd)
        tmp = pd.read_csv(output_file, sep='\t', index_col=0)
    else:  # Alrady have the frequency plots, just want to plot
        tmp = pd.read_csv(input_file, sep='\t', index_col=0)
    if f is None:
        f = plt.figure(dpi=300)
        ax = f.add_subplot(111)
        to_return = False
    else:
        ax = f.gca()
        to_return = True
    ax.plot(tmp.index.values, tmp['A frequency'], linewidth=1)
    ax.plot(tmp.index.values, tmp['C frequency'], linewidth=1)
    ax.plot(tmp.index.values, tmp['T frequency'], linewidth=1)
    ax.plot(tmp.index.values, tmp['G frequency'], linewidth=1)
    if desired_lims is None:
        ax.vlines(0, ax.get_ylim()[0], ax.get_ylim()[1])
    else:
        ax.vlines(0, desired_lims[0], desired_lims[1])
    ax.legend()
    ax.set_xlabel('bp from TSS')
    ax.set_ylabel('Frequency')

    helper_save(output_file)

    if to_return:
        return f
    else:
        return
    #plt.savefig(output_file + '.png')


#######################################
def extract_peak_sequences(bed_file,peak_list,genome,f_save,upstream=1000,downstream=100):
    ''' 
    Extracts genomic sequences given a list of peaks
    Assumes ID is the 4th column
    '''
    #Save the list to bed_file first. Note bed file is index base 0!
    all_peaks = pd.read_csv(bed_file,sep='\t',header=None)
    peaks = all_peaks[all_peaks[3].isin(peak_list)] #The ID is the 4th column
    peaks.to_csv(f_save + 'before_extension.bed',sep='\t', index=False, header=None)
    for ind,val in peaks.iterrows(): 
        center = int((val[1]+val[2])/2)
        if val[5] == '+':
            peaks.set_value(ind, 1, center - upstream + 1)#-1 for bed
            # index
            peaks.set_value(ind, 2, center + downstream + 1)
#             val[1] = val[1] - downstream
#             val[2] = val[2] + upstream
        else:
            peaks.set_value(ind, 1, center-downstream + 1)
            peaks.set_value(ind, 2, center+upstream + 1)

    peaks.to_csv(f_save + '.bed',sep='\t', index=False, header=None)
    #Get sequence info for each peak
    cmd = 'homerTools extract %s %s -fa > %s' % (f_save + '.bed',
                                                 genome, f_save)
    print(cmd)
    os.system(cmd)    
    
    #Get sequence info for each peak
    cmd = 'homerTools extract %s %s -fa > %s' % (f_save +'before_extension.bed', genome,'before_extension_' + f_save)
    print(cmd)
    os.system(cmd)   
    return


#######################################
def query_gene(peak_files,gene):
    """Takes in a list of peak files, looks for peaks near specific gene and rerturns those peaks
    Input:
    * peak_files: List of peak tsv files that have been annotated
    * gene: gene name 
    """

    gene_peaks = pd.DataFrame()
    count = 0
    query = '(?i){g}'.format(g=gene)
    for f in peak_files:
        tmp = pd.read_csv(f,sep='\t',comment='#',index_col=0)
        tmp.dropna(axis=1,how='all',inplace=True)
        tmp.dropna(inplace=True)
        #tmp[tmp['Annotation'].index.str.contains(gene_id)]
        if count == 0:
            gene_peaks = tmp[tmp['Nearest PromoterID'].str.contains(query)]
        else:   
            gene_peaks = gene_peaks.append(tmp[tmp['Nearest PromoterID'].str.contains(gene)])
        count += 1
        print((len(gene_peaks),f))
    return gene_peaks


# if __name__ == '__main__':
#     in_peak_file = '/data/shangzhong/TSS/fq/f05_annoPeaks/merge.anno'
#     filter_peak_file = '/data/shangzhong/TSS/fq/f05_annoPeaks/merge_filter.anno'
#     filter_anno_peak(in_peak_file,filter_peak_file)
                

