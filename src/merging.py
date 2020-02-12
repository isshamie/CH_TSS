from Homer import *
import glob
import json

all_peak_merge_files = glob.glob("f04_peaks/*")
print(all_peak_merge_files)

def merge_peaks(p):
    """
    @param p: The parameters dictionary
    :param outdir:
    :param peaks_dir:
    :param dist:
    :return:
    """

    all_peak_merge_files = glob.glob(p["peaks_dir"])
    if not os.path.exists(p["outdir"]):
        os.mkdir(p["outdir"])
    peak_f = os.path.join(p["outdir"],'samples.merge')


    # Run merging
    merge_peaks(all_peak_merge_files, peak_f, dist=p["dist"])


    anno_f = peak_f + ".anno"
    annotate_peaks(peak_file=peak_f, output_file= anno_f,
                   ref_fa=p["genome"], annotation=p["annotation"])


    seq_f = anno_f + ".fa"
    cmd = 'homerTools extract {f_in} {gen} -fa > {f_out}'.format(
        f_in=anno_f, gen=p["genome"], f_out=seq_f)
    print(cmd)
    os.system(cmd)


    # Save the parameters
    with open('params_used.json', 'w') as fp:
        json.dump(p, fp)
    return


def merge_tags(p):
    ### Parameters setup
    # Create output directory
    """
    :param outdir:
    :param data_folder:
    :param tissues:
    """
    if p["outdir"]is None:
        p["outdir"] = "./data/"
    if p["data_folder"] is None:
        p["data_folder"] = "../data/"

    ### Till here

    ## Cap Tags
    all_tags_start = []
    for t in p["tissues"]:
        for s in ['GROCap', 'START']:
            # curr = glob.glob(j + '*f03_tags/trim*.peak')
            curr_raw = glob.glob(p["data_folder"] + t + '/*/')
            # print(curr_raw)
            for j in curr_raw:
                curr_type = j.split('/')[-2]
                if curr_type == 'GROCap':

                    new = glob.glob(os.path.join(p["data_folder"], s, j,
                                                 'f03_tags/*GROCap*'))
                elif curr_type == 'START':
                    new = np.array(glob.glob(
                        os.path.join(p["data_folder"], s, j, 'f03_tags/*')))
                    if not len(new) == 0:
                        new = (new[(np.array(
                            map(lambda x: 'input' not in x,
                                new)))])  # Remove the input
                else:
                    continue
                all_tags_start.extend(new)

    curr_out = os.path.join(p["outdir"],'tags_TSS_merged')

    cmd = 'makeTagDirectory {out_f} -single -d {files}'.format(
        out_f=curr_out, files=' '.join(all_tags_start))
    print(cmd)
    os.system(cmd)

    ## Input Tags
    all_tags_input = []
    for t in p["tissues"]:
        for s in ['GROCap','START']:
            #curr = glob.glob(j + '*f03_tags/trim*.peak')
            j = glob.glob(os.path.join(p["data_folder"], t,s))
            if len(j) > 0:
                j = j[0]
            else:
                continue
            curr_type = s
            if curr_type == 'GROCap':
                #print(os.path.join(j,'f03_tags/*5GRO*'))
                new = np.array(glob.glob(os.path.join(p["data_folder"],s,j,'f03_tags/*GRO*')))
                if not len(new) == 0:
                    new = (new[(np.array(map(lambda x: 'Cap' not in x,new)))]) #Remove the input
            elif curr_type == 'START':
                new = np.array(glob.glob(os.path.join(p["data_folder"],s,j,'f03_tags/*input*')))
            else:
                continue
            all_tags_input.extend(new)

    curr_out = os.path.join(p["outdir"],'tags_input_merged')
    cmd = 'makeTagDirectory {out_f} -single -d {files}'.format(out_f=curr_out,files = ' '.join(all_tags_input))
    print(cmd)
    os.system(cmd)
    return


def merging(p):
    if p["outdir"]is None:
        p["outdir"] = "./data/"
    if p["data_folder"] is None:
        p["data_folder"] = "../data/"
    if p["peaks_dir"] is None:
        p["peaks_dir"] = "../data/peaks/"
    if p["dist"] is None:
        p["dist"] = "given"

    merge_peaks(p)
    merge_tags(p)
    return
