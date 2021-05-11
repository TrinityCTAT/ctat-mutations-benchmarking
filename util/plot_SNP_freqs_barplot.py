#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd
import matplotlib.backends.backend_pdf
import collections
from collections import Counter, OrderedDict
import sys, os
import glob
import ntpath


## Count snps belonging to each category (TP, FP, FN)
def count_snps(df,snps,param):
    data = dict.fromkeys(snps, 0)
    if param == 'FN':
        vals = df.groupby('Ref_SNP').Ref_SNP.count().to_dict()
        data.update(vals)
    elif param == 'TP' or param == 'FP':
        vals = df.groupby('RNA_SNP').RNA_SNP.count().to_dict()
        data.update(vals)
    return data

## box plot per chr
def plot_bars(lsts,ax):
    for lst in lsts:
        h = lst.get_height()
        ax.text(lst.get_x()+lst.get_width()/2, 1.05*h, '%d'%int(h),
                ha='center', va='bottom', rotation=90, size='smaller')

    return


## SNP frequency 
def snp_freq(file_names, output_dir):
    
    param_lst = ['TP', 'FP', 'FN']

    dfile_names = []
    #file_names = [os.path.join(output_dir, "class_"+ntpath.basename(name)) for name in file_names]
    for i in range(len(file_names)):
        globals()['df%s' % i] = pd.read_csv(file_names[i],sep='\t',index_col=0)
        dfile_names.append(globals()['df%s' % i])

    num_chr = len(np.asarray(list(df0.index)))#.get_values())
    chr_name = np.asarray(df0.index)#.get_values()

    snps = ['A:G', 'T:C', 'G:C', 'C:A', 'G:A', 'C:T', 'T:G', 'C:G', 'G:T', 'A:T', 'A:C', 'T:A']#snps of interest
    
    #convert snps to a dictionary with all values assigned to 0 count
    N = len(snps)
    index = np.arange(N)  # the x locations for the groups
    width = 0.2      # the width of the bars
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    all_data = {}
    lists = ()
    for i in range(len(dfile_names)):
        final_data = collections.Counter()
        for param in param_lst:
            df_subset = dfile_names[i][dfile_names[i].Class==param]#subset the dataframe to include on the param values
            data = count_snps(df_subset,snps,param)
            final_data += Counter(data)

        ### Select SNPs only in snps - INDELs excluded
        fd = {snp:final_data[snp] for snp in snps}
        vals = list(OrderedDict(fd).values())
        labels = tuple(OrderedDict(fd).keys())

        globals()['lists%s' % i] = ax.bar(index+width*i, vals, width)
        lists = lists + (globals()['lists%s' % i],)
    ax.legend(lists,tuple(file_names))
    ax.autoscale(tight=True)
    ax.set_xticks(index+width)
    ax.set_xticklabels(labels,rotation=45)
    for lst in lists:
        plot_bars(lst,ax)
    plt.ylabel('Counts')
    plt.title("SNPs")
    plt.savefig(os.path.join(output_dir, 'plot_snp_freq.png'))
    
    return


def main():

    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs CTAT_mutations benchmarking\n")

    parser.add_argument('--input_scored_table', nargs='+', required = True,
                        help="List of input files containing class column with: TP, FP, FN assignments")

    parser.add_argument("--output_dir", required=True, help="output directory")

    args = parser.parse_args()

    #create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)


    input_filenames = args.input_scored_table
    
    snp_freq(input_filenames, args.output_dir)
    
    
    sys.exit(0)



if __name__ == "__main__":
    
    main()


