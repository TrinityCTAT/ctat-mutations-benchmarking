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


## plots the ROC/Precision Recall curves based on depth values (compares methods)
def comparative_plots_per_depth(file_names, output_dir):
    
     pdf1 = matplotlib.backends.backend_pdf.PdfPages(os.path.sep.join([output_dir, "Precision_recall_depth.pdf"]))
     pdf2 = matplotlib.backends.backend_pdf.PdfPages(os.path.sep.join([output_dir, "ROC_depth.pdf"]))
     
     for f_name in file_names:
         df = pd.read_csv(f_name, sep='\t')
         df_group = df.groupby('ts_min_rna_cov', as_index=False)

         for name,group in df_group:
             sub_df = group
             #print(group)
             sn = sub_df['sn']
             ppv = sub_df['ppv']
             tp = sub_df['tp']
             fp = sub_df['fp']
             fn = sub_df['fn']
             min_cov = list(sub_df['eval_min_rna_cov'])
             point_labels = min_cov

             ############
             ## figure 1 - Precision_recall_depth.pdf

             fig1 = plt.figure(1)
                 
             plt.plot(sn,ppv,label=name, marker='.')
             for i, label in enumerate(point_labels):
                 lst_sn = list(sn)
                 lst_ppv = list(ppv)
                 plt.text(lst_sn[i], lst_ppv[i], label, fontsize=8,color = 'k',fontweight = 'bold' ) 

             fdr = 1-ppv
             plt.xlabel('sn')
             plt.ylabel('ppv')
             plt.title(os.path.basename(f_name))
             plt.legend()

             ###########
             ## figure 2 - ROC_depth.pdf

             fig2 = plt.figure(2)
             tn = 3.2E9-(tp+fn+fp)
             fpr = fp / (fp+tn)
             plt.plot(fpr,sn,label=name,marker='.')
             plt.xticks(fontsize=7.5)
             for i, label in enumerate(point_labels):
                 lst_fpr = list(fpr)
                 lst_sn = list(sn)
                 plt.text(lst_fpr[i], lst_sn[i], label, fontsize=8,color = 'k',fontweight = 'bold' ) 

             plt.xlabel('fpr')
             plt.ylabel('sn')
             plt.title(os.path.basename(f_name))
             plt.legend()
         
         pdf1.savefig(fig1)
         pdf2.savefig(fig2)
         plt.close(fig1)
         plt.close(fig2)
     
     pdf1.close()
     pdf2.close()
     
     return

def make_legend_string(longname):

    substr_size = 3

    newname_pts = list()
    for pts in longname.split("."):
        if len(pts) > substr_size:
            pts = pts[0:substr_size]
            
        newname_pts.append(pts)
    
    longname = ".".join(newname_pts)
    #if len(longname) > 25:
        # take first 7 and last 18
    #    longname = longname[0:7] + ".." + longname[18:]
    

    return longname


## plots the ROC/Precision Recall curves based on the prediction algorithm (compares depth)
def comparative_plots_per_method(file_names, output_dir):
    
    pdf1 = matplotlib.backends.backend_pdf.PdfPages(os.path.join(output_dir, "Precision_recall_method.pdf"))
    pdf2 = matplotlib.backends.backend_pdf.PdfPages(os.path.join(output_dir, "ROC_method.pdf"))
    outdata_filename = os.path.join(output_dir, "accuracy_data_summary.tsv")

    ## read in all data, create single dataframe
    DF = pd.DataFrame()
    for f_name in file_names:
        df = pd.read_csv(f_name,sep='\t')
        df['f_name'] = f_name
        DF = pd.concat([DF,df])

    DF.to_csv(outdata_filename, sep="\t")
    
    ## generate plots
    df_group = DF.groupby('ts_min_rna_cov', as_index=False)
    for name,group in df_group:
        sub_df = group.groupby('f_name', as_index=False)
        for rna_cov, df in sub_df:
            f_name = df.f_name.unique()[0]
            sn = df['sn']
            ppv = df['ppv']
            tp = df['tp']
            fp = df['fp']
            fn = df['fn']
            min_cov = list(df['eval_min_rna_cov'])
            #print(df)

            ###############
            ## plot - Precision_recall_method.pdf
            fig1 = plt.figure(1)
            point_labels = min_cov
            legend_label = make_legend_string(os.path.basename(f_name))
            
            #if len(legend_label) > 35:
            #    legend_label = legend_label[:35]
            

            plt.plot(sn,ppv,label=legend_label, marker='.')
            for i, label in enumerate(point_labels):
                lst_sn = list(sn)
                lst_ppv = list(ppv)
                plt.text(lst_sn[i], lst_ppv[i], label, fontsize=8,color = 'k',fontweight = 'bold' ) 

            fdr = 1-ppv
            plt.xlabel('sn')
            plt.ylabel('ppv')
            plt.title("rna_cov="+str(name))
            #plt.ylim(ymax=1)
            plt.legend()

            ###################
            ## plot ROC_method.pdf
            fig2 = plt.figure(2)
            tn = 3.2E9-(tp+fn+fp)
            fpr = fp / (fp+tn)
            plt.plot(fpr, sn, label=legend_label, marker='.')
            for i, label in enumerate(point_labels):
                lst_fpr = list(fpr)
                lst_sn = list(sn)
                plt.text(lst_fpr[i], lst_sn[i], label, fontsize=8,color = 'k',fontweight = 'bold' ) 
                plt.xlabel('fpr')
                plt.ylabel('sn')
                plt.title("rna_cov="+str(name))
                plt.ylim(ymax=1)
                plt.legend()

        
        pdf1.savefig(fig1)
        pdf2.savefig(fig2)
        plt.close(fig1)
        plt.close(fig2)

    pdf1.close()
    pdf2.close()

    return




def main():

    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs CTAT_mutations benchmarking\n")

    parser.add_argument('--input_accuracy_table', nargs='+', required = True,
                        help="List of input files containing TP, FP, FN, Sn, PPV, etc.")
    parser.add_argument("--output_dir", required=True, help="output directory")

    args = parser.parse_args()

    #create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)


    input_filenames = args.input_accuracy_table
    
    ## ROC and PRs, each method by truth-set defining depth
    comparative_plots_per_depth(input_filenames, args.output_dir)

    ## ROCs and PRs, each truth-set defining depth by method
    comparative_plots_per_method(input_filenames, args.output_dir)
    
        
    sys.exit(0)



if __name__ == "__main__":
    
    main()


