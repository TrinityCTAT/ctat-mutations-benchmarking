#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Vrushali Fangal"
__copyright__ = "Copyright 2014"
__credits__ = [ "Maxwell Brown", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Vrushali Fangal"
__email__ = "vrushali@broadinstitute.org"
__status__ = "Development"


## Import Libraries

import os, sys, csv
import glob  # File
import warnings
import re

warnings.filterwarnings("ignore")

## Import python libraries
import matplotlib.pyplot as plt
import pandas as pd  
import numpy as np   
import argparse    
import seaborn as sns
sns.set(style="whitegrid")


## Check if python 3 is imported
if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)

def filename_to_type_string(longname):
    
    longname = os.path.basename(longname)
    substr_size = 10

    # breast.RF-regressor.vcf.gz.aggregated.dat.filtered.scored.accuracy_stats
    longname = longname.replace(".aggregated.dat.filtered.scored.accuracy_stats", "")
    longname = longname.replace(".vcf", "")
    longname = re.sub(".gz$", "", longname)

    newname_pts = list()
    for pts in longname.split("."):
        
        if len(pts) > substr_size:
            pts = pts[0:substr_size]
            
        newname_pts.append(pts)
    
    longname = ".".join(newname_pts)

    return longname




def plot_f1_score(input_table_filenames, output_dir):
    
    df_accuracy_stats = pd.DataFrame()
    seen = set()
    for filename in input_table_filenames:
        df = pd.read_csv(filename, sep='\t')
        df = df[df['ts_min_rna_cov']==10]
        typename = filename_to_type_string(filename)
        if typename in seen:
            raise RuntimeError("Error, converted {} to {} for unique type but not unique".format(filename, typename))
        seen.add(typename)
        df['Type'] = typename
        df_accuracy_stats = pd.concat([df_accuracy_stats,df],ignore_index=True)
    
    ## Compute F1 score (harmonic mean)
    df_accuracy_stats['F1'] = (2 * df_accuracy_stats['sn'] * df_accuracy_stats['ppv'])/ (df_accuracy_stats['sn'] + df_accuracy_stats['ppv'])
    
    ## Compute accuracy
    df_accuracy_stats['Accuracy'] = ( df_accuracy_stats['tp'])/ ( df_accuracy_stats['tp'] + df_accuracy_stats['fp'] + df_accuracy_stats['fn'] )
    
    ## Mathews correlation coefficient : (TP*TP - FP*FN)/ (TP+FP)*(TP+FN) 
    df_accuracy_stats['MCC'] = ( (df_accuracy_stats['tp'] * df_accuracy_stats['tp']) - (df_accuracy_stats['fp'] * df_accuracy_stats['fn'])) / \
                                 ( (df_accuracy_stats['tp'] + df_accuracy_stats['fp']) * (df_accuracy_stats['tp'] + df_accuracy_stats['fn'] ) )


    ######################################
    ################ F1 #################
    ######################################
    
    ## Rank boxplots
    plt.figure()
    medianprops = {'color': 'magenta', 'linewidth': 2}
    boxprops = {'edgecolor': 'k', 'linestyle': '-', 'color':'lightskyblue'}
    flierprops = dict(color='white',marker='o', markersize=1,markerfacecolor='auto', markeredgecolor='k')
    whiskerprops=dict(color='k',linewidth=1,linestyle='dashed')
    capprops=dict(color='black',linestyle= '-',linewidth=0)
    
    
    ranks = df_accuracy_stats.groupby("Type")['F1'].max().sort_values()[::-1].index
    f1 = df_accuracy_stats[df_accuracy_stats['eval_min_rna_cov']==10].sort_values(by=['F1'])[::-1]
    f1['F1'] = f1['F1']*100
    sns_plot = sns.barplot(y="F1", x="Type", data= f1, color='hotpink',edgecolor=(0.5,0.5,0.5),linewidth=2)
    
    #sns_plot = sns.boxplot(y="F1", x="Type", order=ranks, data=df_accuracy_stats, palette='tab10')
    #sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)

    print(f1[['Type','F1']])
    plt.ylim([np.min(f1['F1'])-3,np.max(f1['F1'])+2])


    for p in sns_plot.patches:
        sns_plot.annotate(format(p.get_height(), '.2f'),
                          (p.get_x() + p.get_width() / 2., p.get_height()),
                          ha = 'center', #va = 'center',
                          xytext = (0, 23), rotation = 90, fontsize=10,
                          textcoords = 'offset points') 

    

    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig = sns_plot.get_figure()

    plt.xlabel('Algorithm')
    plt.title("F1 plot at Coverage 10")
    fig.savefig(os.path.join(output_dir,"F1.pdf"),bbox_inches = "tight")


    ######################################
    ########### Accuracy #################
    ######################################

    ## Rank accuracy 
    plt.figure()
    ranks = df_accuracy_stats.groupby("Type")['Accuracy'].max().sort_values()[::-1].index
    acc = df_accuracy_stats[df_accuracy_stats['eval_min_rna_cov']==10].sort_values(by=['Accuracy'])[::-1]
    acc.to_csv(os.path.join(output_dir, 'accuracy_stats.tsv'), sep='\t')
    #sns_plot = sns.pointplot(y="Accuracy", x="Type", data= acc, join=True, sort=False, markers=["o"], color='k')
    #sns_plot = sns.pointplot(linestyles=["-"],y="Accuracy", x="Type", data= acc, join=True, sort=False, markers=["o"], palette='Set1')
    sns_plot = sns.barplot(y="Accuracy", x="Type", data= acc, palette='Set1')
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig = sns_plot.get_figure()
    plt.title("Accuracy plot at Coverage 10")
    fig.savefig(os.path.join(output_dir,"Accuracy.pdf"),bbox_inches = "tight")


    ######################################
    ################ MCC #################
    ######################################

    ## Rank MCC
    plt.figure()
    ranks = df_accuracy_stats.groupby("Type")['MCC'].max().sort_values()[::-1].index
    sns_plot = sns.boxplot(y="MCC", x="Type", order=ranks, data=df_accuracy_stats, palette='tab10')
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig = sns_plot.get_figure()
    plt.title("MCC plot at Coverage 10")
    fig.savefig(os.path.join(output_dir,"MCC.pdf"), bbox_inches = "tight")


def main():

    ## Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Plots F1 score \n")
    parser.add_argument('--input_accuracy_table', nargs='+', required = True, help="List of input files")
    parser.add_argument('--output_dir', required = True, help="output directory")
    args = parser.parse_args()
    plot_f1_score(args.input_accuracy_table, args.output_dir)

    sys.exit(0)

if __name__ == "__main__":

    main()
