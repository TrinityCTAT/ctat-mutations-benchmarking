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

warnings.filterwarnings("ignore")

## Import python libraries
import matplotlib.pyplot as plt
import pandas as pd  
import numpy as np   
import argparse    
import seaborn as sns
sns.set(style="whitegrid")
import ntpath


## Check if python 3 is imported
if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)


def plot_f1_score(input_table_filenames, output_dir):
    
    df_accuracy_stats = pd.DataFrame()
    for file in input_table_filenames:
        df = pd.read_csv(file, sep='\t')
        df = df[df['rna_cov']==10]
        df['Type'] = '_'.join(ntpath.basename(file).split('_',2)[:2])
        df_accuracy_stats = pd.concat([df_accuracy_stats,df],ignore_index=True)
    
    ## Compute F1 score (harmonic mean)
    df_accuracy_stats['F1'] = (2 * df_accuracy_stats['sn'] * df_accuracy_stats['ppv'])/ (df_accuracy_stats['sn'] + df_accuracy_stats['ppv'])
    
    ## Compute accuracy
    df_accuracy_stats['Accuracy'] = ( df_accuracy_stats['tp'])/ ( df_accuracy_stats['tp'] + df_accuracy_stats['fp'] + df_accuracy_stats['fn'] )
    
    ## Rank boxplots
    plt.figure()
    ranks = df_accuracy_stats.groupby("Type")['F1'].max().sort_values()[::-1].index
    sns_plot = sns.boxplot(y="F1", x="Type", order=ranks, data=df_accuracy_stats, palette='tab10')
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig = sns_plot.get_figure()
    plt.title("F1 plot at Coverage 10")
    fig.savefig(os.path.join(output_dir,"F1.pdf"),bbox_inches = "tight")

    ## Rank accuracy 
    plt.figure()
    ranks = df_accuracy_stats.groupby("Type")['Accuracy'].max().sort_values()[::-1].index
    acc = df_accuracy_stats[df_accuracy_stats['cov']==10].sort_values(by=['Accuracy'])[::-1]
    sns_plot = sns.pointplot(y="Accuracy", x="Type", data= acc, join=True, sort=False, markers=["o"], color='k')
    sns_plot = sns.pointplot(linestyles=["-"],y="Accuracy", x="Type", data= acc, join=True, sort=False, markers=["o"], palette='Set1')
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig = sns_plot.get_figure()
    plt.title("Accuracy plot at Coverage 10")
    fig.savefig(os.path.join(output_dir,"Accuracy.pdf"),bbox_inches = "tight")


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
