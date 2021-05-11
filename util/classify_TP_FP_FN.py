#!/usr/bin/env python


import argparse
import pandas as pd
import sys, os
import numpy as np


def classify_snps(df):
     ## Classifies SNPs as TP, FP and FN based on depth coverage
     df['Class'] = np.NaN

     # in reference but not rna-seq: false negatives
     df.loc[df['RNA_SNP'].isna() & df['Ref_SNP'].notna(),'Class'] = 'FN' #FN

     # in rnaseq, not in refernece: false positives
     df.loc[df['RNA_SNP'].notna() & df['Ref_SNP'].isna(),'Class'] = 'FP' #FP

     # in both, but dont match: false positive (more stringent test... might reconsider this or make optional.)
     df.loc[df['RNA_SNP'].notna() & df['Ref_SNP'].notna() & (df['RNA_SNP'] != df['Ref_SNP']),'Class'] = 'FP' #FP

     # in both and match: true positives
     df.loc[df['RNA_SNP'] == df['Ref_SNP'],'Class'] = 'TP' #TP

     return df




def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs CTAT_mutations benchmarking, Creates the Summary statistics. \n")

    parser.add_argument('--input_table', required=True, help="input file containing aggregated data, requires columns: RNA_SNP and Ref_SNP")
    parser.add_argument('--output_table', required=True, help="output file")

    args = parser.parse_args()


    df = pd.read_csv(args.input_table, sep='\t', low_memory=False)

    df = classify_snps(df)
    
    f_out = open(args.output_table, 'w')
    df.to_csv(f_out,sep='\t',index=False,na_rep='NA')

    sys.exit(0)

if __name__ == "__main__":

    main()
