#!/usr/bin/env python

import numpy as np
import argparse
import pandas as pd
import sys, os


def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Filter variants by attributes. \n")

    parser.add_argument('--input_table', required = True, help="input table of aggregated variant info")
    parser.add_argument('--output_table', required = True, help="output table of aggregated variant info")
        
    parser.add_argument('--min_RNAseq_Depth',  type = int , help = 'Depth value for SNP analysis \n', default=1)

    parser.add_argument("--remove_indels", action='store_true', default=False, help='remove indels')

    parser.add_argument('--min_exome_depth',  type = int , help = 'min exome coverage depth (set for exome analysis, requires Exome_Depth column) \n', default=0)
    
    # options to remove specific variants 
    parser.add_argument('--remove_intersect',
                        nargs='+',
                        required = False,
                        help= "Filter out variants that intersect with rnaediting and/or dbsnp",
                        default = '', choices = ['rnaediting', 'dbsnp'])
    
    args = parser.parse_args()

    
    df = pd.read_csv(args.input_table, sep='\t', low_memory=False)

    df = df[ df['RNAseq_Depth'] >= args.min_RNAseq_Depth ]

    if args.min_exome_depth > 0:
        df = df[ df['Exome_Depth'] >= args.min_exome_depth ]
    
    ## remove indels
    if args.remove_indels:
        df = df[df['RNA_SNP'].map(str).apply(len).lt(4)] #remove indels in rna snp calls
        df = df[df['Ref_SNP'].map(str).apply(len).lt(4)] #remove indels in reference
        #df = df[df['RNA_SNP'].str.len().lt(4)] #remove indels in rna snp calls
        #df = df[df['Ref_SNP'].str.len().lt(4)] #remove indels in reference
    else:
        print("WARNING - indels are not being removed.  Use --remove_indels to engage", file=sys.stderr)


    if 'rnaediting' in args.remove_intersect:
        # retain those that are not found in rnaediting
        df = df[ df['rnaediting_SNP'].isna() ]
            
    if 'dbsnp' in args.remove_intersect:
        # retain those that are not found in dbsnp
        df = df[ df['dbsnp_SNP'].isna() ]


    f_out = open(args.output_table, 'w')
    df.to_csv(f_out,sep='\t',index=False,na_rep='NA')

    sys.exit(0)




if __name__ == "__main__":

    main()
