#!/usr/bin/env python

import numpy as np
import argparse
import pandas as pd
import sys, os
import logging


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)



def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Filter variants by attributes. \n")

    parser.add_argument('--input_table', required = True, help="input table of aggregated variant info")
    parser.add_argument('--output_table', required = True, help="output table of aggregated variant info")
        
    parser.add_argument('--min_RNAseq_Depth',  type = int , help = 'Depth value for SNP analysis \n', default=1)

    parser.add_argument("--snvs_only", action='store_true', default=False, help='restrict to snvs')
    parser.add_argument("--indels_only", action='store_true', default=False, help='restrict to indels')

    parser.add_argument('--min_exome_depth',  type = int , help = 'min exome coverage depth (set for exome analysis, requires Exome_Depth column) \n', default=0)
    
    # options to remove specific variants 
    parser.add_argument('--remove_rna_editing',
                        action = 'store_true',
                        help= "Filter out variants that intersect with rna-editing sites",
                        default = False)
    
    args = parser.parse_args()

    
    df = pd.read_csv(args.input_table, sep='\t', low_memory=False)

    df = df[ df['RNAseq_Depth'] >= args.min_RNAseq_Depth ]

    if args.min_exome_depth > 0:
        df = df[ df['Exome_Depth'] >= args.min_exome_depth ]

    if args.snvs_only and args.indels_only:
        raise RuntimeError("Error, cannot enable both --snvs_only and --indels_only... choose one or the other")
    
    ## remove indels
    if args.snvs_only:
        logger.info("-analyzing snvs only")
        df = df[df['RNA_SNP'].map(str).apply(len).lt(4)] #remove indels in rna snp calls
        df = df[df['Ref_SNP'].map(str).apply(len).lt(4)] #remove indels in reference
     
    elif args.indels_only:
        logger.info("-analyzing indels only")
        df = df[df['RNA_SNP'].map(str).apply(len).ne(3) | df['Ref_SNP'].map(str).apply(len).ne(3)]

    if args.remove_rna_editing:
        # retain those that are not found in rnaediting
        df = df[ df['rnaediting_SNP'].isna() ]
            

    f_out = open(args.output_table, 'w')
    df.to_csv(f_out,sep='\t',index=False,na_rep='NA')

    sys.exit(0)




if __name__ == "__main__":

    main()
