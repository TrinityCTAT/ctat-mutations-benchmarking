#!/usr/bin/env python

import argparse
import pandas as pd
import sys, os
import ntpath


def evaluate_accuracy_statistics(input_table_filenames, minimum_coverage, filter_coverage, output_dir):
    ## Generates statistics for plotting ROC/Precision Recall curves
    # Create a file that for prediction statistics for each variant 

    master_data_dict = {}

    # loop over all the files given 
    for input_table_filename in input_table_filenames:
        dirname, filename = ntpath.split(input_table_filename)
        
        # Set the file output name 
        output_filename = os.path.join(output_dir, filename + ".accuracy_stats")
        print("creating the file: ", output_filename)
        
        master_data = {}
        
        # Read in the file and sort the values by DEPTH
        df = pd.read_csv(input_table_filename, sep='\t', low_memory=False)
        df = df[['CHROMPOS','RNA_SNP', 'Ref_SNP', 'RNAseq_Depth', 'Class']]
        df = df.drop_duplicates(keep='last') ## important, dont want to multi-count chrom positions expanded due to data frame expansions w/ attributes. (ie. dbsnp, rnaediting, etc.)
        
        df = df.sort_values('RNAseq_Depth')

        # loop over each of the minimal coverages 
        for ts_min_rna_cov in minimum_coverage:

            ## defines the truth set (ts).
            tmp_df = df[df['RNAseq_Depth'] >= ts_min_rna_cov].copy()
            
            sn_data = []
            sp_data = []
            # total variants in reference 
            total = len(tmp_df[tmp_df['Ref_SNP'].notna()])
            # loop over each of the filter coverages 
            for eval_min_rna_cov in filter_coverage:

                if eval_min_rna_cov < ts_min_rna_cov:
                    continue
                
                tiny_df = tmp_df[tmp_df['RNAseq_Depth'] >= eval_min_rna_cov].copy()
                data = tiny_df.groupby('Class').Class.count().to_dict()
                tp = data['TP'] if 'TP' in data else 0
                fp = data['FP'] if 'FP' in data else 0
                fn = total - tp
                try:
                    sn = float(tp)/float(tp + fn)
                    ppv = (float(tp)/float(tp+fp))
                    key = str(ts_min_rna_cov) + '-' + str(eval_min_rna_cov)
                    master_data[key] = [ts_min_rna_cov, eval_min_rna_cov, tp, fp, fn, sn, ppv]
                except ZeroDivisionError:
                    print("-div by zero... skipping this calc.", file=sys.stderr)
                
        master_data_dict[output_filename] = master_data

    
    #-----------------------
    ## Write master files
    #-----------------------
    for f_name in master_data_dict:
        f_out = open(f_name,'w')
        f_out.write('ts_min_rna_cov\teval_min_rna_cov\ttp\tfp\tfn\tsn\tppv\n')
        col = ['ts_min_rna_cov','eval_min_rna_cov','tp','fp','fn','sn','ppv']
        master_data = master_data_dict[f_name]
        df = pd.DataFrame(list(master_data.values()))
        df.columns = col
        df = df.sort_values(by=['ts_min_rna_cov', 'eval_min_rna_cov'], ascending=[True, True])
        df.to_csv(f_name,sep='\t',index=False,na_rep='NA')

    return



def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Performs CTAT_mutations benchmarking, Creates the Summary statistics. \n")

    parser.add_argument('--input_tables', nargs='+', required = True, help="List of input files")

    parser.add_argument('--output_dir', required = True, help="output directory")

    parser.add_argument('--truth_min_rna_cov_thresholds',
                        nargs='+',
                        type=int,
                        help="List of minimum rna-seq coverage values to define the truth set (ie. ref snps must be expressed to be included)",
                        default=[1, 3, 5, 10, 20])

    parser.add_argument('--eval_min_rna_cov_thresholds',
                        nargs='+',
                        type=int,
                        help="List of filteration coverage values",
                        default=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60])

    
    args = parser.parse_args()
    
        
    #create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    evaluate_accuracy_statistics(args.input_tables, args.truth_min_rna_cov_thresholds, args.eval_min_rna_cov_thresholds, args.output_dir)

    sys.exit(0)
    


if __name__ == "__main__":

    main()
