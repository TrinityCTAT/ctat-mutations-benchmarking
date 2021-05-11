#!/usr/bin/env python


import os,sys,argparse
import subprocess
import numpy as np
import pandas as pd
import re
import pickle
import csv

csv.field_size_limit(100000000)

def configure_vcf(vcf, column_name, chr_idx=0, coord_idx=1, refallele_idx=3, varallele_idx=4):
     '''
     Configure the vcf file to prepare for merging with other data 
     '''

     print('Loading {} VCF: {}'.format(column_name, vcf))
     
     df = pd.read_csv(vcf, sep='\t', header=None, comment='#',engine='python')
     
     ## Remove "chr" if in front of chromosome number
     df.loc[:,chr_idx] = df.loc[:,chr_idx].map(str)
     df.loc[:,chr_idx] = df.loc[:,0].map(lambda x: x.lstrip('chr'))
     df=df.replace({chr_idx: {"X":"23","Y":"24","MT":"25","M":"25"}})
     df["CHROMPOS"] = df.iloc[:,chr_idx] + ':' +df.iloc[:,coord_idx].map(str)

     ## set column name to SNP ref:alt
     df[column_name] = df.iloc[:,refallele_idx].map(str) + ':' +df.iloc[:,varallele_idx].map(str)
     
     ## use column names from the VCF format.
     df.rename(columns={chr_idx:'CHROM',coord_idx:'POS', refallele_idx:'REF', varallele_idx:'VAR'}, inplace=True)
     
     
     print(df.head())
     
     return df



def add_variant_attribute_from_vcf(df, vcf_filename, column_name, **kwargs):

    # check for pickle file and load 
    pickle_path = os.path.basename(vcf_filename) + ".pickle"
    print(pickle_path)
    
    var_attribute_df = None
    if (os.path.exists(pickle_path)):
        print('Loading Pickle: {}'.format(pickle_path))
        inPickle = open(pickle_path,'rb')
        var_attribute_df = pickle.load(inPickle)
        inPickle.close()
        
    # if not found, read and create the pickle 
    else:
        print('Creating Pickle: {}'.format(pickle_path))
        var_attribute_df = configure_vcf(vcf_filename, column_name, **kwargs)
        outPickle = open(pickle_path, 'wb')
        pickle.dump(var_attribute_df, outPickle)
        outPickle.close()
    
    var_attribute_df = var_attribute_df[['CHROMPOS', column_name]]
    df = pd.merge(df, var_attribute_df, how='left', on='CHROMPOS')
    
    # Check 
    print(df.head())

    return df


def check_if_chr_bam(filename):
     cmd = 'samtools idxstats '+ filename +'| cut -f1'
     output = (subprocess.check_output(cmd, shell=True))
     print(output)
     m_chr = re.findall(r'Y\\n(.*)\\n', str(output))[0].split('\\n')[0] ## seems brittle... //TODO: more rigorous way?
     
     chr_prefix_flag = 'chr' in str(output)

     print('****',chr_prefix_flag, m_chr)
     
     return chr_prefix_flag, m_chr

def add_depth_info(df, bam_filename, column_name):

    print(df.head())
    
    #########################################################################
    # Create a bed file of SNP locations, used to extract coverage depth info.
    ########################################################################
    
    df_bed = df[['CHROM', 'POS']].copy()
    df_bed.columns = ['Chr', 'Pos'] # using different ones for the bed

    print(df_bed.head())

    df_bed['Pos'] = df_bed['Pos'].map(int) # store positions as integers.
    df_bed['Pos-1'] = (df_bed['Pos']-1)
    
    ## Check if "chr" in front of chromosome name in the bam file
    chr_prefix_flag, m_chr_pred = check_if_chr_bam(bam_filename)
        
    if chr_prefix_flag:
        df_bed['Chr'] = 'chr' + df_bed['Chr'].map(str) # add chr at the beginning of chr names
        
    # back-convert X,Y,M from numeric vals used earlier 
    df_bed['Chr'] = df_bed.loc[:,'Chr'].replace('23','X') #replace chr23 with chrX
    df_bed['Chr'] = df_bed.loc[:,'Chr'].replace('24','Y') #replace chr24 with chrY
    df_bed['Chr'] = df_bed.loc[:,'Chr'].replace('25', m_chr_pred) #replace chr25 with chrM


    # sort by chr, pos
    df_bed = df_bed.sort_values(['Chr', 'Pos'], ascending=[True, True])
    
    # make bed file containing target sites for depth calc
    bed_filename = os.path.basename(bam_filename) + ".{}_count.variants_pos.bed".format(len(df_bed))

    if os.path.exists(bed_filename):
        print("-reusing bed file: {}".format(bed_filename))
    else:
        df_bed.to_csv(bed_filename, sep ='\t', index=False, header=False, na_rep='NA', columns=['Chr', 'Pos-1', 'Pos'])
        print('Bed File created: {}'.format(bed_filename))
        
    
    # make depth file
    depth_filename = bed_filename + ".depth"
    if os.path.exists(depth_filename):
        print("-reusing depth file: {}".format(depth_filename))
    else:
        cmd = "samtools depth -b {} {} > {}".format(bed_filename, bam_filename, depth_filename)
        subprocess.check_call(cmd, shell=True)
        print("Depth file created: {}".format(depth_filename))

    
    # Merge depth values with SNP values
    df_depth = pd.read_csv(depth_filename, sep='\t', header=None, low_memory=False)

    if chr_prefix_flag:
        df_depth.loc[:,0] = df_depth.loc[:,0].map(lambda x: x.lstrip('chr'))

    df_depth=df_depth.replace({0: {"X":"23", "Y":"24", "MT":"25", "M":"25"} })
    df_depth["CHROMPOS"] = df_depth.iloc[:,0].map(str) + ':' +df_depth.iloc[:,1].map(str)
    df_depth[column_name] = df_depth.iloc[:,2]

    ## simplify to just the 2 columns we want
    df_depth = df_depth[['CHROMPOS',column_name]]

    ## merge into the original data frame
    df = pd.merge(df, df_depth, how='left', on='CHROMPOS')

    return df


def restrict_regions(input_vcf_file, restrict_regions_bed_file):

    restricted_vcf_file = input_vcf_file + ".restricted.vcf"
    
    cmd = "bcftools filter -T {} {} -o {}".format(restrict_regions_bed_file, input_vcf_file, restricted_vcf_file)

    subprocess.check_call(cmd, shell=True)

    return restricted_vcf_file



def main():

     ## Input arguments
     parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = "Convert vcf file to benchmarking input file\n")
     parser.add_argument('--pred_vcf', required = True, dest = 'pred_vcf', help="input prediction vcf")
     parser.add_argument('--truth_vcf', required = True, help="reference truth set vcf")
     parser.add_argument('--pred_bam', required = True , dest = 'bam_file', help="Alignment file of reads with the reference")     
     parser.add_argument('--out', required = True , dest = 'out_file', help="Output filename")

     parser.add_argument('--exome_bam', help="Alignment file of the exome with the reference")
     
     parser.add_argument("--restrict_regions_bed", help="bed file containing regions to restrict analysis to. (ie. high confidence regions)")
          
     ## filtering options.
     parser.add_argument('--dbsnp', dest = 'dbsnp', help="input dbsnp file")
     parser.add_argument('--rnaediting', dest = 'rnaediting', help="input rnaediting file")
     parser.add_argument('--cosmic', dest = 'cosmic', help="input cosmic file")

     args = parser.parse_args()


     pred_vcf = args.pred_vcf
     truth_vcf = args.truth_vcf

     ## apply region restrictions to input vcfs
     if args.restrict_regions_bed:
         pred_vcf = restrict_regions(pred_vcf, args.restrict_regions_bed)
         truth_vcf = restrict_regions(truth_vcf, args.restrict_regions_bed)
     

     ## Load input vcfs to python df
     
     df_pred = configure_vcf(pred_vcf,'RNA_SNP')  ## load prediction vcf

     df_ref = configure_vcf(truth_vcf,'Ref_SNP')
     
     ## merge dfs
     df = pd.merge(pd.DataFrame(df_pred, columns = ['CHROMPOS','CHROM', 'POS', 'RNA_SNP']),
                   pd.DataFrame(df_ref, columns = ['CHROMPOS', 'CHROM' ,'POS', 'Ref_SNP']),
                   on=['CHROMPOS', 'CHROM', 'POS'] ,how='outer') 
     
          
     #--------------------------------------
     # Merge with the DBSNP vcf file to see common variants 
     #--------------------------------------
     if args.dbsnp:
         df = add_variant_attribute_from_vcf(df, args.dbsnp, 'dbsnp_SNP')
     
     #--------------------------------------
     # Merge with the rnaediting vcf file 
     #--------------------------------------
     if args.rnaediting:
         # rnaediting column idxs are slightly off from expected vcf formatting, requires slight adjustment.
         df = add_variant_attribute_from_vcf(df, args.rnaediting, 'rnaediting_SNP', refallele_idx=2, varallele_idx=3)

     #--------------------------------------
     # Merge with the COSMIC vcf file
     #--------------------------------------
     if args.cosmic:
         df = add_variant_attribute_from_vcf(df, args.cosmic, 'cosmic_SNP')
     
     print('VCFs loaded.')

     ## get rna-seq depth info
     df = add_depth_info(df, args.bam_file, 'RNAseq_Depth')

     if args.exome_bam:
         df = add_depth_info(df, args.exome_bam, "Exome_Depth")


     ## make filtering a separate script.
     
         
     ##filter out the indels
     #if args.filter_indels :
     #     print('Filtering Indels ... ')
     #     mask = (df['RNA_SNP'].astype(str).str.len() <= 3) & (df['Ref_SNP'].astype(str).str.len() <= 3)
     #     df = df.loc[mask]


     df.to_csv(args.out_file,sep = '\t', index=False, na_rep='NA')
     
     print('Done!')
     
     
if __name__ == '__main__':
    main()

