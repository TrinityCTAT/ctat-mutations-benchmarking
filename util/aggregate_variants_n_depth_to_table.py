#!/usr/bin/env python
# -*- coding: utf-8 -*-

######################
# Import Modules 
import os,sys,argparse
import subprocess
import numpy as np
import pandas as pd
import re
import pickle
import csv
import time
import logging
import gzip

import variant_bed_to_read_depth

import pathos.multiprocessing as pm
import pathos

import io
import dask.dataframe as dd

######################

UTILDIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.sep.join([UTILDIR, "../pylib_common"]))
import ctat_util

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)



csv.field_size_limit(100000000)


def open_file_for_reading(filename):
    if re.search("\.gz$", filename):
        return gzip.open(filename, 'rt', encoding="utf-8",
                         errors='ignore')  # t needed for python3 to look like regular text
    else:
        return open(filename, 'rt', encoding="utf-8", errors='ignore')  # regular text file
    


def configure_vcf(vcf, column_name, chr_idx=0, coord_idx=1, refallele_idx=3, varallele_idx=4):
    '''
    Configure the vcf file to prepare for merging with other data 
       Read vcf in as a pandas data frame 
       1. Adjust the chromosome labels 
       2. change SNPs to 'ref:alt' 
    '''

    logger.info('\t\t Loading {} VCF: {}'.format(column_name, vcf))

    fh = open_file_for_reading(vcf)
    df = pd.read_csv(fh, sep='\t', header=None, comment='#',engine='python')
    
    ## Remove "chr" if in front of chromosome number
    df.loc[:,chr_idx] = df.loc[:,chr_idx].map(str)
    df.loc[:,chr_idx] = df.loc[:,0].map(lambda x: x.lstrip('chr'))
    df=df.replace({chr_idx: {"X":"23","Y":"24","MT":"25","M":"25"}})
    df["CHROMPOS"] = df.iloc[:,chr_idx] + ':' +df.iloc[:,coord_idx].map(str)

    ## set column name to SNP ref:alt
    df[column_name] = df.iloc[:,refallele_idx].map(str) + ':' +df.iloc[:,varallele_idx].map(str)
    
    ## use column names from the VCF format.
    df.rename(columns={chr_idx:'CHROM',coord_idx:'POS', refallele_idx:'REF', varallele_idx:'VAR'}, inplace=True)
    
    
    # print(df.head())
    
    return df



def add_variant_attribute_from_vcf(df, vcf_filename, column_name, tmpdir, **kwargs):

    set_chrompos = set(df['CHROMPOS'].tolist())

    tmpfile_vcf = os.path.join(tmpdir, os.path.basename(vcf_filename) + "{}.tmp".format(time.time()) )
    ofh = open(tmpfile_vcf, 'w')
    
    count_found = 0
    
    with ctat_util.open_file_for_reading(vcf_filename) as fh:
        for line in fh:
            if re.match("#", line):
                ofh.write(line)
                continue
            vals = line.split("\t")
            chr_val = vals[0]
            pos_val = vals[1]

            chr_val = chr_val.replace("chr", "")
            
            chrompos = "{}:{}".format(chr_val, pos_val)
            if chrompos in set_chrompos:
                ofh.write(line)
                count_found += 1
                
    ofh.close()

    logger.info("\t\t - identified {}/{} variants in {}".format(count_found, len(set_chrompos), vcf_filename))
    
    var_attribute_df = configure_vcf(tmpfile_vcf, column_name, **kwargs)

    #os.remove(tmpfile_vcf)


    var_attribute_df = var_attribute_df[['CHROMPOS', column_name]]
    df = pd.merge(df, var_attribute_df, how='left', on='CHROMPOS')
    
    # Check 
    logger.info("\t\t - merged on {} from {}".format(column_name, vcf_filename))
    # print(df.head())
    
    return df

def add_variant_attribute_from_vcf(df, vcf_filename, bed_filename, column_name, cpu, **kwargs):
    #--------------------------------------------------------
    # Get the header for the vcf file and hold it as a string
    #-------------------------------------------------------- 
    header = []
    with ctat_util.open_file_for_reading(vcf_filename) as fh:
        for line in fh:
            if line[0] == '#':
                header.append(line)
            else: break
    header = "".join(header)


    def processing(args):
        '''
        Function to use in the multiprocessing 
        Uses Tabix to split the variants by chromosome,
        Then used Betools intersect to get the variants found in both the vcf of interest and the attributes of interest 
        '''
        chrom, vcf, header  = args
        # command for Tabix
        cmd = "tabix {} chr{}".format(vcf, chrom)
        output1 = subprocess.check_output(cmd, shell=True).decode('utf-8')
        output1 = header + output1
        # command for bedtools intersect 
        cmd = "bedtools intersect -a stdin -b {}".format(bed_filename)
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf8', shell=True).communicate(input = output1)[0]

        return output

    # get the list of chromosomes in the vcf 
    chromosomes = df.CHROM.unique()
    # create list of arguments to pass to multiprocessing
    passing_args = [(i, vcf_filename,header) for i in chromosomes]

    #------------------------
    # Run Multiprocessing 
    #------------------------
    p = pathos.pools._ProcessPool(cpu)
    running = p.map_async(processing, passing_args)
    p.close()
    p.join()

    # process the results from multiprocessing as a string, then convert into a pandas data frame 
    results = "".join(running._value)
    StringData = io.StringIO(results)
    df2 = pd.read_csv(StringData, sep="\t")


    #----------------------------------------------------
    # Processing the output variants from multiprocessing 
    #----------------------------------------------------
    def process_vcf(df, column_name):
        '''
        Function for processing the variant data
            Remove "chr" if in front of chromosome number
        '''
        chr_idx=0
        coord_idx=1
        refallele_idx=3
        varallele_idx=4
        df.iloc[:,chr_idx] = df.iloc[:,chr_idx].map(str)
        df.iloc[:,chr_idx] = df.iloc[:,0].map(lambda x: x.lstrip('chr'))
        df = df.replace({chr_idx: {"X":"23","Y":"24","MT":"25","M":"25"}})
        df["CHROMPOS"] = df.iloc[:,chr_idx] + ':' +df.iloc[:,coord_idx].map(str)
        df[column_name] = df.iloc[:,refallele_idx].map(str) + ':' +df.iloc[:,varallele_idx].map(str)
        df = df[["CHROMPOS", column_name]]
        # df = df.set_index("CHROMPOS")
        return df

    # process the data 
    df2 = process_vcf(df2, column_name)
    # merge the two data frames to include the variant annotation information 
    df = dd.merge(df, df2, how='left', on='CHROMPOS') #left_index=True, right_index=True)
    
    logger.info("\t\t - merged on {} from {}".format(column_name, vcf_filename))
    print(df.head())
    
    return df



def old_add_variant_attribute_from_vcf (df, vcf_filename, column_name, **kwargs):

    """
    old code... replaced by the above method: add_variant_attribute_from_vcf
    retained here just in case...
    """
        
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
    # print(df.head())

    return df


def check_if_chr_bam(filename):
    
    if not os.path.exists(filename):
        errmsg = "Error, cannot find file: {}".format(filename)
        logger.critical(errmsg)
        raise RuntimeError(errmsg)
    
    cmd = 'samtools idxstats '+ filename +'| cut -f1'
    output = (subprocess.check_output(cmd, shell=True))
    # logger.info(output)
    m_chr = re.findall(r'Y\\n(.*)\\n', str(output))[0].split('\\n')[0] ## seems brittle... //TODO: more rigorous way?
    
    chr_prefix_flag = 'chr' in str(output)
    
    logger.info("\t\t chr_prefix_flag: {}, m_chr: {}".format(chr_prefix_flag, m_chr))
    
    return chr_prefix_flag, m_chr



def create_bed_file(df, bam_filename, bed_filename):
    #########################################################################
    # Create a bed file of SNP locations, used to extract coverage depth info.
    ########################################################################
    
    logger.info('\t\tCreating BED file: {}'.format(bed_filename))
    # create a dataframe of chromosome and position 
    df_bed = df[['CHROM', 'POS']].copy()
    df_bed.columns = ['Chr', 'Pos'] # using different ones for the bed

    # Add the position information to the dataframe 
    df_bed['Pos'] = df_bed['Pos'].map(int) # store positions as integers.
    df_bed['Pos-1'] = (df_bed['Pos']-1)

    ## Check if "chr" in front of chromosome name in the bam file
    chr_prefix_flag, m_chr_pred = check_if_chr_bam(bam_filename)
    # convert numbers to X,Y and MT 
    df_bed = df_bed.replace({"Chr": {"23":"X", "24":"Y", "25":m_chr_pred} })
    # remvoe the chr in from of the M so not to get double "chr"
    if chr_prefix_flag:
        df_bed['Chr'] = df_bed['Chr'].map(lambda x: x.lstrip('chr'))
    # add chr 
    if chr_prefix_flag:
        df_bed['Chr'] = 'chr' + df_bed['Chr'].map(str) # add chr at the beginning of chr names

    #-----------------------------
    # remove non-chromosome rows 
    #-----------------------------
    removeID = [i for i in df_bed["Chr"].unique() if len(i)>5]
    uniq = df_bed['Chr'].unique()
    # print("\t\tRemoving the following rows:", removeID)
    keep = set(list(uniq)) - set(removeID)
    df_bed = df_bed[df_bed.Chr.isin(keep)]

    # sort by chr, pos
    df_bed = df_bed.sort_values(['Chr', 'Pos'], ascending=[True, True])
    # Save teh created BED file 
    df_bed.to_csv(bed_filename, sep ='\t', index=False, header=False, na_rep='NA', columns=['Chr', 'Pos-1', 'Pos'])
    
    return df_bed


def generate_depths_file(variants_bed_filename, bam_filename, cpu, depth_filename):
    '''
    Generate the depth values for each on the variants in the VCF using given BAM file 
        Uses the locations in. the bed file to find the depths in the BAM file. 
    '''

    variant_bed_to_read_depth.snp_depths(bed_filename = variants_bed_filename, 
                                         bam_filename = bam_filename, 
                                         cpu          = cpu, 
                                         outputfile   = depth_filename)


        

def add_depth_info(df, df_depth, column_name, bam_filename):
    '''
     Merge depth values with SNP values
    '''
    logger.info("\t\t Adding Depths to SNPs")
    
    # Check if "chr" in front of chromosome name in the bam file
    chr_prefix_flag, m_chr_pred = check_if_chr_bam(bam_filename)
    if chr_prefix_flag:
        df_depth.loc[:,0] = df_depth.loc[:,0].map(lambda x: x.lstrip('chr'))

    df_depth=df_depth.replace({0: {"X":"23", "Y":"24", "MT":"25", "M":"25"} })
    df_depth["CHROMPOS"] = df_depth.iloc[:,0].map(str) + ':' + df_depth.iloc[:,1].map(str)
    df_depth[column_name] = df_depth.iloc[:,2]

    ## simplify to just the 2 columns we want
    df_depth = df_depth[['CHROMPOS', column_name]]

    #df_depth.to_csv('tmp.debug.depth_info.csv', sep = '\t', index=False, na_rep='NA') ## debugging
    
    ## merge into the original data frame
    df = pd.merge(df, df_depth, how='left', on='CHROMPOS')

    return df


def restrict_regions(input_vcf_file, restrict_regions_bed_file, tmpdir):

    restricted_vcf_file = os.path.join(tmpdir, os.path.basename(input_vcf_file) + ".restricted.vcf")
    
    if os.path.exists(restricted_vcf_file):
        logger.info(" \t\t Restricted vcf file {} already exists. using it.".format(restricted_vcf_file))
    else:
        cmd = "bcftools filter -T {} {} -o {}".format(restrict_regions_bed_file, input_vcf_file, restricted_vcf_file)
        subprocess.check_call(cmd, shell=True)

    return restricted_vcf_file



def log_exec_time(time_start, msg, logger):
    logger.info(msg + " {:.2f} minutes".format( (time.time()-time_start)/60.0))



def main():

    ## Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = "Convert vcf file to benchmarking input file\n")
    parser.add_argument('--pred_vcf', required = True, dest = 'pred_vcf', help="input prediction vcf")
    parser.add_argument('--truth_vcf', required = True, help="reference truth set vcf")
    parser.add_argument('--pred_bam', required = True , dest = 'bam_filename', help="Alignment file of reads with the reference")     
    parser.add_argument('--out', required = True , dest = 'out_file', help="Output filename")


    parser.add_argument('--tmpdir', required=False, default="tmpdir.{}".format(time.time()), help="directory for temporary files")
    parser.add_argument('--cpu', required = False , dest = 'cpu', default=1, type=int, help="CPU count for multiprocessing.")

    parser.add_argument('--exome_bam', help="Alignment file of the exome with the reference")

    parser.add_argument("--restrict_regions_bed", help="bed file containing regions to restrict analysis to. (ie. high confidence regions)")
      
    ## annotate certain sites
    parser.add_argument('--dbsnp', dest = 'dbsnp', help="input dbsnp file")
    parser.add_argument('--rna_editing', dest = 'rna_editing', help="input rna-editing file")
    parser.add_argument('--cosmic', dest = 'cosmic', help="input cosmic file")

    args = parser.parse_args()


    pred_vcf = args.pred_vcf
    truth_vcf = args.truth_vcf
    bam_filename = args.bam_filename
    cpu = args.cpu

    tmpdir = args.tmpdir
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    logger.info("Agregating Variants")
    ## apply region restrictions to input vcfs
    logger.info(" \tRestricting:")
    if args.restrict_regions_bed:
        logger.info(" \t \t - Restricting {} positions according to regions {}".format(pred_vcf, args.restrict_regions_bed))
        pred_vcf = restrict_regions(pred_vcf, args.restrict_regions_bed, tmpdir)

        logger.info(" \t \t - Restricting {} positions according to regions {}".format(truth_vcf, args.restrict_regions_bed))
        truth_vcf = restrict_regions(truth_vcf, args.restrict_regions_bed, tmpdir)
     

    #---------------------------------------------------------
    ## Load and configure input vcfs, then merge into single dataframe 
    ##     configure: PREDICTED VCF and REFERENCE VCF
    #---------------------------------------------------------
    logger.info(" \tParsing:")
    logger.info(" \t\t - Parsing rna-vcf: {}".format(pred_vcf))
    tstart = time.time()
    df_pred = configure_vcf(pred_vcf,'RNA_SNP')  ## load prediction vcf
    log_exec_time(tstart, " \t\t configuring RNA snp", logger)
    
    logger.info(" \t\t - Parsing truth_vcf: {}".format(truth_vcf))
    tstart = time.time()
    df_ref = configure_vcf(truth_vcf,'Ref_SNP')
    log_exec_time(tstart, "\t\t configure truth snps", logger)
    
    ## merge dfs
    logger.info("\t Merging:")
    logger.info(" \t\t Merging RNA_SNP and Ref_SNP data frames")
    tstart = time.time()
    df = pd.merge(pd.DataFrame(df_pred, columns = ['CHROMPOS','CHROM', 'POS', 'RNA_SNP']),
                pd.DataFrame(df_ref, columns = ['CHROMPOS', 'CHROM' ,'POS', 'Ref_SNP']),
                on=['CHROMPOS', 'CHROM', 'POS'] ,how='outer') 
    log_exec_time(tstart, "\t\t merge variant calls and truth set", logger)

    # print(df.head())
        
    # df.to_csv("df_test1.vcf",sep = '\t', index=False, na_rep='NA')

    logger.info("\t Merging with given VCF's:")

    #------------------
    # make bed file
    #------------------
    logger.info("\t BED file:")
    variants_bed_filename = os.path.join(tmpdir, os.path.basename(bam_filename) + ".{}_count.variants_pos.bed".format(len(df)))
    if os.path.exists(variants_bed_filename):
        logger.info("\t\t - reusing bed file: {}".format(variants_bed_filename))
    else:
        tstart = time.time()
        create_bed_file(df, bam_filename, variants_bed_filename)
        log_exec_time(tstart, "\t\t rnaseq BED added", logger)




    #--------------------------------------
    # Merge with the DBSNP vcf file to see common variants 
    #--------------------------------------
    time_start = time.time()
    if args.dbsnp:
        logger.info("\t\t Adding dbsnp info")
        tstart = time.time()
        df = add_variant_attribute_from_vcf(df            = df, 
                                            vcf_filename  = args.dbsnp, 
                                            bed_filename  = variants_bed_filename,
                                            column_name   = 'dbsnp_SNP',
                                            cpu = cpu)
        log_exec_time(tstart, "\t\t dbsnp addition", logger)

    #--------------------------------------
    # Merge with the rna-editing vcf file 
    #--------------------------------------
    time_start = time.time()
    if args.rna_editing:
        logger.info(" \t\t Adding rna-editing info")
        tstart = time.time()
        df = add_variant_attribute_from_vcf(df            = df, 
                                            vcf_filename  = args.rna_editing, 
                                            bed_filename  = variants_bed_filename, 
                                            column_name   = 'rnaediting_SNP',
                                            cpu = cpu)
        log_exec_time(tstart, "\t\t rna-editing addition", logger)

    #--------------------------------------
    # Merge with the COSMIC vcf file
    #--------------------------------------
    time_start = time.time()
    if args.cosmic:
        logger.info(" \t\t Adding cosmic info")
        tstart = time.time()
        df = add_variant_attribute_from_vcf(df            = df, 
                                            vcf_filename  = args.cosmic, 
                                            bed_filename  = variants_bed_filename, 
                                            column_name   = 'cosmic_SNP',
                                            cpu = cpu)
        log_exec_time(tstart, "\t\t cosmic addition", logger)


    
    #------------------
    # make DEPTH file
    #------------------
    logger.info("\t Depth file:")
    depth_filename = variants_bed_filename + ".depth"
    if os.path.exists(depth_filename):
        logger.info("\t\t - reusing depth file: {}".format(depth_filename))
    else:
        tstart = time.time()
        generate_depths_file(variants_bed_filename, bam_filename, cpu, depth_filename)
        log_exec_time(tstart, "\t rnaseq depth created", logger)

    df_depth_rna = pd.read_csv(depth_filename, sep='\t', header=None, comment='#',engine='python')
    
    #------------------
    # Merge Depths with df
    #------------------
    logger.info("\t Merging Depth file:")
    tstart = time.time()
    df = add_depth_info(df, df_depth_rna, 'RNAseq_Depth', bam_filename)
    log_exec_time(tstart, "\t\t rnaseq depth added", logger)
    

    #------------------
    # Make BED file
    #------------------
    # Exome 
    if args.exome_bam:
        logger.info(" \t Exome Depth")
        exome_depth_filename = os.path.join(tmpdir, os.path.basename(args.exome_bam) + ".{}_count.variants_pos.exome.depth".format(len(df)))
        if os.path.exists(exome_depth_filename):
            logger.info("\t\t - reusing depth file: {}".format(exome_depth_filename))
        else:
            tstart = time.time()
            generate_depths_file(variants_bed_filename, args.exome_bam, cpu, exome_depth_filename)
            log_exec_time(tstart, "\t\t Exome depth created", logger)

        df_depth_exome = pd.read_csv(exome_depth_filename, sep='\t', header=None, comment='#',engine='python')
        
        #------------------
        # Merge Depths with df
        #------------------
        tstart = time.time()
        df = add_depth_info(df, df_depth_exome,'Exome_Depth', bam_filename)
        log_exec_time(tstart, "\t\t Exome depth added", logger)


     ## make filtering a separate script.
     
     ##filter out the indels
     #if args.filter_indels :
     #     print('Filtering Indels ... ')
     #     mask = (df['RNA_SNP'].astype(str).str.len() <= 3) & (df['Ref_SNP'].astype(str).str.len() <= 3)
     #     df = df.loc[mask]


    logger.info("\t\t outputting final aggregated data table")
    tstart = time.time()
    df.to_csv(args.out_file,sep = '\t', index=False, na_rep='NA')
    log_exec_time(tstart, "\t\t aggregatd table reporting", logger)
    # logger.info('Done!')
     
     
if __name__ == '__main__':
    main()

