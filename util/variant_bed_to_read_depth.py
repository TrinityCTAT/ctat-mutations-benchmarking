#!/usr/bin/env python
# -*- coding: utf-8 -*-

######################
# Import Modules 
import os,sys,argparse
import subprocess
import numpy as np
import pandas as pd
import re
import csv
import time
import logging
import pathos.multiprocessing as pm
import pathos
######################


sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../pylib_common"]))

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)



csv.field_size_limit(100000000)

#--------------------------------------------------------------------------------------------------------------------
# Functions associated with getitng the Depth values 
#---------------------------------------------------
def progress_bar(progress_percent):
        #-------------------------
        # Create the Progress bar 
        #-------------------------
        # Create progress bar to monitor the progress of the multiprocessing 
        ## Remove the line from before 
        sys.stdout.write("\r")
        sys.stdout.flush()
        ## print the progress bar and percent 
        if progress_percent == 100:
            sys.stdout.write("[{}{}]{}".format("*" * progress_percent, " "* (100-progress_percent), str(progress_percent)+"%\n"))
        else:
            sys.stdout.write("[{}{}]{}".format("*" * progress_percent, " "* (100-progress_percent), str(progress_percent)+"%"))
        sys.stdout.flush()



def samtools_depths(arg):
    '''
    Find the depths for each location in a bedfile 
        input
            arg: contains a bed file, and a path to a bam file.
        output
            samtools_output: Returns the depths for each location in the bed file 
    '''
    # unpack the inputted arguments 
    subset_df_bed, bam_filename = arg

    bam_index_filename = f"{bam_filename}.bai"
    if not os.path.exists(bam_index_filename):
        cmd = "samtools index {}".format(bam_index_filename)
        subprocess.check_call(cmd, shell=True)
        

    # gather the current Chromosome along with the regions start and end point that the bed file spans 
    chrom = subset_df_bed['Chr'].iloc[0]
    start = subset_df_bed['Pos-1'].iloc[0]
    end = subset_df_bed['Pos'].iloc[-1]

    # join the columns of the bed subset with tabs, then combine the rows by '\n' to make one string
    joined_df_bed = subset_df_bed['Chr'] + "\t" + subset_df_bed['Pos-1'].map(str) + "\t" + subset_df_bed['Pos'].map(str)
    sub_bed = "\n".join(joined_df_bed)

    # set up the region annotation CHROMOSOME:START-STOP
    region = "{}:{}-{}".format(chrom, start, end)

    # Samtools depth command to run in shell using subprocess 
    cmd = "samtools depth -b /dev/stdin -r {} {}".format(region, bam_filename)
    # run the command and collect and return the output 
    samtools_output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf8', shell=True).communicate(input = sub_bed)[0]

    return samtools_output


def generate_depths(df_bed, bam_filename, cpu, depth_filename):
    '''
    Generate the RNA-seq depths using Samtool's samtool depths for each variant 
    use pathos module for multiprocessing on each of the chromosomes
        input
            df_bed:         BED file as pandas dataframe 
            bam_filename:   BAM file used to create the VCF
        output 
            depths for each of the locations in the bed file 
    '''
    logger.info("\t\t Creating Depth file: {}".format(depth_filename))

    bam_index_filename = f"{bam_filename}.bai"
    if not os.path.exists(bam_index_filename):
        cmd = "samtools index {}".format(bam_index_filename)
        subprocess.check_call(cmd, shell=True)

    
    # List of all the chromosomes in the BEDfile
    chrs = list(df_bed['Chr'].unique())
    # count for the number of lines in bed file 
    bed_line_count = len(df_bed)

    def split_chrs(df_bed, i, cpu):
        x = np.array_split(df_bed.loc[df_bed['Chr'] == i], cpu)
        nonempty = [df for df in x if not df.empty]
        return nonempty
    
    # calling_args = [ np.array_split(df_bed.loc[df_bed['Chr'] == i], cpu)  for i in chrs]
    calling_args = [ split_chrs(df_bed, i, cpu)  for i in chrs]
    # make the list flat (make list of lists into a single list)
    passing_args = [(item, bam_filename) for sublist in calling_args for item in sublist]

    #------------------------------------
    # Multiprocessing for each Chromosome
    #------------------------------------
    # get a start time for the multiprocessing for checking purposes 
    start_time = time.time()
    # set up the multiprocessing 
    p = pathos.pools._ProcessPool(cpu)
    running = p.map_async(samtools_depths, passing_args, chunksize = 1)

    # make a progress bar while multiprocessing is running 
    while not running._value.count(None) == 0: 
        # count for job number 
        counts = len(running._value) - running._value.count(None)
        # print(running._value)
        progress_percent = int((counts)/len(running._value) * 100)
        # Print the progress percentage to the terminal
        progress_bar(progress_percent)
        time.sleep(3)

    # close and join the pool 
    p.close()
    p.join()
    # print 100% when done 
    progress_bar(100)

    print("*** %s minutes for depths ***" % ((time.time() - start_time)/60))
    
    # -------------------
    # Make the depth file
    # -------------------

    # fetch the results once finished, and make into a single string for printing 
    final_depth_results = "".join(running._value)
    # write to a file
    new_file = open(depth_filename, "w")
    new_file.write(final_depth_results)
    new_file.close()



def log_exec_time(time_start, msg, logger):
    logger.info(msg + " {:.2f} minutes".format( (time.time()-time_start)/60.0))




def snp_depths(bed_filename, bam_filename, cpu, outputfile):
    '''
    Main function in the depth generating process. 
    Generate the depth values for each on the variants in the VCF using given BAM file 
        Uses the locations in. the bed file to find the depths in the BAM file. 
    '''
            
    logger.info("\t\t - parsing bed file: {}".format(bed_filename))
    df_bed = pd.read_csv(bed_filename, sep='\t', header=None, comment='#',engine='python')
    df_bed.columns = ['Chr', 'Pos-1', 'Pos']
    
    #------------------
    # make DEPTH file
    #------------------
    depth_filename = outputfile
    tstart = time.time()
    df_depth_rna = generate_depths(df_bed, bam_filename, cpu, depth_filename)
    log_exec_time(tstart, "\t\t rnaseq depth created", logger)



def main():
# Input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = "Get samtools depth for variant bed\n")
    parser.add_argument('--var_bed', required = True, help="variant bed file")
    parser.add_argument('--bam', required = True , dest = 'bam_filename', help="Alignment file of reads with the reference")     
    parser.add_argument('--out', required = True , dest = 'out_file', help="Output filename")
    parser.add_argument('--cpu', required = False , dest = 'cpu', default=1, type=int, help="CPU count for multiprocessing.")

    args = parser.parse_args()


    bed_filename = args.var_bed
    bam_filename = args.bam_filename
    cpu = args.cpu
    outputfile = args.out_file

    snp_depths(bed_filename, bam_filename, cpu, outputfile)

if __name__ == '__main__':
    main()
