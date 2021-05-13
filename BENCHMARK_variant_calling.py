#!/usr/bin/env python
# encoding: utf-8



############################
# import modules
############################
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, re, sys, glob
import argparse
import subprocess
import time
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "pylib_common"]))
from Pipeliner import Pipeliner, Command

BASEDIR = os.path.dirname(__file__)
UTILDIR = os.sep.join([BASEDIR, "util"])


def create_output_dir_name(args):
    
    output_dir = str(".".join(["outdir",
                               *["P-" + os.path.basename(pred_vcf) for pred_vcf in args.pred_vcf],
                               "R-" + os.path.basename(args.truth_vcf),
                               "minRcov_{}".format(args.min_RNAseq_Depth)]))
    if args.exome_bam:
        output_dir += ".ExomeMode"

    if args.indels_only:
        output_dir += ".indels_only"
    elif args.snvs_only:
        output_dir += ".snvs_only"

    if args.rna_editing:
        output_dir += ".filt_RNAediting"

            
    return output_dir

        
def audit_command_info(output_dir, argv):

    cmdstr = " ".join(argv)

    cmd_audit_filename = os.path.join(output_dir, "__benchmark_cmd.{}.txt".format(time.time()))

    with open(cmd_audit_filename, 'w') as ofh:
        ofh.write(cmdstr)

    return



def main():
    
    ##############################
    # Command line arguments 
    ##############################

    #add options to inputs
    arg_parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                         description = "Performs CTAT_mutations benchmarking, Creates the Summary statistics. \n")
    
    
    ## required options
    required_opts_group = arg_parser.add_argument_group('required arguments')
    required_opts_group.add_argument('--pred_vcf',
                                     nargs='+',
                                     required = True,
                                     dest = 'pred_vcf',
                                     help="input prediction vcf")
    required_opts_group.add_argument('--truth_vcf', required = True, help="reference truth set vcf")
    required_opts_group.add_argument('--pred_bam', required = True , dest = 'bam_file', help="Alignment file of reads with the reference")
    

    ## general optional args
    optional_opts_group = arg_parser.add_argument_group('optional arguments')
    optional_opts_group.add_argument('--output_dir', required = False , dest = 'output_dir', help="Output directory, default is set to pred_vcf and filtering options applied")
    optional_opts_group.add_argument('--cpu', required = False , dest = 'cpu', default=1, help="CPU count to be applied for multiprocessing.")

    ## exome mode opts
    exome_opts_group = arg_parser.add_argument_group('exome reference mode')
    exome_opts_group.add_argument('--exome_bam', help="Alignment file of the exome with the reference")
    exome_opts_group.add_argument('--min_exome_depth', type = int, help="min exome depth required", default=10)
    

    ## restrict to high conf region opts
    restrict_regions_opts_group = arg_parser.add_argument_group("restrict evaluation regions on genome")
    restrict_regions_opts_group.add_argument("--restrict_regions_bed", help="bed file containing regions to restrict analysis to. (ie. high confidence regions)")

    ## var annotation options used later for filtering or analysis.
    var_annotation_opts_group = arg_parser.add_argument_group("variant annotations to include in variants table")
    var_annotation_opts_group.add_argument('--rna_editing', dest = 'rna_editing', help="input rna-editing vcf file")
    
    ## variant filtering options:
    var_filtering_opts_group = arg_parser.add_argument_group("variant filtering parameters")
    var_filtering_opts_group.add_argument('--min_RNAseq_Depth',  type = int , help = 'Depth value for SNP analysis \n', default=1)
    var_filtering_opts_group.add_argument("--indels_only", action='store_true', default=False, help='indels_only')
    var_filtering_opts_group.add_argument("--snvs_only", action='store_true', default=False, help='snvs_only')
    
    
    # Parse the arguments given 
    args = arg_parser.parse_args()
    

    ########
    ## Begin

    # set the output directory 
    # check if the output directory was given 
    # if not given as an argument, then create one 
    output_dir = args.output_dir
    if not output_dir:
        output_dir = create_output_dir_name(args)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ## write command for later auditing
    audit_command_info(output_dir, sys.argv)

    ## set up pipeliner
    checkpoints_dir = os.path.join(output_dir, "__chkpts")
    pipeliner = Pipeliner(checkpoints_dir)
    

    ############################################
    ''' 
    Aggregate data from vcfs and bam to table

    script used: *aggregate_variants_n_depth_to_table.py*
        input
            pred_vcf  : VCF of interest  
            truth_vcf : Reference VCF
            pred_bam  : BAM file used to  create the VCF of interest
            out       : Output file name with  extension ".aggregated.dat"
            tmpdir    : temp file directory 
            cpu.      : CPU count for multiprocessing
            restrict_regions_bed : Restrict the VCF to regions in the given bed file 
            exome_bam : If exome data is used 
            rna_editing: file holding RNA editing sites for removal 

    '''


    score_files = list()
    accuracy_files = list()

    # run on each vcf given as input. 
    for pred_vcf in args.pred_vcf:

        print("-PROCESSING {}".format(pred_vcf))
        
        pred_vcf_basename = os.path.basename(pred_vcf)

        aggregated_data_table_filename = os.path.join(output_dir, pred_vcf_basename + ".aggregated.dat")

        cmdstr = " ".join([ os.path.join(UTILDIR, "aggregate_variants_n_depth_to_table.py"),
                            "--pred_vcf {}".format(pred_vcf),
                            "--truth_vcf {}".format(args.truth_vcf),
                            "--pred_bam {}".format(args.bam_file),
                            "--out {}".format(aggregated_data_table_filename),
                            "--tmpdir {}".format(output_dir),
                            "--cpu {}".format(args.cpu)  ])
        
        if args.restrict_regions_bed:
            cmdstr += " --restrict_regions_bed {} ".format(args.restrict_regions_bed)

        if args.exome_bam:
            cmdstr += " --exome_bam {} ".format(args.exome_bam)

        if args.rna_editing:
            cmdstr += " --rna_editing {} ".format(args.rna_editing)


        pipeliner.add_commands([Command(cmdstr, "aggregated_{}.ok".format(pred_vcf_basename))])

        ################################
        '''
        Initial filtering of variants by attributes 
            1. Filters variants by RNA-seq depth 
            2. Filters variants by Exome depth 
            3. restrict to snvs or indels (optional)
            4. Remove variants found in rna-editing database 

        script used: *filterVariantTableByAttributes.py*
            input
                input_tabel    : 
                output_tabel   :

        '''

        filtered_data_file = aggregated_data_table_filename + ".filtered"

        cmdstr = " ".join([ os.path.join(UTILDIR, "filterVariantTableByAttributes.py"),
                            "--input_table {}".format(aggregated_data_table_filename), 
                            "--output_table {}".format(filtered_data_file) ])

        if args.snvs_only:
            cmdstr += " --snvs_only "
        elif args.indels_only:
            cmdstr += " --indels_only "
        

        if args.rna_editing:
            cmdstr += " --remove_rna_editing " 

        if args.exome_bam:
            cmdstr += " --min_exome_depth {} ".format(args.min_exome_depth)

        pipeliner.add_commands([Command(cmdstr, "filtered_table_{}.ok".format(pred_vcf_basename))])


        ###################
        ## score TP, FP, FN

        scored_variants_file = filtered_data_file + ".scored"

        cmdstr = " ".join([ os.path.join(UTILDIR, "classify_TP_FP_FN.py"),
                            "--input_table {}".format(filtered_data_file),
                            "--output_table {}".format(scored_variants_file) ])

        pipeliner.add_commands([Command(cmdstr, "scored_variants_{}.ok".format(pred_vcf_basename))])

        # tracking some output files for summaries below.
        score_files.append(scored_variants_file)
        accuracy_stats_file = scored_variants_file + ".accuracy_stats"  ## created by EvaluateAccuracy.py below for each input file.
        accuracy_files.append(accuracy_stats_file)

    ####################
    ## evaluate accuracy
    #ok_files = os.path.join(output_dir,'__chkpts','eval_accuracy_stats_*')

    #for filename in glob.glob(ok_files):
    #    os.remove(filename) 
    
    cmdstr = " ".join([ os.path.join(UTILDIR, "EvaluateAccuracy.py"),
                        "--input_tables {}".format(" ".join(score_files)),
                        "--output_dir {}".format(output_dir) ])
    pipeliner.add_commands([Command(cmdstr, "eval_accuracy_stats_{}.ok".format(time.time()))])

    #pipeliner.add_commands([Command(cmdstr, "eval_accuracy_stats_{}.ok".format(pred_vcf_basename))])

    ##################################################################
    ## always regenerate these plots - timepoints added to checkpoints.
    
    ## generate F1 score plots
    cmdstr = " ".join([ os.path.join(UTILDIR, "Accuracy_statistics.py"),
                        "--input_accuracy_table {}".format(" ".join(accuracy_files)),
                        "--output_dir {}".format(output_dir) ])
    pipeliner.add_commands([Command(cmdstr, "Accuracy_statistics_plots.{}.ok".format(time.time()))])

    ## generate PR and ROC plots

    cmdstr = " ".join([ os.path.join(UTILDIR, "plot_PR_n_ROC.py"),
                        "--input_accuracy_table {}".format(" ".join(accuracy_files)),
                        "--output_dir {}".format(output_dir) ])

    pipeliner.add_commands([Command(cmdstr, "PR_ROC_plots.{}.ok".format(time.time()))])

    ## generate SNP freqs plot
    if not args.indels_only:
        cmdstr = " ".join([ os.path.join(UTILDIR, "plot_SNP_freqs_barplot.py"),
                            "--input_scored_table {}".format(" ".join(score_files)),
                            "--output_dir {}".format(output_dir) ])

        pipeliner.add_commands([Command(cmdstr, "snp_freq_barplots.{}.ok".format(time.time()))])

    
    pipeliner.run()

    
if __name__ == '__main__':
    main()

