#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: marek svoboda
"""

import argparse
from gffutils import create_db
from RNAseqAnalysis import getBaselineData
from BAMfilter import scanBAMfilter


def createDB_(args):
    create_db(args.data,
              args.dbfn,
              merge_strategy = args.merge_strategy,
              verbose = args.verbose)
    
def createTRANS_(args):
    getBaselineData(args.out_transBaselineData,
                    args.out_db,
                    args.bamfile,
                    includeIntrons = args.includeIntrons,
                    weightedCov = args.weightedCov)

def BAMfilter_(args):
    scanBAMfilter(args.covLen,
                  args.minSNRlen,
                  args.bamfile,
                  args.fastafile,
                  args.out_transBaselineData,
                  out_db = args.out_db,
                  base = args.base,
                  mism = args.mism,
                  out_bamfile = args.out_bamfile,
                  tROC = args.tROC,
                  includeIntrons = args.includeIntrons,
                  cbFile = args.cbFile,
                  out_cbFile = args.out_cbFile,
                  weightedCov = args.weightedCov,
                  verbose = args.verbose,
                  nThreads = args.nThreads)


if __name__ == '__main__':
    
    # Create the top-level parser and enable addition of subparsers
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        help = 'Select one of these functions to run:')
    
    
    # Parser for gffutils.create_db
    parser_createDB = subparsers.add_parser(
        'createDB',
        help = 'Creates a reference annotation database from a GFF/GTF file',
        description = '''
        Wrapper around the create_db function from the gffutils package using
        predefined parameters to ensure that the database is created properly
        for the puposes of the BAMfilter.
        ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_createDB.add_argument(
        'data',
        metavar = 'GTF/GFF_FILE',
        help = 'Location of the reference GTF/GFF file')
    parser_createDB.add_argument(
        'dbfn',
        metavar = 'DB_FILE',
        help = 'Location of the output annotation database')
    parser_createDB.add_argument(
        '-m', '--merge_strategy',
        default = 'create_unique',
        help = argparse.SUPPRESS)
    parser_createDB.add_argument(
        '-v', '--verbose',
        action = 'store_true',
        help = 'Extra messages are logged.')
    parser_createDB.set_defaults(func = createDB_)
    
    
    # Parser for getBaselineData
    parser_createTRANS = subparsers.add_parser(
        'createTRANS',
        help = 'Creates a cache file that stores all expressed transcripts',
        description = '''
        Creates a cache file with all transcripts that are expressed as per
        the associated BAM file.
        ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_createTRANS.add_argument(
        'out_db',
        metavar = 'DB_FILE',
        help = 'Location of the reference annotation database')
    parser_createTRANS.add_argument(
        'bamfile',
        metavar = 'BAM_FILE',
        help = 'Location of the sorted and indexed BAM file')
    parser_createTRANS.add_argument(
        'out_transBaselineData',
        metavar = 'TRANS_FILE',
        help = 'Location of the output cache file')
    parser_createTRANS.add_argument(
        '-i', '--introns',
        action = 'store_true',
        dest = 'includeIntrons',
        help = 'Instructs to include intronic coverage')
    parser_createTRANS.add_argument(
        '-w', '--weightedCov',
        action = 'store_false',
        help = argparse.SUPPRESS)
    parser_createTRANS.set_defaults(func = createTRANS_)
    
    
    # Parser for the BAMfilter
    parser_BAMfilter = subparsers.add_parser(
        'BAMfilter',
        help = 'Filters a BAM file to remove internally primed alignments',
        description = '''
        Filters a sorted and indexed BAM file to remove sparse alignments that
        likely resulted from internal priming.
        ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_BAMfilter.add_argument(
        'covLen',
        type = int,
        metavar = 'COVLEN',
        help = 'Maximal coverage distance')
    parser_BAMfilter.add_argument(
        'minSNRlen',
        type = int,
        metavar = 'MINSNRLEN',
        help = 'Minimal SNR length')
    parser_BAMfilter.add_argument(
        'bamfile',
        metavar = 'BAM_FILE',
        help = 'Location of the sorted and indexed BAM file')
    parser_BAMfilter.add_argument(
        'fastafile',
        metavar = 'FASTA_FILE',
        help = 'Location of the reference FASTA file')
    parser_BAMfilter.add_argument(
        'out_transBaselineData',
        metavar = 'TRANS_FILE',
        help = 'Location of a cache file that stores expressed transcripts; if'
        ' it does not already exist, it will be created')
    parser_BAMfilter.add_argument(
        '-d', '--out_db',
        metavar = 'DB_FILE',
        help = 'Location of the reference annotation database; not necessary '
        'if TRANS_FILE has already been created.')
    parser_BAMfilter.add_argument(
        '-b', '--base',
        type = str,
        choices = ['A', 'T', 'C', 'G', 'N'],
        default = 'A',
        metavar = 'N',
        help = 'DNA base constituting SNRs to be searched for')
    parser_BAMfilter.add_argument(
        '-m', '--mismatches',
        type = int,
        default = 0,
        dest = 'mism',
        metavar = 'MISM',
        help = 'The maximum number of mismatches allowed in an SNR')
    parser_BAMfilter.add_argument(
        '-o', '--out_bamfile',
        metavar = 'OUT_BAM_FILE',
        help = 'Location of the resulting filtered bamfile; if not provided, '
        'the name of the input bamfile name is used with ".filtered.bam" '
        'appended.')
    parser_BAMfilter.add_argument(
        '-t', '--tROC',
        type = dict,
        default = {},
        help = argparse.SUPPRESS)
    parser_BAMfilter.add_argument(
        '-i', '--introns',
        action = 'store_true',
        dest = 'includeIntrons',
        help = 'Instructs to include intronic coverage')
    parser_BAMfilter.add_argument(
        '-c', '--cbFile',
        metavar = 'CB_FILE',
        help = 'If the scumi package is being used for alignment counting, '
        'provide the location of the file with cell barcode count (output of '
        'scumi merge_fastq)')
    parser_BAMfilter.add_argument(
        '-e', '--out_cbFile',
        metavar = 'OUT_CB_FILE',
        help = 'If the scumi package is being used for alignment counting, '
        'provide the location for the modified cell barcode count file to be '
        'used as input for scumi count_umi; if not provided in the presence of'
        ' cbFile, cbFile name is used with ".filtered.tsv" appended')
    parser_BAMfilter.add_argument(
        '-w', '--weightedCov',
        action = 'store_false',
        help = argparse.SUPPRESS)
    parser_BAMfilter.add_argument(
        '-v', '--verbose',
        action = 'store_true',
        help = 'Extra messages are logged.')
    parser_BAMfilter.add_argument(
        '-p', '--processes',
        type = int,
        dest = 'nThreads',
        metavar = 'P',
        help = '''
        Set the maximum number of processes to use besides the main process.
        If not provided, serial processing is used. Note that parallel
        processing will produce the following warning for each temp file, which
        can be safely ignored:
            
        "[E::idx_find_and_load] Could not retrieve index file for <file>"
        ''')
    parser_BAMfilter.set_defaults(func = BAMfilter_)    
    
    
    # Parse the arguments and populate the namespace; execute the code
    args = parser.parse_args()
    args.func(args)

