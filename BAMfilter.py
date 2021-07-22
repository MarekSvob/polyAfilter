#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:11:28 2020

@author: marek
"""
import os
import gc
import pysam
import logging
import multiprocessing
import multiprocessing.pool

import numpy as np
import pandas as pd
import regex as re
from functools import partial
from scumi import scumi as scu
from Bio import SeqIO
from random import seed, random
from collections import defaultdict

from RNAseqAnalysis import getBaselineData, getTransEndSensSpec, \
    sortSNRsByLenStrdRef
from SNRanalysis import flattenIntervals, getOverlaps
from SNRdetection import loadPKL


logging.basicConfig(level = logging.INFO,
                    format = '%(asctime)s - %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def cbFileFilter(toRemove, cbFile, out_cbFile, verbose):
    """Function that subtracts the appropriate numbers from the CB counts,
    including the T frequencies in UMIs. Note that where related, parts of this
    function's code were adapted from the scumi module
    [https://bitbucket.org/jerry00/scumi-dev/src/master/].

    Parameters
    ----------
    toRemove : (set)
        The set of bamfile alignments to remove.
    cbFile : (str), optional
        Location of the scumi cell barcode count file.
    out_cbFile : (str), optional
        Location of a new cell barcode count file. If None, the cbFile name
        is used with '.filtered.tsv' appended.
    verbose : (bool)
        Indicates whether extra messages are logged.

    Returns
    -------
    None.
    """
    
    if verbose:
        logger.info('Counting UMIs per cell to be removed from the CB file...')
        
    if len(toRemove):
        # Note: the following code has been adapted from the scumi module
        # Get one read as an example
        firstRead = next(iter(toRemove))
        # Determine which barcodes are present
        barcodes = set()
        for barcode in ['CB_', 'UB_']:
            if barcode in firstRead.query_name:
                barcodes.add(barcode)
        # Construct the barcode parser based on the barcodes present
        barcodeParser = '.*'
        if 'CB_' in barcodes:
            barcodeParser += ':CB_(?P<CB>[A-Z\-]+)'
        if 'UB_' in barcodes:
            barcodeParser += ':UB_(?P<UB>[A-Z\-]+)'
        # Skip counting barcodes if none present
        if barcodeParser == '.*':
            logger.error('Error: no cell barcodes or UMIs.')
            return
        barcodeParser += ':*'
        barcodeParser = re.compile(barcodeParser)
        
        # Initiate the counter of UMIs per cell, including T frequency
        cbDF = pd.read_csv(cbFile, sep = '\t', index_col = 0, header = 0)
        cbCounter = defaultdict(
            partial(np.zeros, shape = cbDF.shape[1], dtype = np.uint32))
        
        for read in toRemove:
            # Extract the respective barcodes & count them for each read
            match = barcodeParser.match(read.query_name)
            cb = scu._extract_tag(match, 'CB')
            umi = scu._extract_tag(match, 'UB')
            cbCounter[cb] += [x == 'T' for x in 'T' + umi]
        
        # Create a df from the results
        rmDF = pd.DataFrame.from_dict(cbCounter, orient = 'index')
        rmDF.index.name = cbDF.index.name
        rmDF.columns = cbDF.columns
        # Subtract the rmDF from the cbDF & save
        diffDF = cbDF.subtract(rmDF, fill_value = 0)
        # Remove rows with 0 in cb_counts or with any negative numbers anywhere
        diffDF = diffDF[diffDF['cb_count'] > 0]
        diffDF = diffDF[(diffDF >= 0).all(1)].astype(np.uint32)
        # Sort by cb_count
        diffDF.sort_values(by = 'cb_count', ascending = False, inplace = True)
        # Report the results
        nCells = rmDF.shape[0]
        nUMIs = sum(rmDF['cb_count'])
    else:
        diffDF = pd.read_csv(cbFile, sep = '\t', index_col = 0, header = 0)
        nCells = 0
        nUMIs = 0
        
    # Save the file
    if out_cbFile is None:
        out_cbFile = cbFile + '.filtered.tsv'
    diffDF.to_csv(out_cbFile, sep = '\t')

    if verbose:
        logger.info(f'{nUMIs:,d} reads across {nCells:,d} cell barcodes '
                    'removed from the CB file.')


def getAlignmentsToRemove(strd, refName, eachTransStart):
    """Helper function to run on each thread for a multithreaded alignment
    identification and write as temporary files to be merged later. Note that
    some of the variables used are set to be global by child_initialize().
    
    Parameters
    ----------
    strd : (bool)
        Strand of the reference; +(True) / -(False)
    refName : (str)
        Name of the reference to scan.
    eachTransStart : (tuple)
        Tuple of trans start pieces such that
        (((start, end), (start, end), ...), ((start, end), (start, end), ...))
    
    Returns
    -------
    None.
    """
    
    if verbose:
        logger.info('Identifying alignments to be removed on reference '
                    f'{refName}{"+" if strd else "-"}...')
    flatStarts = []
    for tStart in eachTransStart:
        flatStarts.extend(tStart)
    flatStarts = flattenIntervals(flatStarts)
    # Connect to the bam file for each thread separately
    bam = pysam.AlignmentFile(bamfile, 'rb')
    refLen = bam.get_reference_length(refName)
    
    # Flatten all the expressed transcript starts on this strd/ref
    # minSNRlen = None means that all non-end coverage is removed,
    #  regardless of SNR pieces
    if minSNRlen is None:
        SNRpieceOverlaps = flatStarts
    else:
        # Flatten all the relevant SNR pieces
        SNRpieces = []
        for length, SNRsByStrdRef in SNRsByLenStrdRef.items():
            if length >= minSNRlen:
                for SNR in SNRsByStrdRef[strd][refName]:
                    if strd:
                        start = max(0, SNR.start - covLen)
                        end = SNR.start
                    else:
                        start = SNR.end
                        end = min(refLen, SNR.end + covLen)
                    # Make sure that this is a non-0 interval
                    #  (when SNR is at the beginning of the ref)
                    if start != end:
                        SNRpieces.append((start, end))
        SNRpieces = flattenIntervals(SNRpieces)
        # Find overlaps between SNRpieces & tStarts
        SNRpieceOverlaps = getOverlaps(flatStarts, SNRpieces)
    
    # Go over all the pieces to fetch the alignments & write them into a temp
    tempFile = f'{out_bamfile}_{strd}_{refName}.bam'
    counter = 0
    bamTEMP = pysam.AlignmentFile(tempFile, 'wb', template = bam)
    for start, end in SNRpieceOverlaps:
        for alignment in bam.fetch(contig = refName, start = start, stop = end):
            # Ensure the alignment is on the correct strd before adding
            if alignment.is_reverse != strd:
                # Make sure that alignments filtered indeed map onto the area
                #  identified, not just overlap it with their start-end range
                #  (in case of splicing)
                mapping = getOverlaps(alignment.get_blocks(), [(start, end)])
                if mapping != []:
                    bamTEMP.write(alignment)
                    counter += 1
    bam.close()
    bamTEMP.close()
    
    if verbose:
        logger.info(f'Identified {counter:,d} alignments to be removed on '
                    f'reference {refName}{"+" if strd else "-"} and saved to '
                    f'{tempFile}.')


def child_initialize(_SNRsByLenStrdRef, _covLen, _minSNRlen, _bamfile,
                     _verbose, _out_bamfile):
    """Helper function to initialize and share variables between parallel
    processes. Taken from
    https://stackoverflow.com/questions/25825995/
    python-multiprocessing-only-one-process-is-running
    
    Parameters
    ----------
    _SNRsByLenStrdRef : (dict)
        { length : { strd : { ref : [ SNRs ] } } } or, alternatively when
        sortedSNRs = False: { length : [ SNRs ] }
    _covLen : (int)
        The assumed distance between the beginning of alignment coverage and
        the priming event.
    _minSNRlen : (int, None)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    _bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    _verbose : (bool)
        Indicates whether extra messages are logged. The default is False.
    _out_bamfile : (str)
        The location of the target bamfile to save temporary files in the same
        folder.

    Returns
    -------
    None.
    """
    global SNRsByLenStrdRef, covLen, minSNRlen, bamfile, verbose, out_bamfile
    
    SNRsByLenStrdRef = _SNRsByLenStrdRef
    covLen = _covLen
    minSNRlen = _minSNRlen
    bamfile = _bamfile
    verbose = _verbose
    out_bamfile = _out_bamfile
    

class AngelProcess(multiprocessing.Process):
    """Class for embedded parallel processes adapted from
    https://stackoverflow.com/questions/6974695/
    python-process-pool-non-daemonic
    """
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


class MyPool(multiprocessing.pool.Pool):
    """Class for embedded parallel processes adapted from
    https://stackoverflow.com/questions/6974695/
    python-process-pool-non-daemonic
    """
    Process = AngelProcess
    

def parallel_wrapper(covLen, minSNRlen, bamfile, out_SNRsByLenStrdRef,
                     out_transBaselineData, out_db, out_bamfile, tROC,
                     includeIntrons, sortedSNRs, weightedCov, verbose,
                     nThreads):
    """Helper function to hold SNRsByLenStrdRef in RAM only temporarily.
    
    Parameters
    ----------
    eachTransStartByStrdRef : (dict)
        { strd : { refName : [ ( (start, end), ... ) ] } } for start exon
        pieces for each covered transcript.
    covLen : (int)
        The assumed distance between the beginning of alignment coverage and
        the priming event.
    minSNRlen : (int, None)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    out_SNRsByLenStrdRef : (str)
        Location of the SNRS to be loaded, such that:
        { length : { strd : { ref : [ SNRs ] } } } or, alternatively when
        sortedSNRs = False: { length : [ SNRs ] }
    out_transBaselineData : (str)
        Path to the saved baseline data, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    out_bamfile : (str, None)
        Location of the resulting filtered bamfile. If None, the name of the
        input bamfile name is used with '.filtered.bam' appended. The default
        is None.
    tROC : (dict)
        { endLen : ( TP, FN, TN, FP, eachTransStartByStrdRef ) }.
    includeIntrons : (bool)
        Indicates whether to consider intron coverage.
    sortedSNRs : (bool)
        Indicates whether the SNRs are already sorted by length, strd, and ref.
    weightedCov : (bool)
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Alternatively, the sensitivity/specificity
        relationship is purely binary (flat) and the covered bases (1) have the
        same weight as non-covered ones (0).
    verbose : (bool)
        Indicates whether extra messages are logged.
    nThreads : (int),
        Set the maximum number of threads to use. If None, serial processing is
        used but not enforced (underlying processes may use multiple threads).

    Returns
    -------
    strdRefs : list
        [ (strd, ref) ]
    """
    
    # Get transcript starts if N/A for this specific length
    if covLen not in tROC:
        BLdata = getBaselineData(out_transBaselineData, out_db, bamfile,
                                 includeIntrons, weightedCov = weightedCov)
        tROC[covLen] = getTransEndSensSpec(
            covLen, bamfile, BLdata, includeIntrons, getSensSpec = False,
            weightedCov = weightedCov)    
    # Extract the transcript starts for this coverage length
    eachTransStartByStrdRef = tROC[covLen][4]
    # If relevant, load and pre-sort SNRs by len, strd & ref
    if minSNRlen is None:
        SNRsByLenStrdRef = None
    else:
        logger.info(f'Loading SNRs from {out_SNRsByLenStrdRef}...')
        SNRsByLenStrdRef = loadPKL(out_SNRsByLenStrdRef)
        if not sortedSNRs:
            SNRsByLenStrdRef = sortSNRsByLenStrdRef(SNRsByLenStrdRef)        
    
    logger.info(f'Identifying alignments {covLen:,d} bp upstream of non-'
                f'terminal SNR{minSNRlen}+ to be removed across {nThreads} '
                'parallel processes and writing temporary files...')
    # Create a pool of processes with shared variables
    pool = multiprocessing.Pool(processes = nThreads,
                                initializer = child_initialize,
                                initargs = (SNRsByLenStrdRef, covLen,
                                            minSNRlen, bamfile, verbose,
                                            out_bamfile))
    # Identify the alignments to be removed by strand / ref
    strdRefs = []
    for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
        for refName, eachTransStart in eachTransStartByRef.items():
            strdRefs.append((strd, refName))
            pool.apply_async(func = getAlignmentsToRemove,
                             args = (strd, refName, eachTransStart))
    # Close the pool
    pool.close()
    # Join the processes
    pool.join()
    
    return strdRefs
    

def BAMfilter(covLen, minSNRlen, bamfile, out_SNRsByLenStrdRef,
              out_transBaselineData, out_db, out_bamfile = None, tROC = {},
              includeIntrons = False, cbFile = None, out_cbFile = None,
              sortedSNRs = True, weightedCov = True, verbose = False,
              nThreads = None):
    """Function that filters an indexed BAM file to remove non-canonical
    alignments that likely resulted from poly(dT) priming onto genomically
    encoded polyA single nucleotide repeats (SNRs), as opposed to mRNA polyA
    tails. Note that parallel processing will produce the following warning
    for each temp file, which can be safely ignored:
        "[E::idx_find_and_load] Could not retrieve index file for <file>"

    Parameters
    ----------
    covLen : (int)
        The assumed distance between the beginning of alignment coverage and
        the priming event.
    minSNRlen : (int, None)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    out_SNRsByLenStrdRef : (str)
        Location of the SNRS to be loaded, such that:
        { length : { strd : { ref : [ SNRs ] } } } or, alternatively when
        sortedSNRs = False: { length : [ SNRs ] }
    out_transBaselineData : (str)
        Path to the saved baseline data, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    out_bamfile : (str, None)
        Location of the resulting filtered bamfile. If None, the name of the
        input bamfile name is used with '.filtered.bam' appended. The default
        is None.
    tROC : (dict), optional
        { endLen : ( TP, FN, TN, FP, eachTransStartByStrdRef ) }. The default
        is {}.
    includeIntrons : (bool), optional
        Indicates whether to consider intron coverage. The default is False.
    cbFile : (str, None), optional
        Location of the scumi cell barcode count file, if any. The default is
        None.
    out_cbFile : (str, None), optional
        Location of a new scumi cell barcode count file, if any. If None, the
        cbFile name is used with '.filtered.tsv' appended. The default is None.
    sortedSNRs : (bool), optional
        Indicates whether the SNRs are already sorted by length, strd, and ref.
        The default is True.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Only included for the purposes of using
        cached baseline files; otherwise has no impact on the results of this
        function. The default is True.
    verbose : (bool), optional
        Indicates whether extra messages are logged. The default is False.
    nThreads : (int), optional
        Set the maximum number of threads to use besides the main thread.
        If None, serial processing is used. The default is None.

    Returns
    -------
    None.
    """    
    # Initiate the correct output file name
    if out_bamfile is None:
        out_bamfile = bamfile + '.filtered.bam'
    # If the BAM file already exists, do not overwrite
    if os.path.isfile(out_bamfile):
        logger.info(f'The filtered BAM file "{out_bamfile}" already exists.')
        return
        
    # Initialize the set of alignments to be removed
    toRemove = set()
    
    if nThreads is None:
        # Get transcript starts if N/A for this specific length
        if covLen not in tROC:
            BLdata = getBaselineData(out_transBaselineData, out_db, bamfile,
                                     includeIntrons, weightedCov = weightedCov)
            tROC[covLen] = getTransEndSensSpec(
                covLen, bamfile, BLdata, includeIntrons, getSensSpec = False,
                weightedCov = weightedCov)    
        # Extract the transcript starts for this coverage length
        eachTransStartByStrdRef = tROC[covLen][4]
        # If relevant, load and pre-sort SNRs by len, strd & ref
        if minSNRlen is not None:
            logger.info(f'Loading SNRs from {out_SNRsByLenStrdRef}...')
            SNRsByLenStrdRef = loadPKL(out_SNRsByLenStrdRef)
            if not sortedSNRs:
                SNRsByLenStrdRef = sortSNRsByLenStrdRef(SNRsByLenStrdRef)
        
        logger.info(f'Identifying alignments {covLen:,d} bp upstream of non-'
                    f'terminal SNR{minSNRlen}+ to be removed in 1 serial '
                    'process...')
        # Connect to the bam file
        bam = pysam.AlignmentFile(bamfile, 'rb')
        
        SNRpieceOverlapsByStrdRef = defaultdict(lambda: defaultdict(dict))
        
        # Go over each strand and ref separately to obtain overlaps with SNRs
        for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
            for refName, eachTransStart in eachTransStartByRef.items():
                if verbose:
                    logger.info('Getting overlaps of SNRs on transcripts on re'
                                f'ference {refName}{"+" if strd else "-"}...')
                if not strd:
                    refLen = bam.get_reference_length(refName)
                # Flatten all the expressed transcript starts on this strd/ref
                flatStarts = []
                for tStart in eachTransStart:
                    flatStarts.extend(tStart)
                flatStarts = flattenIntervals(flatStarts)
                # minSNRlen = None means that all non-end coverage is removed,
                #  regardless of SNR pieces
                if minSNRlen is None:
                    SNRpieceOverlapsByStrdRef[strd][refName] = flatStarts
                else:
                    # Flatten all the relevant SNR pieces
                    SNRpieces = []
                    for length, SNRsByStrdRef in SNRsByLenStrdRef.items():
                        if length >= minSNRlen:
                            for SNR in SNRsByStrdRef[strd][refName]:
                                if strd:
                                    start = max(0, SNR.start - covLen)
                                    end = SNR.start
                                else:
                                    start = SNR.end
                                    end = min(refLen, SNR.end + covLen)
                                # Make sure that this is a non-0 interval
                                #  (when SNR is at the beginning of the ref)
                                if start != end:
                                    SNRpieces.append((start, end))
                    SNRpieces = flattenIntervals(SNRpieces)
                    # Find overlaps between SNRpieces & tStarts
                    SNRpieceOverlapsByStrdRef[strd][refName] = getOverlaps(
                        flatStarts, SNRpieces)
                    
        # Remove SNRsByLenStrdRef to free up RAM for the alignments
        if minSNRlen is not None:
            del SNRsByLenStrdRef
            gc.collect()
        
        # Go over all the pieces, by strd & ref, to fetch the alignments
        for strd, SNRpieceOverlapsByRef in SNRpieceOverlapsByStrdRef.items():
            for refName, SNRpieceOverlaps in SNRpieceOverlapsByRef.items():
                if verbose:
                    logger.info('Identifying alignments to be removed on refe'
                                f'rence {refName}{"+" if strd else "-"}...')
                for start, end in SNRpieceOverlaps:
                    for alignment in bam.fetch(
                            contig = refName, start = start, stop = end):
                        # Ensure the alignment is on the correct strd
                        if alignment.is_reverse != strd:
                            # Make sure that filtered alignments indeed map to
                            #  the area identified, not just overlap it with
                            #  their start-end range (in case of splicing)
                            mapping = getOverlaps(alignment.get_blocks(),
                                                  [(start, end)])
                            if mapping != []:
                                toRemove.add(alignment)
        bam.close()
    
    else:
        # Set the number of threads for NumExpr (used by pandas & numpy)
        os.environ['NUMEXPR_MAX_THREADS'] = str(nThreads)
        # Start 1 parallel process to load SNRsByLenStrdRef only temporarily
        pool = MyPool(processes = 1)
        result = pool.apply_async(
            func = parallel_wrapper,
            args = (covLen, minSNRlen, bamfile, out_SNRsByLenStrdRef,
                    out_transBaselineData, out_db, out_bamfile, tROC,
                    includeIntrons, sortedSNRs, weightedCov, verbose,
                    nThreads))
        # Close the pool
        pool.close()
        # Join the processes
        pool.join()
        # Obtain the list of strd-ref pairs processed; also helps uncover any
        #  exceptions that came up during that process
        strdRefs = result.get()
        
        # Merge the temporary files into a single set & remove the files
        logger.info(f'Merging {len(strdRefs)} temporary files into memory and '
                    'deleting...')
        for strand, ref in strdRefs:
            tempfile = f'{out_bamfile}_{strand}_{ref}.bam'
            bamTEMP = pysam.AlignmentFile(tempfile, 'rb')
            for alignment in bamTEMP.fetch(until_eof = True):
                toRemove.add(alignment)
            bamTEMP.close()
            os.remove(tempfile)
    
    toRemoveN = len(toRemove)
    # Filter the cbFile, if any
    if cbFile:
        cbFileFilter(toRemove, cbFile, out_cbFile, verbose)   
    
    # Create the bamfile and add the reads not in the toRemove set
    # Manually reopen the bam file to avoid multiple iterators issues
    bamIN = pysam.AlignmentFile(bamfile, 'rb')   
    # Get the total number of reads in the indexed bam file
    total = bamIN.mapped + bamIN.unmapped
    logger.info(f'Writing the filtered BAM file, to exclude {toRemoveN:,d} '
                f'alignments out of {total:,d} total...')
    # Write the new bam file
    bamOUT = pysam.AlignmentFile(out_bamfile, 'wb', template = bamIN)
    nAll = 0
    included = 0
    for alignment in bamIN.fetch(until_eof = True):
        nAll += 1
        if alignment not in toRemove:
            bamOUT.write(alignment)
            included += 1
            
        if not nAll % 10000000 and nAll and verbose:
            logger.info(f'Processed {nAll:,d} input alignments, of which '
                        f'{nAll-included:,d} were excluded...')
            
    # Close the files
    bamIN.close()
    bamOUT.close()
    logger.info(f'A filtered BAM file with {included:,d}/{nAll:,d} alignments '
                f'has been created. [Excluded {nAll-included:,d} alignments '
                f'{covLen:,d} bp upstream of non-terminal SNR{minSNRlen}+.]')


def propBAMfilter(remProp, bamfile, setSeed = 1, out_bamfile = None,
                  cbFile = None, out_cbFile = None, verbose = False):
    """Function that RANDOMLY filters an indexed BAM file to remove a given
    proportion of total alignments (note that only *mapped* alignments will be
    removed).

    Parameters
    ----------
    remProp : (float)
        The proportion of *all* alignments to be removed. (However, only
        mapped alignments will be removed.)
    bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    setSeed : (int), optional
        Ensures reproducibility of pseudorandom functions. The default is 1.
    out_bamfile : (str, None)
        Location of the resulting filtered bamfile. If None, the name of the
        input bamfile name is used with '.filtered.bam' appended. The default
        is None.
    cbFile : (str, None), optional
        Location of the scumi cell barcode count file, if any. The default is
        None.
    out_cbFile : (str, None), optional
        Location of a new scumi cell barcode count file, if any. If None, the
        cbFile name is used with '.filtered.tsv' appended. The default is None.
    verbose : (bool), optional
        Indicates whether extra messages are logged. The default is False.

    Returns
    -------
    None.
    """    
    # Initiate the correct output file name
    if out_bamfile is None:
        out_bamfile = bamfile + '.filtered.bam'
    # If the BAM file already exists, do not overwrite
    if os.path.isfile(out_bamfile):
        logger.info(f'The filtered BAM file "{out_bamfile}" already exists.')
        return
    
    # Set the seed for reproducibility
    seed(a = setSeed)
    # Initiate the set of removed alignments (only used for cbFile filtering)
    removed = set()
    
    # Create the bamfile and add the reads randomly kept
    bamIN = pysam.AlignmentFile(bamfile, 'rb')
    # Obtain the probability for the *mapped* (as opposed to all) alignments
    #  to be kept
    # Get the total # of alignments in the BAM file; more details:
    #  https://github.com/pysam-developers/pysam/issues/968
    total = bamIN.mapped + bamIN.unmapped
    # Calculate the # of alignmens to be removed
    nAligToRemove = int(remProp * total)
    # Calculate the proportion of mapped alignments to be kept
    keepMapProb = 1 - nAligToRemove / (total - bamIN.nocoordinate)

    logger.info(f'Writing the filtered BAM file, to exclude approximately '
                f'{nAligToRemove:,d} mapped alignments out of {total:,d} '
                'total...')
    # Write the new bam file, including all unmapped reads and a proportion of
    #  mapped reads
    bamOUT = pysam.AlignmentFile(out_bamfile, 'wb', template = bamIN)
    nAll = 0
    included = 0
    for alignment in bamIN.fetch(until_eof = True):
        nAll += 1
        if alignment.is_unmapped or random() < keepMapProb:
            bamOUT.write(alignment)
            included += 1
        elif cbFile:
            removed.add(alignment)
            
        if not nAll % 10000000 and nAll and verbose:
            logger.info(f'Processed {nAll:,d} input alignments, of which '
                        f'{nAll-included:,d} were excluded...')
            
    # Close the files
    bamIN.close()
    bamOUT.close()
    logger.info(f'A filtered BAM file with {included:,d}/{nAll:,d} alignments '
                f'has been created. [Randomly excluded {nAll-included:,d} '
                'alignments.]')

    # Filter the cbFile, if any
    if cbFile:
        cbFileFilter(removed, cbFile, out_cbFile, verbose)


def findSNRpieces(covLen, seq, strd, base = 'A', mism = 0, minSNRlen = 5):
    """A version of findSNRsWmisms() that detects SNR pieces directly on each
    individual strand and does not store any information about SNRs themselves.
    A universal detection of SNRs of maximum length with a given number of
    mismatches allowed in non-terminal positions. This function is meant to run
    once per processing thread. Conceptually, once the maximum # of mismatches
    is encountered after at least 1 block of bases, and the total sum of bp in
    the current block(s) is higher than minSNRlen, a new SNR is added.

    Parameters
    ----------
    covLen : (int)
        The length of the pieces to detect.
    seq : (str)
        The ACTG reference sequence.
    strd : (bool)
        The strand to scan for SNRs; T (+) or F (-).
    base : (str), optional
        DNA base constituting SNRs to be searched for: 'A', 'T', 'C', 'G', 'N'.
        The default is 'A'.
    mism : (int), optional
        The maximum number of mismatches allowed. The default is 0.
    minSNRlen : (int), optional
        The minimal length of the SNRs saved; limited to save storage space
        needed. The default is 5.
    
    Returns
    -------
    SNRpieces : (list)
        [ (start, end) ]
    """
    
    def mayAddPiece():
        """A helper function that adds an SNR if the length condition is met.
        Similar to the helper function in findSNRsWmisms() but adds SNR pieces
        instead of SNRs.
        """
        nonlocal SNRpieces, lastBlockSaved
        
        first = blocks[0][0]
        last = blocks[-1][-1]
        length = last - first
        # If minSNRlen has been met, add the SNR piece into the output list and
        #  mark the last block as saved
        if length >= minSNRlen:
            if strd:
                start = max(0, first - covLen)
                end = first
            else:
                start = last
                end = min(refLen, last + covLen)
            # Make sure that this is a non-0 interval (when SNR is at the
            #  beginning of the ref)
            if start != end:
                SNRpieces.append((start, end))
            
            lastBlockSaved = True
    
    # Handling of non-caps in the sequence
    capsDict = {'A':{'a', 'A'},
                'T':{'t', 'T'},
                'C':{'c', 'C'},
                'G':{'g', 'G'},
                'N':{'n', 'N'}}
    # Handling of complementary bases for '-' strd
    compDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    
    SNRpieces = []
    # If this scan of a '-' strand, switch to the complementary base
    #  and save the sequence length (irrelevant on '+' strand)
    if not strd:
        base = compDict[base]
        refLen = len(seq)
    
    # Allow for either caps to be included
    eitherCaps = capsDict[base]
    # Initiate the mismatch counter
    nMism = 0
    # Initiate the list of valid base blocks
    blocks = []
    # Initiate the indicator of being in a block
    blocked = False
    # Initiate the indicator of having saved an SNR piece with the last block
    lastBlockSaved = True
    # Scanning the sequence, while adding and removing blocks
    for i, bp in enumerate(seq):
        # If this is the base but block has not been started, start one
        #  (if currently in a block, just pass onto the next base - the two
        #  conditions below cannot be merged into one)
        if bp in eitherCaps:
            if not blocked:
                start = i
                blocked = True
        else:
            # If a block just ended, add it
            if blocked:
                blocks.append((start, i))
                lastBlockSaved = False
                blocked = False
            # Count the mismatch only if inside a potential SNR
            if blocks:
                nMism += 1
                # If mism has been exceeded, an SNR may be added;
                #  always remove the 1st block & decrease the mism
                if nMism > mism:
                    # Only save the SNR if the most recent block hasn't
                    #  been saved in an SNR yet (if it has, it implies
                    #  that this would not be the longest SNR possible)
                    if not lastBlockSaved:
                        mayAddPiece()
                    # Deduct all the mismatches before the 2nd block
                    if len(blocks) > 1:
                        nMism -= blocks[1][0] - blocks[0][-1]
                    else:
                        nMism = 0
                    # Remove the 1st block
                    blocks = blocks[1:]
    # Once the sequence scanning is done, the last SNR may be added
    #  even if the mism # has not been reached
    # If the sequence ended in a block, add it
    if blocked:
        blocks.append((start, i+1))
        lastBlockSaved = False
    # If the last block has not been saved, an SNR may be added
    if blocks and not lastBlockSaved:
        mayAddPiece()
    
    return SNRpieces


def scanForAlignmentsToRemove(covLen, refName, strd, minSNRlen, bamfile,
                              out_bamfile, fastafile, eachTransStart,
                              base = 'A', mism = 0, verbose = True):
    """Helper function to run on each thread for a multithreaded scanning and
    alignment identification and writing into temporary files to be merged
    later.
    
    Parameters
    ----------
    covLen : (int)
        The assumed distance between the beginning of alignment coverage and
        the priming event.
    refName : (str)
        Name of the reference to scan.
    strd : (bool)
        Strand of the reference; +(True) / -(False)
    minSNRlen : (int, None)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    out_bamfile : (str)
        Location of the resulting filtered bamfile.
    fastafile : (str)
        Location of the reference genome fasta file.
    eachTransStart : (tuple)
        Tuple of trans start pieces such that
        (((start, end), (start, end), ...), ((start, end), (start, end), ...))
    base : (str), optional
        DNA base constituting SNRs to be searched for: 'A', 'T', 'C', 'G', 'N'.
        The default is 'A'.
    mism : (int), optional
        The maximum number of mismatches in an SNR allowed. The default is 0.
    verbose : (bool), optional
        Indicates whether extra messages are logged. The default is False.
    
    Returns
    -------
    None.
    """        
    
    # Flatten all the expressed transcript starts on this strd/ref
    flatStarts = []
    for tStart in eachTransStart:
        flatStarts.extend(tStart)
    flatStarts = flattenIntervals(flatStarts)
    # minSNRlen = None means that all non-end coverage is removed,
    #  regardless of SNR pieces
    if minSNRlen is None:
        SNRpieces = flatStarts
    else:
        # Scan for and flatten all the relevant SNR pieces
        if verbose:
            logger.info('Looking for segments on reference '
                        f'{refName}{"+" if strd else "-"}...')
        recDict = SeqIO.index(fastafile, 'fasta')
        SNRpieces = findSNRpieces(
            covLen = covLen, seq = str(recDict[refName].seq),
            strd = strd, base = base, mism = mism,
            minSNRlen = minSNRlen)
        SNRpieces = flattenIntervals(SNRpieces)
        # Find overlaps between SNRpieces & tStarts
        SNRpieces = getOverlaps(flatStarts, SNRpieces)
    # Go over all the pieces to fetch the alignments & write them into a temp
    if verbose:
        logger.info('Identifying alignments to be removed on reference '
                    f'{refName}{"+" if strd else "-"}...')
    counter = 0
    tempFile = f'{out_bamfile}_{strd}_{refName}.bam'
    # Connect to the bam file and create a new one
    bam = pysam.AlignmentFile(bamfile, 'rb')
    bamTEMP = pysam.AlignmentFile(tempFile, 'wb', template = bam)
    for start, end in SNRpieces:
        for alignment in bam.fetch(
                contig = refName, start = start, stop = end):
            # Ensure the alignment is on the correct strd
            if alignment.is_reverse != strd:
                # Make sure that filtered alignments indeed map to
                #  the area identified, not just overlap it with
                #  their start-end range (in case of splicing)
                mapping = getOverlaps(alignment.get_blocks(),
                                      [(start, end)])
                if mapping != []:
                    bamTEMP.write(alignment)
                    counter += 1
    
    bam.close()
    bamTEMP.close()
    
    if verbose:
        logger.info(f'Identified {counter:,d} alignments to be removed on '
                    f'reference {refName}{"+" if strd else "-"} and saved to '
                    f'{tempFile}.')


def scanBAMfilter(covLen, minSNRlen, bamfile, fastafile, out_transBaselineData,
                  out_db, base = 'A', mism = 0, out_bamfile = None, tROC = {},
                  includeIntrons = False, cbFile = None, out_cbFile = None,
                  weightedCov = True, verbose = False, nThreads = None):
    """Slower version of BAMfilter() that, in order to save RAM, does not
    require loading pre-scanned SNRs; it obtains the intervals upstream of SNRs
    directly by scanning the fasta reference instead.
    Filters an indexed BAM file to remove non-canonical alignments that likely
    resulted from poly(dT) priming onto genomically encoded polyA single
    nucleotide repeats (SNRs), as opposed to mRNA polyA tails. Note that
    parallel processing will produce the following warning for each temp file,
    which can be safely ignored:
        "[E::idx_find_and_load] Could not retrieve index file for <file>"
    
    Parameters
    ----------
    covLen : (int)
        The assumed distance between the beginning of alignment coverage and
        the priming event.
    minSNRlen : (int, None)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    fastafile : (str)
        Location of the reference genome fasta file.
    out_transBaselineData : (str)
        Path to the saved baseline data, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    base : (str), optional
        DNA base constituting SNRs to be searched for: 'A', 'T', 'C', 'G', 'N'.
        The default is 'A'.
    mism : (int), optional
        The maximum number of mismatches in an SNR allowed. The default is 0.
    out_bamfile : (str, None), optional
        Location of the resulting filtered bamfile. If None, the name of the
        input bamfile name is used with '.filtered.bam' appended. The default
        is None.
    tROC : (dict), optional
        { endLen : ( TP, FN, TN, FP, eachTransStartByStrdRef ) }. The default
        is {}.
    includeIntrons : (bool), optional
        Indicates whether to consider intron coverage. The default is False.
    cbFile : (str, None), optional
        Location of the scumi cell barcode count file, if any. The default is
        None.
    out_cbFile : (str, None), optional
        Location of a new scumi cell barcode count file, if any. If None, the
        cbFile name is used with '.filtered.tsv' appended. The default is None.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Only included for the purposes of using
        cached baseline files; otherwise has no impact on the results of this
        function. The default is True.
    verbose : (bool), optional
        Indicates whether extra messages are logged. The default is False.
    nThreads : (int), optional
        Set the maximum number of threads to use besides the main thread.
        If None, serial processing is used. The default is None.
    """
    
    # Initiate the correct output file name
    if out_bamfile is None:
        out_bamfile = bamfile + '.filtered.bam'
    # If the BAM file already exists, do not overwrite
    if os.path.isfile(out_bamfile):
        logger.info(f'The filtered BAM file "{out_bamfile}" already exists.')
        return
        
    # Initialize the set of alignments to be removed
    toRemove = set()
    
    # Get transcript starts if N/A for this specific length
    if covLen not in tROC:
        BLdata = getBaselineData(out_transBaselineData, out_db, bamfile,
                                 includeIntrons, weightedCov = weightedCov)
        tROC[covLen] = getTransEndSensSpec(
            covLen, bamfile, BLdata, includeIntrons, getSensSpec = False,
            weightedCov = weightedCov)    
    # Extract the transcript starts for this coverage length
    eachTransStartByStrdRef = tROC[covLen][4]
    
    if nThreads is None:
        logger.info(f'Identifying alignments {covLen:,d} bp upstream of non-'
                    f'terminal SNR{minSNRlen}+ with {mism} mismatches to be '
                    'removed in 1 serial process...')
        # Connect to the fasta reference if needed
        if minSNRlen is not None:
            recDict = SeqIO.index(fastafile, 'fasta')
        
        # Connect to the bam file
        bam = pysam.AlignmentFile(bamfile, 'rb')
        
        # Go over each strand and ref separately to obtain overlaps with SNRs
        for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
            for refName, eachTransStart in eachTransStartByRef.items():
                # Flatten all the expressed transcript starts on this strd/ref
                flatStarts = []
                for tStart in eachTransStart:
                    flatStarts.extend(tStart)
                flatStarts = flattenIntervals(flatStarts)
                # minSNRlen = None means that all non-end coverage is removed,
                #  regardless of SNR pieces
                if minSNRlen is None:
                    SNRpieces = flatStarts
                else:
                    # Scan for and flatten all the relevant SNR pieces
                    if verbose:
                        logger.info('Looking for segments on reference '
                                    f'{refName}{"+" if strd else "-"}...')
                    SNRpieces = findSNRpieces(
                        covLen = covLen, seq = str(recDict[refName].seq),
                        strd = strd, base = base, mism = mism,
                        minSNRlen = minSNRlen)
                    SNRpieces = flattenIntervals(SNRpieces)
                    # Find overlaps between SNRpieces & tStarts
                    SNRpieces = getOverlaps(flatStarts, SNRpieces)
                # Fetch the alignments
                if verbose:
                    logger.info('Identifying alignments to be removed on refe'
                                f'rence {refName}{"+" if strd else "-"}...')
                for start, end in SNRpieces:
                    for alignment in bam.fetch(
                            contig = refName, start = start, stop = end):
                        # Ensure the alignment is on the correct strd
                        if alignment.is_reverse != strd:
                            # Make sure that filtered alignments indeed map to
                            #  the area identified, not just overlap it with
                            #  their start-end range (in case of splicing)
                            mapping = getOverlaps(alignment.get_blocks(),
                                                  [(start, end)])
                            if mapping != []:
                                toRemove.add(alignment)
        bam.close()
    
    else:                
        logger.info(f'Identifying alignments {covLen:,d} bp upstream of non-'
                    f'terminal SNR{minSNRlen}+ to be removed across {nThreads} '
                    'parallel processes and writing temporary files...')
        # Set the number of threads for NumExpr (used by pandas & numpy)
        os.environ['NUMEXPR_MAX_THREADS'] = str(nThreads)
        # Create a pool of processes with shared variables
        pool = multiprocessing.Pool(processes = nThreads)
        # Identify the alignments to be removed by strand / ref
        strdRefs = []
        for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
            for refName, eachTransStart in eachTransStartByRef.items():
                strdRefs.append((strd, refName))
                pool.apply_async(func = scanForAlignmentsToRemove,
                                 args = (covLen, refName, strd, minSNRlen,
                                         bamfile, out_bamfile, fastafile,
                                         eachTransStart, base, mism, verbose))
        # Close the pool
        pool.close()
        # Join the processes
        pool.join()        
        
        # Merge the temporary files into 1 in-memory set & remove the files
        logger.info(f'Merging {len(strdRefs)} temporary files into memory and '
                    'deleting...')
        for strand, ref in strdRefs:
            tempfile = f'{out_bamfile}_{strand}_{ref}.bam'
            bamTEMP = pysam.AlignmentFile(tempfile, 'rb')
            for alignment in bamTEMP.fetch(until_eof = True):
                toRemove.add(alignment)
            bamTEMP.close()
            os.remove(tempfile)
    
    toRemoveN = len(toRemove)
    # Filter the cbFile, if any
    if cbFile:
        cbFileFilter(toRemove, cbFile, out_cbFile, verbose)   
    
    # Create the bamfile and add the reads not in the toRemove set
    # Manually reopen the bam file to avoid multiple iterators issues
    bamIN = pysam.AlignmentFile(bamfile, 'rb')   
    # Get the total number of reads in the indexed bam file
    total = bamIN.mapped + bamIN.unmapped
    logger.info(f'Writing the filtered BAM file, to exclude {toRemoveN:,d} '
                f'alignments out of {total:,d} total...')
    # Write the new bam file
    bamOUT = pysam.AlignmentFile(out_bamfile, 'wb', template = bamIN)
    nAll = 0
    included = 0
    for alignment in bamIN.fetch(until_eof = True):
        nAll += 1
        if alignment not in toRemove:
            bamOUT.write(alignment)
            included += 1
            
        if not nAll % 10000000 and nAll and verbose:
            logger.info(f'Processed {nAll:,d} input alignments, of which '
                        f'{nAll-included:,d} were excluded...')
            
    # Close the files
    bamIN.close()
    bamOUT.close()
    logger.info(f'A filtered BAM file with {included:,d}/{nAll:,d} alignments '
                f'has been created. [Excluded {nAll-included:,d} alignments '
                f'{covLen:,d} bp upstream of non-terminal SNR{minSNRlen}+.]')
    

def fetchToKeepScanToRemove(covLen, refName, strd, minSNRlen, bamfile,
                            out_bamfile, fastafile, eachTransStart,
                            allTransEnds, base = 'A', mism = 0,
                            verbose = True):
    """Function equivalent to scanForAlignmentsToRemove() but additionally
    specifically identifies the alignments to keep.
    Helper function to run on each thread for a multithreaded scanning and
    alignment identification and writing into temporary files to be merged
    later.
    
    Parameters
    ----------
    covLen : (int)
        The assumed distance between the beginning of alignment coverage and
        the priming event.
    refName : (str)
        Name of the reference to scan.
    strd : (bool)
        Strand of the reference; +(True) / -(False)
    minSNRlen : (int, None)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    out_bamfile : (str)
        Location of the resulting filtered bamfile.
    fastafile : (str)
        Location of the reference genome fasta file.
    eachTransStart : (tuple)
        Tuple of trans start pieces such that
        (((start, end), (start, end), ...), ((start, end), (start, end), ...))
    allTransEnds: (list)
        A list of trans end pieces flattened across the given reference:
        [(start, end), (start, end), ...]
    base : (str), optional
        DNA base constituting SNRs to be searched for: 'A', 'T', 'C', 'G', 'N'.
        The default is 'A'.
    mism : (int), optional
        The maximum number of mismatches in an SNR allowed. The default is 0.
    verbose : (bool), optional
        Indicates whether extra messages are logged. The default is False.
    
    Returns
    -------
    None.
    """        
    
    # Flatten all the expressed transcript starts on this strd/ref
    flatStarts = []
    for tStart in eachTransStart:
        flatStarts.extend(tStart)
    flatStarts = flattenIntervals(flatStarts)
    # minSNRlen = None means that all non-end coverage is removed,
    #  regardless of SNR pieces
    if minSNRlen is None:
        SNRpieces = flatStarts
    else:
        # Scan for and flatten all the relevant SNR pieces
        if verbose:
            logger.info('Looking for segments on reference '
                        f'{refName}{"+" if strd else "-"}...')
        recDict = SeqIO.index(fastafile, 'fasta')
        SNRpieces = findSNRpieces(
            covLen = covLen, seq = str(recDict[refName].seq),
            strd = strd, base = base, mism = mism,
            minSNRlen = minSNRlen)
        SNRpieces = flattenIntervals(SNRpieces)
        # Find overlaps between SNRpieces & tStarts
        SNRpieces = getOverlaps(flatStarts, SNRpieces)
    # Go over all the pieces to fetch the alignments & write them into a temp
    if verbose:
        logger.info('Identifying alignments to be removed on reference '
                    f'{refName}{"+" if strd else "-"}...')
    
    tempFile = f'{out_bamfile}_{strd}_{refName}.bam'
    
    # Initialize the set of alignments to be removed
    toRemove = set()
    # Connect to the bam file and create a new one
    bam = pysam.AlignmentFile(bamfile, 'rb')

    for start, end in SNRpieces:
        for alignment in bam.fetch(
                contig = refName, start = start, stop = end):
            # Ensure the alignment is on the correct strd
            if alignment.is_reverse != strd:
                # Make sure that filtered alignments indeed map to
                #  the area identified, not just overlap it with
                #  their start-end range (in case of splicing)
                mapping = getOverlaps(alignment.get_blocks(),
                                      [(start, end)])
                if mapping != []:
                    toRemove.add(alignment)
    toRemoveN = len(toRemove)
    
    # Identify all terminal transcripts that at least overlap the transcript
    #  ends to protect them
    if verbose:
        logger.info('Identifying alignments to be protected on reference '
                    f'{refName}{"+" if strd else "-"}...')
    toProtect = set()
    for start, end in allTransEnds:
        for alignment in bam.fetch(
                contig = refName, start = start, stop = end):
            # Ensure the alignment is on the correct strd
            if alignment.is_reverse != strd:
                # Make sure that filtered alignments indeed map to the area
                #  identified, not just overlap it with their start-end range
                #  (in case of splicing)
                mapping = getOverlaps(alignment.get_blocks(),
                                      [(start, end)])
                if mapping != []:
                    toProtect.add(alignment)
    toProtectN = len(toProtect)
    
    bamTEMP = pysam.AlignmentFile(tempFile, 'wb', template = bam)
    bam.close()
    
    if verbose:
        logger.info('Protecting terminal alignments on reference '
                    f'{refName}{"+" if strd else "-"}...')
    toRemove -= toProtect
    
    removedN = len(toRemove)
    for alignment in toRemove:
        bamTEMP.write(alignment)
    bamTEMP.close()
    
    if verbose:
        logger.info(f'Alignments to be removed from {refName}'
                    f'{"+" if strd else "-"} were saved to {tempFile}:\n'
                    f'Identified to remove: {toRemoveN:,d}\n'
                    f'Identified to protect: {toProtectN:,d}\n'
                    f'Removed: {removedN:,d} (Protected: '
                    f'{toRemoveN-removedN:,d})')
    

def spBAMfilter(covLen, minSNRlen, bamfile, fastafile, out_transBaselineData,
                out_db, base = 'A', mism = 0, out_bamfile = None, tROC = {},
                includeIntrons = False, cbFile = None, out_cbFile = None,
                weightedCov = True, verbose = False, nThreads = None):
    """Version of scanBAMfilter() that, in addition to protecting (removing)
    the terminal end pieces of each expressed transcript, also specifically
    protects all transcripts that even just overlap them (not those that are
    just completely contained within them).
    Filters an indexed BAM file to remove non-canonical alignments that likely
    resulted from poly(dT) priming onto genomically encoded polyA single
    nucleotide repeats (SNRs), as opposed to mRNA polyA tails. Note that
    parallel processing will produce the following warning for each temp file,
    which can be safely ignored:
        "[E::idx_find_and_load] Could not retrieve index file for <file>"
    
    Parameters
    ----------
    covLen : (int)
        The assumed distance between the beginning of alignment coverage and
        the priming event.
    minSNRlen : (int, None)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    bamfile : (str)
        Location of the sorted and indexed bamfile to be filtered.
    fastafile : (str)
        Location of the reference genome fasta file.
    out_transBaselineData : (str)
        Path to the saved baseline data, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    base : (str), optional
        DNA base constituting SNRs to be searched for: 'A', 'T', 'C', 'G', 'N'.
        The default is 'A'.
    mism : (int), optional
        The maximum number of mismatches in an SNR allowed. The default is 0.
    out_bamfile : (str, None), optional
        Location of the resulting filtered bamfile. If None, the name of the
        input bamfile name is used with '.filtered.bam' appended. The default
        is None.
    tROC : (dict), optional
        { endLen : ( TP, FN, TN, FP, eachTransStartByStrdRef,
                    allTransEndsByStrdRef ) }.
        The default is {}.
    includeIntrons : (bool), optional
        Indicates whether to consider intron coverage. The default is False.
    cbFile : (str, None), optional
        Location of the scumi cell barcode count file, if any. The default is
        None.
    out_cbFile : (str, None), optional
        Location of a new scumi cell barcode count file, if any. If None, the
        cbFile name is used with '.filtered.tsv' appended. The default is None.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Only included for the purposes of using
        cached baseline files; otherwise has no impact on the results of this
        function. The default is True.
    verbose : (bool), optional
        Indicates whether extra messages are logged. The default is False.
    nThreads : (int), optional
        Set the maximum number of threads to use besides the main thread.
        If None, serial processing is used. The default is None.
    """
    
    # Initiate the correct output file name
    if out_bamfile is None:
        out_bamfile = bamfile + '.filtered.bam'
    # If the BAM file already exists, do not overwrite
    if os.path.isfile(out_bamfile):
        logger.info(f'The filtered BAM file "{out_bamfile}" already exists.')
        return
        
    # Initialize the set of alignments to be removed
    toRemove = set()
    
    # Get transcript starts if N/A for this specific length
    if covLen not in tROC:
        BLdata = getBaselineData(out_transBaselineData, out_db, bamfile,
                                 includeIntrons, weightedCov = weightedCov)
        tROC[covLen] = getTransEndSensSpec(
            covLen, bamfile, BLdata, includeIntrons, getSensSpec = False,
            weightedCov = weightedCov)    
    # Extract the transcript starts & ends for this coverage length
    eachTransStartByStrdRef = tROC[covLen][4]   # By transcript
    allTransEndsByStrdRef = tROC[covLen][5]     # Already flattened across ref
    
    if nThreads is None:
        logger.info(f'Identifying alignments {covLen:,d} bp upstream of non-'
                    f'terminal SNR{minSNRlen}+ with {mism} mismatches to be '
                    'removed in 1 serial process...')
        # Connect to the fasta reference if needed
        if minSNRlen is not None:
            recDict = SeqIO.index(fastafile, 'fasta')
        
        # Connect to the bam file
        bam = pysam.AlignmentFile(bamfile, 'rb')
        
        # Initiate the set of all alignments to be protected
        toProtect = set()
        
        # Go over each strand and ref separately to obtain overlaps with SNRs
        for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
            for refName, eachTransStart in eachTransStartByRef.items():
                # Flatten all the expressed transcript starts on this strd/ref
                flatStarts = []
                for tStart in eachTransStart:
                    flatStarts.extend(tStart)
                flatStarts = flattenIntervals(flatStarts)
                # minSNRlen = None means that all non-end coverage is removed,
                #  regardless of SNR pieces
                if minSNRlen is None:
                    SNRpieces = flatStarts
                else:
                    # Scan for and flatten all the relevant SNR pieces
                    if verbose:
                        logger.info('Looking for segments on reference '
                                    f'{refName}{"+" if strd else "-"}...')
                    SNRpieces = findSNRpieces(
                        covLen = covLen, seq = str(recDict[refName].seq),
                        strd = strd, base = base, mism = mism,
                        minSNRlen = minSNRlen)
                    SNRpieces = flattenIntervals(SNRpieces)
                    # Find overlaps between SNRpieces & tStarts
                    SNRpieces = getOverlaps(flatStarts, SNRpieces)
                # Fetch the alignments
                if verbose:
                    logger.info('Identifying alignments to be removed on refe'
                                f'rence {refName}{"+" if strd else "-"}...')
                for start, end in SNRpieces:
                    for alignment in bam.fetch(
                            contig = refName, start = start, stop = end):
                        # Ensure the alignment is on the correct strd
                        if alignment.is_reverse != strd:
                            # Make sure that filtered alignments indeed map to
                            #  the area identified, not just overlap it with
                            #  their start-end range (in case of splicing)
                            mapping = getOverlaps(alignment.get_blocks(),
                                                  [(start, end)])
                            if mapping != []:
                                toRemove.add(alignment)
                # Identify all terminal transcripts that at least overlap the
                #  transcript ends to protect them
                for start, end in allTransEndsByStrdRef[strd][refName]:
                    for alignment in bam.fetch(
                            contig = refName, start = start, stop = end):
                        # Ensure the alignment is on the correct strd
                        if alignment.is_reverse != strd:
                            # Make sure that filtered alignments indeed map to
                            #  the area identified, not just overlap it with
                            #  their start-end range (in case of splicing)
                            mapping = getOverlaps(alignment.get_blocks(),
                                                  [(start, end)])
                            if mapping != []:
                                toProtect.add(alignment)
        bam.close()
        toRemove -= toProtect
    
    else:                
        logger.info(f'Identifying alignments {covLen:,d} bp upstream of non-'
                    f'terminal SNR{minSNRlen}+ to be removed across {nThreads} '
                    'parallel processes and writing temporary files...')
        # Set the number of threads for NumExpr (used by pandas & numpy)
        os.environ['NUMEXPR_MAX_THREADS'] = str(nThreads)
        # Create a pool of processes with shared variables
        pool = multiprocessing.Pool(processes = nThreads)
        # Identify the alignments to be removed by strand / ref
        strdRefs = []
        for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
            for refName, eachTransStart in eachTransStartByRef.items():
                strdRefs.append((strd, refName))
                pool.apply_async(func = fetchToKeepScanToRemove,
                                 args = (covLen, refName, strd, minSNRlen,
                                         bamfile, out_bamfile, fastafile,
                                         eachTransStart,
                                         allTransEndsByStrdRef[strd][refName],
                                         base, mism, verbose))
        # Close the pool
        pool.close()
        # Join the processes
        pool.join()        
        
        # Merge the temporary files into 1 in-memory set & remove the files
        logger.info(f'Merging {len(strdRefs)} temporary files into memory and '
                    'deleting...')
        for strand, ref in strdRefs:
            tempfile = f'{out_bamfile}_{strand}_{ref}.bam'
            bamTEMP = pysam.AlignmentFile(tempfile, 'rb')
            toRemove.update({a for a in bamTEMP.fetch(until_eof = True)})
            bamTEMP.close()
            os.remove(tempfile)
    
    toRemoveN = len(toRemove)
    # Filter the cbFile, if any
    if cbFile:
        cbFileFilter(toRemove, cbFile, out_cbFile, verbose)   
    
    # Create the bamfile and add the reads not in the toRemove set
    # Manually reopen the bam file to avoid multiple iterators issues
    bamIN = pysam.AlignmentFile(bamfile, 'rb')   
    # Get the total number of reads in the indexed bam file
    total = bamIN.mapped + bamIN.unmapped
    logger.info(f'Writing the filtered BAM file, to exclude {toRemoveN:,d} '
                f'alignments out of {total:,d} total...')
    # Write the new bam file
    bamOUT = pysam.AlignmentFile(out_bamfile, 'wb', template = bamIN)
    nAll = 0
    included = 0
    for alignment in bamIN.fetch(until_eof = True):
        nAll += 1
        if alignment not in toRemove:
            bamOUT.write(alignment)
            included += 1
            
        if not nAll % 10000000 and nAll and verbose:
            logger.info(f'Processed {nAll:,d} input alignments, of which '
                        f'{nAll-included:,d} were excluded...')
            
    # Close the files
    bamIN.close()
    bamOUT.close()
    logger.info(f'A filtered BAM file with {included:,d}/{nAll:,d} alignments '
                f'has been created. [Excluded {nAll-included:,d} alignments '
                f'{covLen:,d} bp upstream of non-terminal SNR{minSNRlen}+.]')

