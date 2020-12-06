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
    # Import the modules that are only necessary if a cbFile was provided
    import numpy as np
    import pandas as pd
    import regex as re
    from collections import defaultdict
    from functools import partial
    from scumi import scumi as scu
    
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
    # Sort by cb_count & save
    diffDF.sort_values(by = 'cb_count', ascending = False, inplace = True)
    if out_cbFile is None:
        out_cbFile = cbFile + '.filtered.tsv'
    diffDF.to_csv(out_cbFile, sep = '\t')
    
    # Report the results
    nCells = rmDF.shape[0]
    nUMIs = sum(rmDF['cb_count'])
    if verbose:
        logger.info(f'{nUMIs:,d} reads across {nCells:,d} cell barcodes '
                    'removed from the CB file.')


def getAlignmentsToRemove(strd, refName, eachTransStart):
    """Helper function to run on each thread for a multithreaded alignment
    identification and write as temporary files to be merged later.
    
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
    https://stackoverflow.com/questions/25825995/python-multiprocessing-only-one-process-is-running
    
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
    https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
    """
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)


class MyPool(multiprocessing.pool.Pool):
    """Class for embedded parallel processes adapted from
    https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
    """
    Process = AngelProcess
    

def parallel_wrapper(covLen, minSNRlen, bamfile, out_SNRsByLenStrdRef,
                     out_transBaselineData, out_db, out_bamfile, tROC,
                     includeIntrons, sortedSNRs, verbose, nThreads):
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
                                 includeIntrons)
        tROC[covLen] = getTransEndSensSpec(covLen, bamfile, BLdata,
                                           includeIntrons, getSensSpec = False)    
    # Extract the transcript starts for this coverage length
    eachTransStartByStrdRef = tROC[covLen][4]
    # Load and, if relevant, pre-sort SNRs by len, strd & ref
    SNRsByLenStrdRef = loadPKL(out_SNRsByLenStrdRef)
    if minSNRlen is not None and not sortedSNRs:
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
              sortedSNRs = True, verbose = False, nThreads = None):
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
    verbose : (bool), optional
        Indicates whether extra messages are logged. The default is False.
    nThreads : (int), optional
        Set the maximum number of threads to use. If None, serial processing is
        used but not enforced (underlying processes may use multiple threads).
        The default is None.

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
                                     includeIntrons)
            tROC[covLen] = getTransEndSensSpec(covLen, bamfile, BLdata,
                                               includeIntrons, getSensSpec = False)    
        # Extract the transcript starts for this coverage length
        eachTransStartByStrdRef = tROC[covLen][4]
        # Load and, if relevant, pre-sort SNRs by len, strd & ref
        SNRsByLenStrdRef = loadPKL(out_SNRsByLenStrdRef)
        if minSNRlen is not None and not sortedSNRs:
            SNRsByLenStrdRef = sortSNRsByLenStrdRef(SNRsByLenStrdRef)
        
        logger.info(f'Identifying alignments {covLen:,d} bp upstream of non-'
                    f'terminal SNR{minSNRlen}+ to be removed in 1 serial '
                    'process...')
        # Connect to the bam file
        bam = pysam.AlignmentFile(bamfile, 'rb')
        
        # Go over each strand and ref separately
        for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
            for refName, eachTransStart in eachTransStartByRef.items():
                if verbose:
                    logger.info('Identifying alignments to be removed on refe'
                                f'rence {refName}{"+" if strd else "-"}...')
                refLen = bam.get_reference_length(refName)
                # Flatten all the expressed transcript starts on this strd/ref
                flatStarts = []
                for tStart in eachTransStart:
                    flatStarts.extend(tStart)
                flatStarts = flattenIntervals(flatStarts)
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
                                SNRpieces.append((start, end))
                    SNRpieces = flattenIntervals(SNRpieces)
                    # Find overlaps between SNRpieces & tStarts
                    SNRpieceOverlaps = getOverlaps(flatStarts, SNRpieces)
                
                # Remove SNRsByLenStrdRef to free up RAM for the alingments
                del SNRsByLenStrdRef
                gc.collect()
                
                # Go over all the pieces to fetch the alignments
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
        # Start a single parallel process to get rid of 
        pool = MyPool(processes = 1)
        result = pool.apply_async(
            func = parallel_wrapper,
            args = (covLen, minSNRlen, bamfile, out_SNRsByLenStrdRef,
                    out_transBaselineData, out_db, out_bamfile, tROC,
                    includeIntrons, sortedSNRs, verbose, nThreads))
        # Close the pool
        pool.close()
        # Join the processes
        pool.join()
        # Obtain the list of strd-ref pairs processed
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
    if cbFile and toRemoveN:
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
            
        if verbose and nAll and not nAll % 10000000:
            logger.info(f'Processed {nAll:,d} input alignments, of which '
                        f'{nAll-included:,d} were excluded...')
            
    # Close the files
    bamIN.close()
    bamOUT.close()
    logger.info(f'A filtered BAM file with {included:,d}/{nAll:,d} alignments '
                f'has been created. [Excluded {nAll-included:,d} alignments '
                f'{covLen:,d} bp upstream of non-terminal SNR{minSNRlen}+.]')
    