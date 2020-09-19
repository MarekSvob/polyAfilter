#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:11:28 2020

@author: marek
"""
import os
import pysam
import logging

from RNAseqAnalysis import getBaselineData, getTransEndSensSpec, \
    sortSNRsByLenStrdRef
from SNRanalysis import flattenIntervals, getOverlaps


logging.basicConfig(level = logging.INFO,
                    format = '%(asctime)s - %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def BAMfilter(lenToSNRs, covLen, minSNRlen, bamfile, out_transBaselineData, 
              out_db, out_bamfile = None, tROC = {}, includeIntrons = False,
              cbFile = None, out_cbFile = None):
    """Function that filters a BAM file to remove non-canonical reads that
    likely resulted from poly(dT) priming onto genomically encoded polyA single
    nucleotide repeats (SNRs), as opposed to mRNA polyA tails.

    Parameters
    ----------
    lenToSNRs : (dict)
        { SNR length : [ SNR ] }
    covLen : (int)
        The assumed distance between the priming event and read coverage.
    minSNRlen : (int)
        The shortest SNR length to consider. If None, remove all non-end
        coverage.
    bamfile : (str)
        Location of the bamfile to be filtered.
    out_transBaselineData : (str)
        Path to the saved baseline data, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    out_bamfile : (str)
        Location of the resulting filtered bamfile. If None, the name of the
        input bamfile name is used with '.filtered.bam' appended. The default
        is None.
    tROC : (dict), optional
        { endLen : ( Sensitivity, Specificity, Accuracy ) }. The default is {}.
    includeIntrons : (bool), optional
        Indicates whether to consider intron coverage. The default is False.
    cbFile : (str), optional
        Location of the scumi cell barcode count file, if any. The default is
        None.
    out_cbFile : (str), optional
        Location of a new scumi cell barcode count file, if any. If None, the
        cbFile name is used with '.filtered.tsv' appended. The default is None.

    Returns
    -------
    None.
    """
    
    # Initiate the correct output file name
    if out_bamfile is None:
        out_bamfile = bamfile + '.filtered.bam'
    # If the BAM file already exists, do not overwrite
    if os.path.isfile(out_bamfile):
        logger.info('The filtered BAM file already exists.')
        return
    # Get transcript starts if N/A for this specific length
    if covLen not in tROC:
        BLdata = getBaselineData(out_transBaselineData, out_db, bamfile,
                                 includeIntrons)
        tROC[covLen] = getTransEndSensSpec(covLen, bamfile, BLdata,
                                           includeIntrons, getSensSpec = False)    
    # Extract the transcript starts for this coverage length
    eachTransStartByStrdRef = tROC[covLen][2]
    # If relevant, pre-sort SNRs by len, strd & ref    
    if minSNRlen is not None:
        SNRsByLenStrdRef = sortSNRsByLenStrdRef(lenToSNRs)
    # Initialize the set of reads to remove
    toRemove = set()
    # Connect to the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # Go over each strand and ref separately
    for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
        for refName, eachTransStart in eachTransStartByRef.items():
            logger.info('Identifying the reads to be removed on reference '
                        f'{refName}{"+" if strd else "-"}...')
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
            for start, end in SNRpieceOverlaps:
                for read in bam.fetch(contig = refName, start = start,
                                      stop = end):
                    # Make sure the read is on the correct strand before adding
                    if read.is_reverse != strd:
                        # Make sure that reads filtered do indeed map onto the
                        #  area identified, not just overlap it with thir
                        #  start-end range (in case of splicing)
                        mapping = getOverlaps(read.get_blocks(), [(start, end)])
                        if mapping != []:
                            toRemove.add(read)
    bam.close()
    toRemoveN = len(toRemove)
    # Filter the cbFile, if any
    if cbFile and toRemoveN:
        cbFileFilter(toRemove, cbFile, out_cbFile)
        
    logger.info(f'Writing the filtered BAM file, excluding {toRemoveN:,d} '
                'reads...')
    # Create the bamfile and add the reads not in the toRemove set
    nReads = 0
    
    bamIN = pysam.AlignmentFile(bamfile, 'rb')
    bamOUT = pysam.AlignmentFile(out_bamfile, 'wb', template = bam)
    for read in bamIN.fetch(until_eof = True):
        if read not in toRemove:
            nReads += 1
            bamOUT.write(read)
    # Close the files
    bamIN.close()
    bamOUT.close()
    logger.info(f'A filtered BAM file with {nReads} reads has been created.')
        
    
def cbFileFilter(toRemove, cbFile, out_cbFile):
    """Function that subtracts the appropriate numbers from the CB counts,
    including the T frequencies in UMIs. Note that where related, parts of this
    function's code were adapted from the scumi module
    [https://bitbucket.org/jerry00/scumi-dev/src/master/].

    Parameters
    ----------
    toRemove : (set)
        The set of bamfile reads to remove.
    cbFile : (str), optional
        Location of the scumi cell barcode count file.
    out_cbFile : (str), optional
        Location of a new cell barcode count file. If None, the cbFile name
        is used with '.filtered.tsv' appended.

    Returns
    -------
    None.
    """
    
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
    logger.info(f'{nUMIs:,d} reads across {nCells:,d} cell barcodes '
                'removed from the CB file.')