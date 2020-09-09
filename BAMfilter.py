#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 09:11:28 2020

@author: marek
"""
import os
import pysam
from datetime import datetime

from RNAseqAnalysis import getBaselineData, getTransEndSensSpec, \
    sortSNRsByLenStrdRef
from SNRanalysis import flattenIntervals, getOverlaps


def BAMfilter(lenToSNRs, covLen, minSNRlen, bamfile, out_transBaselineData, 
              out_db, out_bamfile, tROC = {}, includeIntrons = False):
    
    
    # If the BAM file already exists, do not overwrite
    if os.path.isfile(out_bamfile):
        print('The filtered BAM file already exists.')
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
        for refName, eachTransStart in eachTransStartByRef.itmes:
            print('{} - Identifying the reads to be removed on reference {}' \
                  '{}...'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                 refName, '+' if strd else '-'))
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
                            start = SNR.start - covLen if strd else SNR.end
                            end = SNR.start if strd else SNR.end + covLen
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
    
    print('{} - Filtering out reads from the BAM file...'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    # Go over all aligned reads in the bam file, removing the relevant reads
    toKeep = [read for read in bam.fetch() if read not in toRemove]
    print('{} - Sorting the BAM file...'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    # Reads must be sorted for indexing
    toKeep.sort(key = lambda x:
                (x.reference_id, x.reference_start, x.reference_length))
    
    print('{} - Writing the filtered BAM file...'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    # Create the bamfile to add the reads
    filtBAM = pysam.AlignmentFile(out_bamfile, 'wb', template = bam)
    # Add the reads in order
    for read in toKeep:
        filtBAM.write(read)
    # Close the files
    filtBAM.close()
    bam.close()
    
    print('{} - Indexing the filtered BAM file...'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    pysam.index(filtBAM)
    
    print('{} - A filtered bam file has been created and indexed.'.format(
        datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
        
    
        