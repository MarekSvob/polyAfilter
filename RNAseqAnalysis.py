#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 17:00:37 2020

@author: marek
"""
import os
import pysam
import collections
import numpy as np
from SNRanalysis import getStrandedFeats, getStrandedFeatsByGene
from SNRdetection import savePKL, loadPKL

expCovConc = 0
expCovDisc = 0


def filterReads(
        filt_bamfile,
        out_db,
        out_strandedFeats,
        bamfile,
        featType = 'exon',
        concordant = True
        ):
    """Function to create a bam file from a subset of reads that are aligned
    to a specific feature type on the same or opposite strand.

    Parameters
    ----------
    filt_bamfile : (str)
        The path to the filtered bamfile.
    out_db : (str)
        Path to the file with the database.
    out_strandedFeats : (str)
        Path to the file with stranded feats.
    bamfile : (str)
        The path to the original bamfile.
    featType : TYPE, optional
        Get only reads overlapping this type of feature. The default is 'exon'.
    concordant : (bool), optional
        Indicates whether the mapping reads should be on the same (True) or
        opposite (False) strand wrt/ the feature. The default is True.

    Returns
    -------
    None.
    """
    
    # If the file already exists, announce and do not do anything
    if os.path.isfile(filt_bamfile):
        print('The filtered file already exists.')
        if not os.path.isfile('{}.bai'.format(filt_bamfile)):
            print('!!! Before continuing, index it in terminal using command '\
                  '"samtools index {}"'.format(filt_bamfile))
        return    
    
    # This will usually just load the file
    flatStrandedFeats = getStrandedFeats(
        out_strandedFeats,
        out_db,
        (featType, )
        )
    
    # Connect to the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # Initiate the varuables needed
    # This needs to start off as a set to avoid double-adding reads aligned to
    #  two features of the same type
    filt_reads = set()
    current_ref = ''
    strds = ('+', '-')
    
    # Get only reads on the same/other strand (depends if concordant or not)
    #  as the feature in question
    for strd in strds:
        # Cycle over the feats under the appropriate key
        for feat in sorted(flatStrandedFeats[strd][featType]):
            if current_ref != feat[0]:
                current_ref = feat[0]
                print(
                    'Processing {} strand of reference {}'.format(
                        strd,
                        current_ref
                        )
                    )
            # Gffutils features have 1-based [x,y] closed intervals,
            #  while pysam fetch has 0-based [x,y) half-open intervals
            for read in bam.fetch(
                    contig = feat[0],
                    start = feat[1] - 1,
                    end = feat[2]
                    ):
                # Make sure the read is on the correct strand before adding:
                if concordant == (read.is_reverse == (strd == '-')):
                    filt_reads.add(read)
    
    # Transform the set into a list for sorting
    filt_reads_list = list(filt_reads)
    # Reads must be sorted for indexing
    filt_reads_list.sort(
        key = lambda x: (x.reference_id, x.reference_start, x.reference_length)
        )
    
    # Create the bamfile to add the reads
    filt_bam = pysam.AlignmentFile(filt_bamfile, 'wb', template = bam)
    # Add the reads in order
    for read in filt_reads_list:
        filt_bam.write(read)
    # Close the files
    filt_bam.close()
    bam.close()
    
    print('A filtered bam file has been created.')
    print('!!! Before continuing, index it in terminal using command '\
          '"samtools index {}"'.format(filt_bamfile))


def getExpectedCoverage(
        out_db,
        out_strandedFeats,
        bamfile,
        concordant = True,
        baselineFeat = 'transcript'
        ):
    """Get the per-base expected RNA-seq coverage across the given featuretype.
    
    Parameters
    ----------
    out_db : (str)
        Path to the file with the database.
    out_strandedFeats : (str)
        Path to the file with stranded feats.
    bamfile : (str)
        Path to the (pre-filtered) bamfile.
    concordant : (bool)
        Indicator of which saved value to retrieve.
    baselineFeat : (str), optional
        The featureType whose total length is used for normalization. The
        default is 'transcript'.

    Returns
    -------
    expected : (float)
        The expected per-base coverage from the RNA-seq data.
    """
    
    global expCovConc, expCovDisc
    
    # Try to retrieve the value if already calculated in this session
    if concordant and expCovConc != 0:
        expected = expCovConc
        return expected
    if not concordant and expCovDisc != 0:
        expected = expCovDisc
        return expected
    
    # In most cases, this should just load a file
    flatStrandedFeats = getStrandedFeats(
        out_strandedFeats,
        out_db,
        (baselineFeat, )
        )
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # First, get the total (absolute) coverage by the reads on both strands in
    #  the bam file in question
    totalCoverage = 0
    for ref in bam.references:
        print('Getting total coverage for reference {}'.format(ref))
        totalCoverage += np.sum(
            bam.count_coverage(ref, quality_threshold = 0),
            dtype = 'int'
            )
    bam.close()
    
    # Get the total length of the feat that the SNRs come from
    totalLength = 0
    strds = ('+', '-')
    for strd in strds:
        print('Scanning {}{}s for total length...'.format(strd, baselineFeat))
        for feat in flatStrandedFeats[strd][baselineFeat]:
            # Note that pysam is 0-based; my vars are 1-based (like the gtf)
            totalLength += feat[2] - feat[1] + 1
    
    # Normalize the coverage by #SNRs, total coverage, and total length
    expected = totalCoverage / totalLength
    
    # Save the value for later use in this session
    if concordant:
        expCovConc = expected
    else:
        expCovDisc = expected    
    
    return expected


def getCovPerSNR(
        out_SNRCovGeneLen,
        out_SNROutl,
        lenToSNRs,
        out_db,
        out_strandedFeats,
        exonic_bamfile,
        SNRlengths = range(200),
        window = 2000,
        concordant = True,
        SNRfeat = 'transcript',
        log10Filter = 3
        ):
    """Obtain RNA-seq coverage in the vicinity of SNRs of given length.
    
    Parameters
    ----------
    out_SNRCovLen : (str)
        Path to the file with per-SNR coverage by length.
    out_SNRCovGene : (str)
        Path to the file with per-SNR coverage by gene. Note that only non-0
        peaks are added.
    out_SNROutl : (str)
        Path to the file with per-SNR coverage by length for each outlier.
    lengthToSNRs : (dict)
        { length : [SNR] }
    out_db : (str)
        Path to the file with the database.
    out_strandedFeats : (str)
        Path to the file with stranded feats.
    exonic_bamfile : (str)
        Path to the (pre-filtered) bamfile.
    SNRlengths : (iter), optional
        An iterable (list, tuple, etc.) containing the SNR lengths of interest.
        The default is itertools.count().
    window : (int), optional
        Total size of the window around the SNR covered. The default is 2000.
    concordant : (bool), optional
        Switch between concordant & discordant coverage. The default is True.
    SNRfeat : (str), optional
        The feature by which SNRs are selected. The default is 'transcript'.
    log10Filter : (int), optional
        If coverage at any base of an SNR exceeds this power of 10 multiple of
        the expected coverage, the SNR is discarded as an outlier. The default
        is 3.

    Returns
    -------
    covByGeneLen : (dict)
        { gene : { length : ( SNRcount, SNRzeroCount, np.array ) } }
    outlierPeaks : (dict)
        { SNRlength : [ ( SNR, np.array ) ] }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_SNRCovGeneLen) and os.path.isfile(out_SNROutl):
            covByGeneLen = loadPKL(out_SNRCovGeneLen)
            outlierPeaks = loadPKL(out_SNROutl)
            return covByGeneLen, outlierPeaks

    # Get the expected coverage
    expCov = getExpectedCoverage(
        out_db,
        out_strandedFeats,
        exonic_bamfile,
        concordant = concordant,
        SNRfeat = SNRfeat
        )
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(exonic_bamfile, 'rb')
    
    # Initialize the dict to collect all the data
    covByGeneLen = collections.defaultdict(collections.defaultdict(tuple))
    outlierPeaks = collections.defaultdict(list)
    totalCount = 0
    filteredOut = 0
    zeroCovs = 0
    # Go over the SNRs by length, starting with the longest
    for length in sorted(lenToSNRs.keys(), reverse = True):
        if length in SNRlengths:
            # For each length, initialize the count & array for window coverage
            print(
                'Going over {} SNR{}s...'.format(
                    'concordant' if concordant else 'discordant',
                    length
                    )
                )
            # Note that SNRs are 0-based
            for SNR in lenToSNRs[length]:
                # Only consider SNRs in transcripts, as those are the
                #  sources of internal (concordant or discordant) priming
                if (concordant and SNRfeat in SNR.concFeats) \
                    or (not concordant and SNRfeat in SNR.discFeats):
                    # Get the mid point of the SNR wrt/ the reference
                    mid = (SNR.end + SNR.start) // 2
                    start = int(mid - window / 2)
                    stop = int(mid + window / 2)
                    # Include correctins for the start & end if the window
                    #  falls out of the reference size range
                    if start < 0:
                        corrStart = 0
                    else:
                        corrStart = start
                    refLen = bam.get_reference_length(SNR.record)
                    if stop > refLen:
                        corrStop = refLen
                    else:
                        corrStop = stop
                    # Get the coverage summed over A/T/C/G; count only
                    #  reads on the same (conc) or opposite (disc) strand
                    refCov = np.sum(
                        bam.count_coverage(
                            contig = SNR.record,
                            start = corrStart,
                            stop = corrStop,
                            quality_threshold = 0,
                            read_callback = lambda r: \
                                r.is_reverse != SNR.strand
                            ),
                        axis = 0
                        )
                    totalCount += 1
                    # Extract the appropriate genes
                    if concordant:
                        genesOfInterest = SNR.concGenes
                    else:
                        genesOfInterest = SNR.discGenes
                    # If this is zero coverage, only add the counts
                    if sum(refCov) == 0:
                        zeroCovs += 1
                        # A custom defaultdict workaround
                        for gene in genesOfInterest:
                            if covByGeneLen[gene][length] == ():
                                covByGeneLen[gene][length] = (
                                    1,
                                    1,
                                    np.zeros((window), 'L')
                                    )
                            else:
                                covByGeneLen[gene][length][0] += 1
                                covByGeneLen[gene][length][1] += 1
                    else:                    
                        # If the window was out of ref range, fill in the rest
                        if corrStart != start:
                            refCov = np.append(
                                np.zeros((corrStart - start), 'L'),
                                refCov
                                )
                        if corrStop != stop:
                            refCov = np.append(
                                refCov,
                                np.zeros((stop - corrStop), 'L')
                                )
                        # Normalize and if needed, flip the coverage to be in
                        #  the 5'->3' orientation wrt/ the transcript feature
                        if SNR.strand == concordant:
                            refCov = refCov / expCov
                        else:
                            refCov = refCov[::-1] / expCov
                        # If coverage at any base around this SNR exceeds the
                        #  threshold multiple of the expected coverage, put
                        #  the SNR's coverage aside into outliers
                        if max(refCov) / expCov > 10**log10Filter:
                            filteredOut += 1
                            outlierPeaks[length].append((SNR, refCov))
                        # Otherwise add the coverage by gene & length
                        else:
                            for gene in genesOfInterest:
                                if covByGeneLen[gene][length] == ():
                                    covByGeneLen[gene][length] = (1, 0, refCov)
                                else:
                                    covByGeneLen[gene][length][0] += 1
                                    covByGeneLen[gene][length][2] += refCov
    bam.close()
    savePKL(out_SNRCovGeneLen, covByGeneLen)
    savePKL(out_SNROutl, outlierPeaks)
    print(
        'Filtered out {:.2%} outliers. {:.2%} SNRs had nfo coverage.'.format(
            filteredOut / totalCount,
            zeroCovs / totalCount
            )
        )
    
    return covByGeneLen, outlierPeaks


def getCovPerTran(out_TranCov, concordant = True):
    
    
    # If the file already exists, simply load
    if os.path.isfile(out_SNRCovGeneLen) and os.path.isfile(out_SNROutl):
            covByGene = loadPKL(out_SNRCovGeneLen)
            return covByGene
    
    # Get the expected coverage
    expCov = getExpectedCoverage(
        out_db,
        out_strandedFeats,
        exonic_bamfile,
        concordant = concordant,
        SNRfeat = 'transcript'
        )
    
    # If done before, this will only return the feats
    featsByGene = getStrandedFeatsByGene(
        out_featsByGene,
        out_db,
        featType = 'transcript'
        )
    
    return