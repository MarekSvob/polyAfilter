#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 17:00:37 2020

@author: marek
"""
import os
import pysam
import gffutils
import collections
import numpy as np
from SNRanalysis import getStrandedFeats, getGeneLenToSNRs
from SNRdetection import savePKL, loadPKL
from IPython.display import clear_output

expCovByConcFeat = collections.defaultdict(lambda:collections.defaultdict(int))

expectedExCov = 0
covTransByStrdRef = collections.defaultdict(
    lambda: collections.defaultdict(list)
    )
covExonsByStrdRef = collections.defaultdict(
    lambda: collections.defaultdict(list)
    )

class NoncGene:
    """An object that stores information about a gene identified to be
    non-canonically covered.
    """    
    
    def __init__(self, geneF, cov, trans, SNRs, endings):
        self.geneFeat = geneF   # gffutils.Feature, gff db gene in question
        self.cov = cov          # np.array, the bp-wise coverage of this gene
        self.trans = trans      # list, [ Transcript ] of overlapping Ts
        self.SNRs = SNRs        # dict, { SNRlength : ( float, [ SNR ] ) }
        self.endings = endings  # list, [ (start, end) ] of exon-wise endings
        
    def __str__(self):
        # Displays the base, number of mismatches, and location
        return '{} NoncGene @{}-{}({})'.format(
            self.geneFeat.id,
            self.geneFeat.start,
            self.geneFeat.end,
            self.geneFeat.strand
            )

class Transcript:
    """An object that stores information about a transcript.
    """    
    
    def __init__(self, geneID, start, end, exons):
        self.geneID = geneID    # str, geneID of the gene the T belongs to
        self.start = start      # int, start position on the chromosome
        self.end = end          # int, end position on the chromosome
        self.exons = exons      # list, [ (start, end) ] of child Es; 0-based
        
    def __str__(self):
        # Displays the base, number of mismatches, and location
        return '{} Transcript @{}-{}({})'.format(
            self.geneID,
            self.start,
            self.end,
            self.strand
            )


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
    
    # Get only reads on the same/other strand (depends if concordant or not)
    #  as the feature in question
    for strd in (True, False):
        # Cycle over the feats under the appropriate key
        for feat in sorted(flatStrandedFeats[strd][featType]):
            if current_ref != feat[0]:
                current_ref = feat[0]
                print(
                    'Processing {} strand of reference {}'.format(
                        '+' if strd else '-',
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
                if concordant == (read.is_reverse != strd):
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
    
    global expCovByConcFeat
    
    # Try to retrieve the value if already calculated in this session
    if expCovByConcFeat[concordant][baselineFeat] != 0:
        expected = expCovByConcFeat[concordant][baselineFeat]
        return expected
    
    # In most cases, this should just load a file
    flatStrandedFeats = getStrandedFeats(
        out_strandedFeats,
        out_db,
        (baselineFeat, )
        )
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')

    totalCoverage = 0
    totalLength = 0

    # Get the total length of the feats and the total coverage over them
    for strd in (True, False):
        print(
            'Scanning {}{}s for total coverage and length...'.format(
                '+' if strd else '-',
                baselineFeat
                )
            )
        for feat in flatStrandedFeats[strd][baselineFeat]:
            # Note that pysam is 0-based; my vars are 1-based (like the gtf)
            totalLength += feat[2] - feat[1] + 1
            totalCoverage += np.sum(
                bam.count_coverage(
                    contig = feat[0],
                    start = feat[1] - 1,      # 1-based => 0-based
                    stop = feat[2],
                    quality_threshold = 0,
                    read_callback = lambda r: r.is_reverse != strd
                    ),
                dtype = 'int'
                )
    bam.close()
    
    # Normalize the coverage by #SNRs, total coverage, and total length
    expected = totalCoverage / totalLength
    
    # Save the value for later use in this session
    expCovByConcFeat[concordant][baselineFeat] = expected
    
    return expected


def getCovPerSNR(
        out_SNRCovLen,
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
    covByLen : (dict)
        { length : ( SNRcount, SNRzeroCount, np.array ) }
    outlierPeaks : (dict)
        { SNRlength : [ ( SNR, np.array ) ] }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_SNRCovLen) and os.path.isfile(out_SNROutl):
            covByLen = loadPKL(out_SNRCovLen)
            outlierPeaks = loadPKL(out_SNROutl)
            return covByLen, outlierPeaks

    # Get the expected coverage
    expCov = getExpectedCoverage(
        out_db,
        out_strandedFeats,
        exonic_bamfile,
        concordant = concordant,
        baselineFeat = SNRfeat
        )
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(exonic_bamfile, 'rb')
    
    # Initialize the dict to collect all the data
    covByLen = {}
    outlierPeaks = collections.defaultdict(list)
    totalCount = 0
    filteredOut = 0
    zeroCovs = 0
    # Go over the SNRs by length, starting with the longest
    for length in sorted(lenToSNRs.keys(), reverse = True):
        if length in SNRlengths:
            # For each length, initialize the count & array for window coverage
            print(
                'Going over {} coverage of SNR{}s...'.format(
                    'concordant' if concordant else 'discordant',
                    length
                    )
                )
            SNRcount = 0
            zeros = 0
            coverage = np.zeros((window), 'L')
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
                    # If this is zero coverage, only add the counts
                    if sum(refCov) == 0:
                        zeroCovs += 1
                        zeros += 1
                        SNRcount += 1
                    # Otherwise add the coverage found
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
                            refCov = refCov
                        else:
                            refCov = refCov[::-1]
                        # If coverage at any base around this SNR exceeds the
                        #  threshold multiple of the expected coverage, put
                        #  the SNR's coverage aside into outliers
                        if max(refCov) / expCov > 10**log10Filter:
                            filteredOut += 1
                            outlierPeaks[length].append((SNR, refCov / expCov))
                        # Otherwise add the coverage by length
                        else:
                            coverage += refCov
                            SNRcount += 1
            # At the end of each length normalize & add everything to dict
            covByLen[length] = (SNRcount, zeros, coverage / expCov)
    bam.close()
    savePKL(out_SNRCovLen, covByLen)
    savePKL(out_SNROutl, outlierPeaks)
    print(
        'Filtered out {:.2%} outliers and {:.2%} SNRs had no coverage.'.format(
            filteredOut / totalCount,
            zeroCovs / totalCount
            )
        )
    
    return covByLen, outlierPeaks


def getCovPerTran(
        out_TranCov,
        out_db,
        out_strandedFeats,
        exonic_bamfile,
        window = 2000
        ):
    """Function that aggregates all per-transcript coverage.

    Parameters
    ----------
    out_TranCov : (str)
        Path to the output of this function.
    out_db : (str)
        Path to the file with the database.
    out_strandedFeats : (str)
        Path to the file with stranded feats.
    exonic_bamfile : (str)
        Path to the (pre-filtered) bamfile.
    window : (int), optional
        Total size of the window around the SNR covered. The default is 2000.

    Returns
    -------
    normCov : (np.array)
        Aggregate normalized coverage around the 3' end of transcripts.
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_TranCov):
        normCov = loadPKL(out_TranCov)
        return normCov
    
    # Get the expected coverage
    expCov = getExpectedCoverage(
        out_db,
        out_strandedFeats,
        exonic_bamfile,
        concordant = True,
        baselineFeat = 'transcript'
        )
    
    # In most cases, this should just load a file
    flatStrandedFeats = getStrandedFeats(
        out_strandedFeats,
        out_db,
        ('transcript', )
        )
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(exonic_bamfile, 'rb')
    
    # Initiate the variables needed
    transCount = 0
    coverage = np.zeros((window), 'L')
    
    for strd in (True, False):
        print(
            'Going over coverage of {}transcripts'.format('+' if strd else '-')
            )
        for feat in flatStrandedFeats[strd]['transcript']:
            if strd:
                mid = feat[2]
            else:
                mid = feat[1]
            start = int(mid - window / 2)
            stop = int(mid + window / 2)
            # Include correctins for the start & end if the window
            #  falls out of the reference size range
            if start < 0:
                corrStart = 0
            else:
                corrStart = start
            refLen = bam.get_reference_length(feat[0])
            if stop > refLen:
                corrStop = refLen
            else:
                corrStop = stop
            # Get the coverage summed over A/T/C/G; count only concordant reads
            refCoverage = np.sum(
                bam.count_coverage(
                    contig = feat[0],
                    start = corrStart,
                    stop = corrStop,
                    quality_threshold = 0,
                    read_callback = lambda r: r.is_reverse != strd
                    ),
                axis = 0
                )
            transCount += 1
            # If the window was out of the ref range, fill in the rest
            if corrStart == 0:
                refCoverage = np.append(
                    np.zeros((0 - start), 'L'),
                    refCoverage
                    )
            if corrStop != stop:
                refCoverage = np.append(
                    refCoverage,
                    np.zeros((stop - corrStop), 'L')
                    )
            # If needed, flip the coverage to be in the 5'->3'
            #  orientation wrt/ the transcript feature direction
            if not strd:
                refCoverage = refCoverage[::-1]
            # Add the SNR
            coverage += refCoverage
    # Normalize the coverage only once at the end
    normCov = coverage / (transCount * expCov)
    
    bam.close()
    savePKL(out_TranCov, normCov)
    
    return normCov


def integrateEnding(nEnd, oldEnds):
    """Helper function to merge overlapping endings, similar to integrateFeats. 

    Parameters
    ----------
    nEnd : (tuple)
        ( start, end )
    oldEnds : (list)
        [ (start, end ) ]

    Returns
    -------
    oldEnds : (list)
        [ ( start, end ) ] such that the endings do not overlap.
    """
    
    # Go over the list of previously saved endings, check for overlap and save
    #  all the overlaps found in a list.
    overlaps = []
    for oEnd in oldEnds:
        # If (the new End's beginning is within the old End)
        #  or (the new End's end is within the old End)
        #  or (the new End spans the old End)
        #  ), add old End to the overlaps list
        if (nEnd[0] >= oEnd[0] and nEnd[0] <= oEnd[1]) \
            or (nEnd[1] >= oEnd[0] and nEnd[1] <= oEnd[1]) \
            or (nEnd[0] < oEnd[0] and nEnd[1] > oEnd[1]):
                overlaps.append(oEnd)
    #  If overlaps have been found, merge with the new one
    if overlaps != []:
        # Initialize the start & end of the merged ending using the new End
        #  (not in the overlaps list)
        start = nEnd[0]; end = nEnd[1]
        # Go over all the overlaps & update the start & end as necessary
        for e in overlaps:
            if e[0] < start:
                start = e[0]
            if e[1] > end:
                end = e[1]
            # When done considering this ending, remove it from the master list
            oldEnds.remove(e)
        # When done, add the new merged ending to the list
        oldEnds.append((start, end))
    # If no overlaps have been found, simply add the new ending
    else:
        oldEnds.append(nEnd)
    
    return sorted(oldEnds)


def removeOverlaps(SNRendings, Tendings):
    """Helper function to remove portions of SNRendings (exons) that are
    overlapped by Tendings (endpieces). All 0-based.

    Parameters
    ----------
    SNRendings : list
        [ ( start, end ) ]
    Tendings : list
        [ ( start, end ) ]

    Returns
    -------
    SNRendings : list
        [ ( start, end ) ]
    """

    newEndings = []
    for SNRe in SNRendings:
        start = SNRe[0]
        end = SNRe[1]
        # Look for any overlap with Tendings
        for Te in Tendings:
            # Treat each case of overlap separately:
            # Complete overlap / both start & end within
            if start >= Te[0] and end <= Te[1]:
                break
            # Shorten start if only start within
            elif start >= Te[0] and start <= Te[1]:
                start = Te[1]
            # Shorten end if only end within
            elif end >= Te[0] and end <= Te[1]:
                end = Te[0]
            # If a Tending is entirely within an SNRending, treat recursively,
            #  with the two overhangs each as a new separate SNRending
            elif start < Te[0] and end > Te[1]:
                newEndings.extend(
                    removeOverlaps([(start, Te[0]), (Te[1], end)], Tendings)
                    )
                break
        # If none found, just add as is
        else:
            newEndings.append((start, end))

    return sorted(newEndings)


def getNonCanCovGenes(
        out_NonCanCovGenes,
        lenToSNRs,
        out_db,
        bamfile,
        lastBP = 250,
        covThreshold = 0.05,
        minLen = 0
        ):
    """Function to scan each gene with SNRs and determine what proportion of
    its non-canonical RNA-seq coverage may be accounted for by SNRs.

    Parameters
    ----------
    out_NonCanCovGenes : (str)
        Path to the output of this funciton.
    lenToSNRs : (dict)
        { SNRlength : [ SNR ] }
    out_db : (str)
        Path to the saved GTF/GFF database.
    bamfile : (str)
        Path to the (filtered) bam file.
    lastBP : (int), optional
        The extent to which to look for coverage at the exon-wise end of the
        transcript. The default is 250.
    covThreshold : (float), optional
        The minimal proportion of non-canonical coverage accounted for by SNRs
        to include the gene. The default is 0.05.
    minLen : int, optional
        Minimal SNR length that needs to be present for a gene to be
        considered. The default is 0.

    Returns
    -------
    nonCanCovGenes : (list)
        Genes found to have a significant proportion of coverage outside of
        the canonically expected area.
    """
    # Note that all except gffutils feats is 0-based    
    # If the file already exists, simply load
    if os.path.isfile(out_NonCanCovGenes):
        nonCanCovGenes = loadPKL(out_NonCanCovGenes)
        return nonCanCovGenes
    
    nonCanCovGenes = []
    # Sort SNRs by gene
    geneLenToSNRs = getGeneLenToSNRs(lenToSNRs, concordant = True)
    db = gffutils.FeatureDB(out_db, keep_order = True)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    progressTracker = 0
    prog = 0
    # For each gene that has an SNR detected
    for geneID, glenToSNRs in geneLenToSNRs.items():
        newProg = round(progressTracker / len(geneLenToSNRs), 3)
        if newProg != prog:
            prog = newProg
            print('Looking for genes with non-canonical'\
                  ' coverage... ({:.1%})'.format(prog))
            clear_output(wait = True)
        progressTracker += 1
        # If no SNRs of sufficient length are found, skip
        if all(length < minLen for length in glenToSNRs.keys()):
            continue
        # Get the total bp-wise coverage for the gene
        geneFeat = db[geneID]
        gene0start = geneFeat.start - 1     # Conversion to 0-based
        gene0end = geneFeat.end             # Conversion to 0-based
        geneBPcov = np.sum(
            bam.count_coverage(
                contig = geneFeat.seqid,
                start = gene0start,
                stop = gene0end,
                quality_threshold = 0,
                read_callback = lambda r:
                    r.is_reverse == (geneFeat.strand == '-')
                ),
            axis = 0
            )
        # If there's no gene coverage, skip
        geneCoverage = sum(geneBPcov)
        if geneCoverage == 0:
            continue
        transcripts = []
        Tendings = []
        strd = geneFeat.strand == '+'
        # Go over ALL transcripts overlapping this gene (even from other genes)
        #  add up the coverage of exon-wise last X bp
        for trans in db.features_of_type(
                featuretype = 'transcript',
                limit = (geneFeat.seqid, gene0start, gene0end),
                strand = geneFeat.strand
                ):
            # Extract the exons as (start, end); 1- => 0-based
            exons = sorted(
                [(e.start - 1, e.end) for e in db.children(
                    trans,
                    featuretype = 'exon'
                    )],
                reverse = strd      # Order from last to 1st wrt/ the strand
                )
            trans0start = trans.start - 1       # Conversion to 0-based
            trans0end = trans.end               # Conversion to 0-based
            [geneID] = trans.attributes['gene_id']
            # Save all transcripts (lite) on the same strand as the gene
            transcripts.append(
                Transcript(geneID, trans.seqid, trans0start, trans0end, exons)
                )
            # Determine the transcripts's exon-wise end
            remaining = lastBP
            covStart = None
            # Go over exons from last to first wrt/ the strand
            for exon in exons:
                exLen = exon[1] - exon[0]
                if exLen < remaining:
                    remaining -= exLen
                else:
                    covStart = exon[1] - remaining if strd \
                        else exon[0] + remaining
                    break
            # If we have exhausted all exons, take the beginning of the tran
            else:
                covStart = trans0start if strd else trans0end
            # Only add the ending if the transcripts's exon-wise end is not
            #  beyond the gene end (in case it is from another gene)
            if strd and covStart < gene0end:
                newEnding = (
                    max(covStart, gene0start),
                    min(trans0end, gene0end)
                    )
            elif not strd and covStart > gene0start:
                newEnding = (
                    max(trans0start, gene0start),
                    min(covStart, gene0end)
                    )
            else:
                continue
            Tendings = integrateEnding(newEnding, Tendings)
        # Get the coverage of transcripts endings
        endCoverage = 0
        for covS, covE in Tendings:
            endCoverage += np.sum(
                bam.count_coverage(
                    contig = geneFeat.seqid,
                    start = covS,
                    stop = covE,
                    quality_threshold = 0,
                    read_callback = lambda r:
                        r.is_reverse == (geneFeat.strand == '-')
                    )
                )
        totalSNRcovProp = 0
        lenToCovSNRs = {}
        # Extract and merge the SNR endings by length
        for length, snrs in glenToSNRs.items():
            SNRendings = []
            if length >= minLen:
                for snr in snrs:
                    if strd and snr.start > gene0start:
                        newSNRe = (
                            max(snr.start - lastBP, gene0start),
                            min(snr.start, gene0end)
                            )
                    elif not strd and snr.end < gene0end:
                        newSNRe = (
                            max(snr.end, gene0start),
                            min(snr.end + lastBP, gene0end)
                            )
                    else:
                        continue
                    # Add the new ending to the list without overlap
                    SNRendings = integrateEnding(newSNRe, SNRendings)
                # Remove the portions of SNRendings overlapped by Tendings
                SNRendings = removeOverlaps(SNRendings, Tendings)
                # Get the SNR coverage attributable to SNRs of this length
                snrCoverage = 0
                for covS, covE in SNRendings:
                    snrCoverage += np.sum(
                        bam.count_coverage(
                            contig = geneFeat.seqid,
                            start = covS,
                            stop = covE,
                            quality_threshold = 0,
                            read_callback = lambda r:
                                r.is_reverse == (geneFeat.strand == '-')
                            )
                        )
                prop = snrCoverage / (geneCoverage - endCoverage)
                totalSNRcovProp += prop
                lenToCovSNRs[length] = (prop, snrs)
        # Add this as an interesting gene only if the total proportion of
        #  non-canonical coverage attributable to SNRs exceeds the threshold
        #  Note that the SNR coverage between lengths may overlap, so this is
        #  only an estimate
        if totalSNRcovProp >= covThreshold:
            nonCanCovGenes.append(
                NoncGene(
                    geneFeat, geneBPcov, transcripts, lenToCovSNRs, Tendings
                    )
                )
    bam.close()
    savePKL(out_NonCanCovGenes, nonCanCovGenes)
    
    return nonCanCovGenes


def getBaselineData(out_db, bamfile):
    """Function to obtain the baseline data as an input for calculation of
    sensitivity, specificity, and accuracy. Note that here, only expressed
    (i.e., with any coverage) transcripts are considered.
    
    Parameters
    ----------
    out_db : (str)
        Path to the saved GTF/GFF database.
    bamfile : (str)
        Path to the (filtered) bam file.

    Returns
    -------
    expectedExCov : float
        Expected coverage across exons flattened from expressed transcripts.
    coveredTranscripts : list
        [ Transcript ] of transcripts with non-0 coverage
    """
    
    global expectedExCov, covTransByStrdRef, covExonsByStrdRef
    
    # If done before, simply return the existing results; otherwise process
    if expectedExCov != 0:
        return expectedExCov, covTransByStrdRef, covExonsByStrdRef
    
    db = gffutils.FeatureDB(out_db, keep_order = True)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    references = bam.references
    totalExCov = 0
    totalExLen = 0
    
    for strd in (True, False):
        for refname in references:
            # Initiate the flattened exon list for this reference
            for trans in db.features_of_type(
                featuretype = 'transcript',
                limit = (refname, 1, bam.get_reference_length(refname)),
                strand = '+' if strd else '-'
                ):
                # Process this transcript only if it has ANY coverage
                transCov = np.sum(
                    bam.count_coverage(
                        contig = trans.seqid,
                        start = trans.start - 1,      # 1-based => 0-based
                        stop = trans.end,
                        quality_threshold = 0,
                        read_callback = lambda r: r.is_reverse != strd
                        )
                    )
                # Only transcripts with ANY coverage are considered here
                if transCov == 0:
                    continue
                # Extract exons from the transcript, last-to-first
                exons = sorted(
                    [(e.start - 1, e.end) for e in db.children(
                        trans,
                        featuretype = 'exon'
                        )],
                    reverse = strd
                    )
                trans0start = trans.start - 1       # Conversion to 0-based
                trans0end = trans.end               # Conversion to 0-based
                [geneID] = trans.attributes['gene_id']
                # Save this transcript to the list of covered transcripts
                covTransByStrdRef[strd][refname].append(
                    Transcript(geneID, trans0start, trans0end, exons)
                    )
                # Integrate these exons in the ref-wide flattened exon list
                for exon in exons:
                    covExonsByStrdRef[strd][refname] = integrateEnding(
                        exon,
                        covExonsByStrdRef[strd][refname]
                        )
            # Get the total coverage & length of this ref's exons for baseline
            #  Note: only exons from transcripts w/ ANY coverage are considered
            for ex in covExonsByStrdRef[strd][refname]:
                totalExLen += ex[1] - ex[0]
                totalExCov += np.sum(
                    bam.count_coverage(
                        contig = refname,
                        start = ex[0],      # already 0-based
                        stop = ex[1],
                        quality_threshold = 0,
                        read_callback = lambda r: r.is_reverse != strd
                        ),
                    dtype = 'int'
                    )
    
    # Get the expected coverage
    expectedExCov = totalExCov / totalExLen
    
    return expectedExCov, covTransByStrdRef, covExonsByStrdRef
    

def getTransEndSensSpec(
        out_db,
        bamfile,
        endLenMin = 50,
        endLenMax = 250,
        results = {}
        ):
    """Function to provide data for a ROC curve for various transcript end
    lengths.
    
    Parameters
    ----------
    out_db : (str)
        Path to the saved GTF/GFF database.
    bamfile : (str)
        Path to the (filtered) bam file.
    endLenMin : (int), optional
        The minimal length of exon-wise transcript end to consider.
        The default is 50.
    endLenMax : (int), optional
        The minimal length of exon-wise transcript end to consider.
        The default is 150.
    results : (dict)
        { endLen : (sensitivity, specificity, accuracy) }. The default is {}.

    Returns
    -------
    results : (dict)
        { endLen : (sensitivity, specificity, accuracy) }
    """
    # Notes:
    #
    # Sensitivity = percentage of total exonic coverage captured by the
    #  transcript ends of a given length
    #  - Get total read coverage of all flattened exons
    #  - Get total read coverage of all exon-wise transcript ends
    #  = ENDS / TOTAL
    # Specificity = percentage of non-covered exon area not overlapped
    #  by transcript ends
    #  - Sum the number of flattened exon bases without coverage
    #  - Sum the number of flattened exon bases that do not overlap ends
    #  = NON-ENDS / NON-COVERAGE
    # 
    # Do this for multiple end lengths & make ROC
    # !!! Only consider genes with coverage (expressed)!
    #
    # In order to get accuracy, the following is needed: TP, FP, TN, FN
    # This is binarized wrt/ the expected coverage: >= / <
    # This information also can be used to calculate Sens and Spec
    #
    # The 
    
    # Get the baseline variables
    expectedExCov, covTransByStrdRef, covExonsByStrdRef = getBaselineData(
        out_db,
        bamfile
        )
    
    bam = pysam.AlignmentFile(bamfile, 'rb')
    references = bam.references
    
    # Initiate progress tracking
    refCounter = 0
    totalRefs = len(references) * 6     # For each strand (2) and length (3)
    # Do the end measurements for each potential end lengths min, mid, max to
    #  get the flattened transcript exon-wise end pieces across each ref
    for endLength in (
            endLenMin,
            np.mean((endLenMin, endLenMax), dtype = int),
            endLenMax
            ):
        # "POSITIVE" ~ Inside exon-wise transcript ends
        TP = 0  # Number of bp >= exp in exon-wise ends
        FP = 0  # Number of bp < exp in exon-wise ends
        # "NEGATIVE" ~ Outside exon-wise transcript ends
        TN = 0  # Number of bp < exp outside exon-wise ends
        FN = 0  # Number of bp >= exp outside exon-wise ends
        # Go over each strand separately
        for strd in (True, False):
            # Work for each reference separately
            for refname in references:
                # Track progress & announce in the output
                refCounter += 1
                print('Processing transcript ends of length {}...'\
                      ' ({:.2%})'.format(endLength, refCounter / totalRefs))
                clear_output(wait = True)
                refEndPieces = []
                # If not extracted before, initialize exons for this ref
                if baselineData is None:
                    refExons = []
                # Iterate through each transcript
                for trans in db.features_of_type(
                    featuretype = 'transcript',
                    limit = (refname, 1, bam.get_reference_length(refname)),
                    strand = '+' if strd else '-'
                    ):
                    # Process this transcript only if it has ANY coverage
                    transCov = np.sum(
                        bam.count_coverage(
                            contig = trans.seqid,
                            start = trans.start - 1,      # 1-based => 0-based
                            stop = trans.end,
                            quality_threshold = 0,
                            read_callback = lambda r: r.is_reverse != strd
                            )
                        )
                    if transCov == 0:
                        continue
                    # Extract exons from the transcript, last-to-first
                    exons = sorted(
                        [(e.start - 1, e.end) for e in db.children(
                            trans,
                            featuretype = 'exon'
                            )],
                        reverse = strd
                        )
                    # If not done before, add all relevant exons to the
                    #  flattened ref exon list
                    if baselineData is None:
                        for exon in exons:
                            refExons = integrateEnding(exon, refExons)
                    # Extract & save exon-wise endpieces for this transcript
                    remaining = endLength
                    # For each exon, check whether it is longer than the
                    #  remaining endLength
                    for exon in exons:
                        exLen = exon[1] - exon[0]
                        # If not, add the entire exon
                        if exLen < remaining:
                            refEndPieces = integrateEnding(exon, refEndPieces)
                            remaining -= exLen
                        # Otherwise add only the remaining piece
                        else:
                            refEndPieces = integrateEnding(
                                (exon[1] - remaining, exon[1]) if strd \
                                    else (exon[0], exon[0] + remaining),
                                refEndPieces
                                )
                            break
                
                # For each reference (on each strand), extract and add the
                #  total exon coverage and non-covBP (if not done before),
                #  the end coverage, and the non-endBP
                if baselineData is None:
                    # Save the refExons, if not done before
                    exonsByRef[refname] = refExons
                    for ex in refExons:
                        # Get bp-wise coverage for all flattened exons
                        exCov = np.sum(
                            bam.count_coverage(
                                contig = refname,
                                start = ex[0],      # already 0-based
                                stop = ex[1],
                                quality_threshold = 0,
                                read_callback = lambda r: r.is_reverse != strd
                                ),
                            axis = 0
                            )
                        totalExCov += np.sum(exCov)
                        totalExNonCovBP += np.count_nonzero(exCov == 0)
                # Get the exon-wise end piece coverage
                for refEndpiece in refEndPieces:
                    endsCov += np.sum(
                        bam.count_coverage(
                            contig = refname,
                            start = refEndpiece[0],      # already 0-based
                            stop = refEndpiece[1],
                            quality_threshold = 0,
                            read_callback = lambda r: r.is_reverse != strd
                            )
                        )
                # Remove the portions of flattened exons overlapped by
                #  endpieces and save as starts
                refStarts = removeOverlaps(exonsByRef[refname], refEndPieces)
                # Get the number of non-covered BPs across exon-wise starts
                for refStart in refStarts:
                    startsCov = np.sum(
                        bam.count_coverage(
                            contig = refname,
                            start = refStart[0],      # already 0-based
                            stop = refStart[1],
                            quality_threshold = 0,
                            read_callback = lambda r: r.is_reverse != strd
                            ),
                        axis = 0
                        )
                    nonCovStartBP += np.count_nonzero(startsCov == 0)
                
            # Save the exons by ref for this strand, if not done before
            if baselineData is None:
                exonsByStrdRef[strd] = exonsByRef
        # When done going over both strands, if not done before, save the
        #  baseline data so that this does not have to be done again for other
        #  end lengths
        if baselineData is None:
               baselineData = (totalExCov, totalExNonCovBP, exonsByStrdRef)
        # Save the results for each length
        results[endLength] = (
            endsCov / totalExCov,               # Sensitivity
            nonCovStartBP / totalExNonCovBP     # Specificity
            )
    bam.close()

    return results