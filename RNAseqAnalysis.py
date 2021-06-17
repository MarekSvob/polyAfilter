#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 17:00:37 2020

@author: marek
"""
import os
import pysam
import gffutils
import logging
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from IPython.display import clear_output
from functools import partial
from sklearn import metrics

from SNRanalysis import getFlatFeatsByTypeStrdRef, getSNRsByGeneLen, \
    flattenIntervals, getOverlaps, removeOverlaps
from SNRdetection import savePKL, loadPKL


logging.basicConfig(level = logging.INFO,
                    format = '%(asctime)s - %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


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
            self.geneFeat.id, self.geneFeat.start, self.geneFeat.end,
            self.geneFeat.strand)

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
        return '{} Transcript @{}-{}({})'.format(self.geneID, self.start,
                                                 self.end, self.strand)


def filterReads(filt_bamfile, out_db, out_strandedFeats, bamfile,
                featType = 'exon', concordant = True):
    """Function to create and index a bam file from a subset of reads that are
    aligned to a specific feature type on the same or opposite strand.

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
        logger.info('The filtered bam file already exists.')
        if not os.path.isfile('{}.bai'.format(filt_bamfile)):
            pysam.index(filt_bamfile)
            logger.info('The bam file has been indexed.')
        return    
    
    # This will usually just load the file
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats,
                                          out_db,
                                          (featType, ))
    
    # Connect to the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # Initiate the variables needed
    # This needs to start off as a set to avoid double-adding reads aligned to
    #  two features of the same type
    filt_reads = set()
    
    # Get only reads on the same/other strand (depends if concordant or not)
    #  as the feature in question
    for strd, featsByStrdRef in flatFeats[featType].items():
        # Cycle over the feats under the appropriate key
        for ref, feats in featsByStrdRef.items():
            logger.info(f'Processing {"+" if strd else "-"} strand of '
                        f'reference {ref}...')
            for start, end in feats:
                for read in bam.fetch(contig = ref, start = start, stop = end):
                    # Make sure the read is on the correct strand before adding
                    if concordant == (read.is_reverse != strd):
                        filt_reads.add(read)
    
    # Transform the set into a list for sorting
    filt_reads_list = list(filt_reads)
    # Reads must be sorted for indexing
    filt_reads_list.sort(
        key = lambda x:(x.reference_id, x.reference_start, x.reference_length))
    
    # Create the bamfile to add the reads
    filt_bam = pysam.AlignmentFile(filt_bamfile, 'wb', template = bam)
    # Add the reads in order
    for read in filt_reads_list:
        filt_bam.write(read)
    # Close the files
    filt_bam.close()
    bam.close()
    
    pysam.index(filt_bamfile)
    logger.info('A filtered bam file has been created and indexed.')
    

def getExpectedCoverage(out_db, out_strandedFeats, bamfile, concordant = True,
                        baselineFeat = 'transcript'):
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
    
    # In most cases, this should just load a file
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db,
                                          (baselineFeat, ))
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')

    totalCoverage = 0
    totalLength = 0

    # Get the total length of the feats and the total coverage over them
    for strd, featsByRef in flatFeats[baselineFeat].items():
        logger.info(f'Scanning {"+" if strd else "-"}{baselineFeat}s for total'
                    ' coverage and length...')
        for ref, feats in featsByRef.items():
            for start, end in feats:
                totalLength += end - start
                totalCoverage += np.sum(
                    bam.count_coverage(
                        contig = ref, start = start, stop = end,
                        quality_threshold = 0,
                        read_callback = lambda r:
                            r.is_reverse != (strd == concordant)),
                    dtype = 'int')
    bam.close()
    
    # Normalize the coverage by #SNRs, total coverage, and total length
    expected = totalCoverage / totalLength
    
    return expected


def getCovPerSNR(out_SNRCovLen, out_SNROutl, lenToSNRs, out_db,
                 out_strandedFeats, exonic_bamfile, SNRlengths = range(200),
                 window = 2000, concordant = True, SNRfeat = 'transcript',
                 log10Filter = 3):
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
        out_db, out_strandedFeats, exonic_bamfile, concordant = concordant,
        baselineFeat = SNRfeat)
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(exonic_bamfile, 'rb')
    
    # Initialize the dict to collect all the data
    covByLen = {}
    outlierPeaks = defaultdict(list)
    totalCount = 0
    filteredOut = 0
    zeroCovs = 0
    # Go over the SNRs by length, starting with the longest
    for length in sorted(lenToSNRs.keys(), reverse = True):
        if length in SNRlengths:
            # For each length, initialize the count & array for window coverage
            logger.info('Going over '
                        f'{"concordant" if concordant else "discordant"} '
                        f'coverage of SNR{length}s...')
            SNRcount = 0
            zeros = 0
            coverage = np.zeros((window), 'L')
            # Note that SNRs are 0-based
            for SNR in lenToSNRs[length]:
                # Only consider SNRs in transcripts, as those are the
                #  sources of internal (concordant or discordant) priming
                if ((concordant and SNRfeat in SNR.concFeats)
                        or (not concordant and SNRfeat in SNR.discFeats)):
                    # Get the mid point of the SNR wrt/ the reference
                    mid = (SNR.end + SNR.start) // 2
                    start = int(mid - window / 2)
                    stop = int(mid + window / 2)
                    # Include corrections for the start & end if the window
                    #  falls out of the reference size range
                    refLen = bam.get_reference_length(SNR.record)
                    corrStart = max(0, start)
                    corrStop = min(refLen, stop)
                    # Get the coverage summed over A/T/C/G; count only
                    #  reads on the same (conc) or opposite (disc) strand
                    refCov = np.sum(
                        bam.count_coverage(contig = SNR.record,
                                           start = corrStart,
                                           stop = corrStop,
                                           quality_threshold = 0,
                                           read_callback = lambda r:
                                               r.is_reverse != SNR.strand),
                        axis = 0)
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
                                np.zeros((corrStart - start), 'L'), refCov)
                        if corrStop != stop:
                            refCov = np.append(
                                refCov, np.zeros((stop - corrStop), 'L'))
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
    logger.info(f'Filtered out {filteredOut / totalCount:.2%} outliers and '
                f'{zeroCovs / totalCount:.2%} SNRs had no coverage.')
    
    return covByLen, outlierPeaks


def getCovPerSNRlen(out_SNRCovLen, fasta, out_db, out_strandedFeats, bamfile,
                    base = 'A', minSNRlen = 5, mism = 0, window = 2000,
                    concordant = True, SNRfeat = 'transcript'):
    """Obtain RNA-seq coverage in the vicinity of SNRs of given length by
    scanning the reference genome.
    
    Parameters
    ----------
    out_SNRCovLen : (str)
        Path to the file with per-SNR coverage by length.
    fasta : (str)
        Path to the reference fasta file.
    out_db : (str)
        Path to the file with the database.
    out_strandedFeats : (str)
        Path to the file with stranded feats.
    bamfile : (str)
        Path to the (pre-filtered) bamfile.
    base : (str), optional
        DNA base constituting SNRs to be searched for: 'A', 'T', 'C', 'G', 'N'.
        The default is 'A'.
    minSNRlen : (int), optional
        The shortest SNR length to consider. The default is 5.
    mism : (int), optional
        The maximum number of mismatches in an SNR allowed. The default is 0.
    window : (int), optional
        Total size of the window around the SNR covered. The default is 2000.
    concordant : (bool), optional
        Switch between concordant & discordant coverage. The default is True.
    SNRfeat : (str), optional
        The feature by which SNRs are selected and to which the read depth will
        be normalized. The default is 'transcript'.

    Returns
    -------
    normCovByLen : (dict)
        { length : ( SNRcount, SNRzeroCount, np.array ) }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_SNRCovLen):
        normCovByLen = loadPKL(out_SNRCovLen)
        return normCovByLen
    
    def mayAddCoverage():
        """A helper function that adds an SNR coverage if the length condition
        is met. Similar to the helper function in findSNRsWmisms() but adds SNR
        coverage instead of SNRs.
        """
        nonlocal covByLen, lastBlockSaved, totalCount, zeroCovs
        
        # Add the feature start to work with reference coordinates
        first = blocks[0][0] + fStart
        last = blocks[-1][-1] + fStart
        length = last - first
        # If minSNRlen has been met, add the SNR coverage into the output dict
        #  and mark the last block as saved
        if length >= minSNRlen:
            # Determine the range of the coverage sought
            edge = first if strd else last
            start = int(edge - window / 2)
            stop = int(edge + window / 2)
            # Include corrections for the start & end if the window falls out
            #  of the reference size range
            corrStart = max(0, start)
            corrStop = min(stop, refLen)
            # Get the coverage summed over A/T/C/G; count only reads on the
            #  same (conc) or opposite (disc) strand
            refCov = np.sum(
                bam.count_coverage(contig = ref,
                                   start = corrStart,
                                   stop = corrStop,
                                   quality_threshold = 0,
                                   read_callback = lambda r:
                                       r.is_reverse != (strd == concordant)),
                axis = 0)
            # If this is zero coverage, only add the counts
            if sum(refCov) == 0:
                zeroCovs += 1
                # If the length does not yet exist, initiate, otherwise add
                if length not in covByLen.keys():
                    covByLen[length]['SNRcount'] = 1
                    covByLen[length]['zeros'] = 1
                    covByLen[length]['coverage'] = np.zeros(window, 'L')
                else:
                    covByLen[length]['SNRcount'] += 1
                    covByLen[length]['zeros'] += 1
            # Otherwise add the coverage found
            else:                    
                # If the window was out of ref range, fill in the rest with 0s
                if corrStart != start:
                    refCov = np.append(np.zeros((corrStart - start), 'L'),
                                       refCov)
                if corrStop != stop:
                    refCov = np.append(refCov,
                                       np.zeros((stop - corrStop), 'L'))
                # If needed, flip the coverage to be in the 5'->3' orientation
                #  wrt/ the transcript feature
                refCov = refCov if strd else refCov[::-1]
                # If the length does not yet exist, initiate, otherwise add
                if length not in covByLen.keys():
                    covByLen[length]['SNRcount'] = 1
                    covByLen[length]['zeros'] = 0
                    covByLen[length]['coverage'] = refCov
                else:
                    covByLen[length]['SNRcount'] += 1
                    covByLen[length]['coverage'] += refCov
            # Add to the total and confirm
            totalCount += 1
            lastBlockSaved = True
    
    # Handling of non-caps in the sequence
    capsDict = {'A':{'a', 'A'},
                'T':{'t', 'T'},
                'C':{'c', 'C'},
                'G':{'g', 'G'},
                'N':{'n', 'N'}}
    # Handling of complementary bases for '-' strd
    compDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    
    # Load the file with stranded feats
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db,
                                          (SNRfeat, ))
    # Get the expected coverage
    expCov = getExpectedCoverage(out_db, out_strandedFeats, bamfile,
                                 concordant = concordant,
                                 baselineFeat = SNRfeat)
    # Connect to the fasta reference
    recDict = SeqIO.index(fasta, 'fasta')
    # Attach the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # Initialize the dict and counts to collect all the data
    covByLen = defaultdict(dict)
    totalCount = 0
    zeroCovs = 0
    # Go over the feats by strand and ref
    for strd, featsByRef in flatFeats[SNRfeat].items():
        # Set the base of interest based on the strand
        boi = base if strd else compDict[base]
        # Allow for either caps to be included
        eitherCaps = capsDict[boi]
        for ref, feats in featsByRef.items():
            logger.info('Going over '
                        f'{"concordant" if concordant else "discordant"} '
                        f'coverage around SNR{minSNRlen}+ on {SNRfeat}s on '
                        f'{ref}{"+" if strd else "-"}...')
            seq = str(recDict[ref].seq)
            refLen = len(seq)
            for fStart, fEnd in feats:
                # Initiate the mismatch counter
                nMism = 0
                # Initiate the list of valid base blocks
                blocks = []
                # Initiate the indicator of being in a block
                blocked = False
                # Initiate the indicator of having saved an SNR coverage
                #  with the last block
                lastBlockSaved = True
                # Scanning the sequence, while adding and removing blocks
                featSeq = seq[fStart:fEnd]
                for i, bp in enumerate(featSeq):
                    # If this is boi but block has not been started, start one
                    #  (if currently in a block, just pass onto the next base
                    #  - the two conditions below cannot be merged into one)
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
                                # Only save the SNR if the most recent block
                                #  hasn't been saved in an SNR yet (if it has,
                                #  it implies that this would not be the
                                #  longest SNR possible)
                                if not lastBlockSaved:
                                    mayAddCoverage()
                                # Deduct all mismatches before the 2nd block
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
                    mayAddCoverage()
    bam.close()
    
    # Normalize by the expected coverage for the given feature
    normCovByLen = {}
    for length, covDict in covByLen.items():
        normCovByLen[length] = {'SNRcount': covDict['SNRcount'],
                                'zeros': covDict['zeros'],
                                'coverage': covDict['coverage'] / expCov}
    
    logger.info(f'Detected coverage for {totalCount:,d} SNRs on {SNRfeat}s, '
                f'of which {zeroCovs:,d} ({zeroCovs / totalCount:.2%}) '
                'had no coverage.')
    savePKL(out_SNRCovLen, normCovByLen)
    
    return normCovByLen


def getCovPerTran(out_TranCov, out_db, out_strandedFeats, exonic_bamfile,
                  window = 2000, concordant = True, normFeat = 'exon'):
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
    concordant : (bool), optional
        Determines which strand the alignments come from (fwd vs. rev).
    normFeat : (str), optional
        The feature where SNRs will be sought and to which the coverage will be
        normalized.

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
    expCov = getExpectedCoverage(out_db, out_strandedFeats, exonic_bamfile,
                                 concordant = concordant,
                                 baselineFeat = normFeat)
    
    # In most cases, this should just load a file
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db,
                                          ('transcript', ))
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(exonic_bamfile, 'rb')
    
    # Initiate the variables needed
    transCount = 0
    coverage = np.zeros((window), 'L')
    
    for strd, featsByRef in flatFeats['transcript'].items():
        logger.info(f'Going over coverage of {"+" if strd else "-"}'
                    'transcripts...')
        for ref, feats in featsByRef.items():
            refLen = bam.get_reference_length(ref)
            for feat in feats:
                if strd:
                    mid = feat[1]
                else:
                    mid = feat[0]
                start = int(mid - window / 2)
                stop = int(mid + window / 2)
                # Include corrections for the start & end if the window
                #  falls out of the reference size range
                corrStart = max(0, start)
                corrStop = min(refLen, stop)
                # Get the coverage summed over A/T/C/G; count only concordant
                refCoverage = np.sum(
                    bam.count_coverage(contig = ref,
                                       start = corrStart,
                                       stop = corrStop,
                                       quality_threshold = 0,
                                       read_callback = lambda r:
                                           r.is_reverse != (
                                               strd == concordant)),
                    axis = 0)
                transCount += 1
                # If the window was out of the ref range, fill in the rest
                if corrStart == 0:
                    refCoverage = np.append(np.zeros((0 - start), 'L'),
                                            refCoverage)
                if corrStop != stop:
                    refCoverage = np.append(refCoverage,
                                            np.zeros((stop - corrStop), 'L'))
                # If needed, flip the coverage to be in the 5'->3'
                #  orientation wrt/ the transcript feature direction
                if not strd:
                    refCoverage = refCoverage[::-1]
                # Add the transcript coverage
                coverage += refCoverage
    # Normalize the coverage only once at the end
    normCov = coverage / (transCount * expCov)
    
    bam.close()
    savePKL(out_TranCov, normCov)
    
    return normCov


def getNonCanCovGenes(out_NonCanCovGenes, lenToSNRs, out_db, bamfile,
                      lastBP = 250, covThreshold = 0.05, minLen = 0):
    """Function to scan each gene with SNRs and determine what proportion of
    its non-canonical RNA-seq coverage (exon-wise ends) may be accounted for
    by SNRs.

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
    geneLenToSNRs = getSNRsByGeneLen(lenToSNRs, concordant = True)
    db = gffutils.FeatureDB(out_db, keep_order = True)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    progressTracker = 0
    prog = 0
    nGenes = len(geneLenToSNRs)
    # For each gene that has an SNR detected
    for geneID, glenToSNRs in geneLenToSNRs.items():
        progressTracker += 1
        newProg = round(progressTracker / nGenes, 3)
        if newProg != prog:
            prog = newProg
            logger.info('Looking for genes with non-canonical coverage... '
                        f'({prog:.1%})')
            clear_output(wait = True)
        # If no SNRs of sufficient length are found, skip
        if all(length < minLen for length in glenToSNRs.keys()):
            continue
        # Get the total bp-wise coverage for the gene
        geneFeat = db[geneID]
        gene0start = geneFeat.start - 1     # Conversion to 0-based
        gene0end = geneFeat.end             # Conversion to 0-based
        geneBPcov = np.sum(
            bam.count_coverage(contig = geneFeat.seqid,
                               start = gene0start,
                               stop = gene0end,
                               quality_threshold = 0,
                               read_callback = lambda r:
                                   r.is_reverse == (geneFeat.strand == '-')),
            axis = 0)
        # If there's no gene coverage, skip
        geneCoverage = sum(geneBPcov)
        if geneCoverage == 0:
            continue
        transcripts = []
        Tendings = []
        strd = geneFeat.strand == '+'
        # Go over ALL transcripts overlapping this gene (even from other genes)
        #  on the same strand and add up the coverage of exon-wise last X bp
        for trans in db.features_of_type(
                featuretype = 'transcript',
                limit = (geneFeat.seqid, geneFeat.start, geneFeat.end), # 1-b
                strand = geneFeat.strand):
            # Extract the exons as (start, end); 1- => 0-based
            exons = sorted([(e.start - 1, e.end) for e
                            in db.children(trans, featuretype = 'exon')],
                           reverse = strd) # Order from last to 1st wrt/ strand
            trans0start = trans.start - 1       # Conversion to 0-based
            trans0end = trans.end               # Conversion to 0-based
            [t_geneID] = trans.attributes['gene_id']
            # Save all transcripts (lite) on the same strand as the gene
            transcripts.append(Transcript(t_geneID, trans0start, trans0end,
                                          exons))
            # Determine the transcript's exon-wise end
            remaining = lastBP
            # Go over exons from last to first wrt/ the strand
            for exon in exons:
                exLen = exon[1] - exon[0]
                if exLen < remaining:
                    remaining -= exLen
                else:
                    covStart = (exon[1] - remaining if strd
                                else exon[0] + remaining)
                    break
            # If we have exhausted all exons, take the beginning of the tran
            else:
                covStart = trans0start if strd else trans0end
            # Only add the ending if the transcripts's exon-wise end is not
            #  entirely beyond the gene end (in case it is from another gene)
            if strd and covStart < gene0end:
                newEnding = (max(covStart, gene0start),
                             min(trans0end, gene0end))
            elif not strd and covStart > gene0start:
                newEnding = (max(trans0start, gene0start),
                             min(covStart, gene0end))
            else:
                continue
            Tendings.append(newEnding)
        # Flattend the endings
        Tendings = flattenIntervals(Tendings)
        # Get the coverage of transcripts endings
        endCoverage = 0
        for covS, covE in Tendings:
            endCoverage += np.sum(
                bam.count_coverage(
                    contig = geneFeat.seqid, start = covS, stop = covE,
                    quality_threshold = 0, read_callback = lambda r:
                        r.is_reverse == (geneFeat.strand == '-')))
        totalSNRcovProp = 0
        lenToCovSNRs = {}
        # Extract and merge the SNR endings by length
        for length, snrs in glenToSNRs.items():
            SNRendings = []
            if length >= minLen:
                for snr in snrs:
                    if strd and snr.start > gene0start:
                        newSNRe = (max(snr.start - lastBP, gene0start),
                                   min(snr.start, gene0end))
                    elif not strd and snr.end < gene0end:
                        newSNRe = (max(snr.end, gene0start),
                                   min(snr.end + lastBP, gene0end))
                    else:
                        continue
                    # Add the new ending to the list without overlap
                    SNRendings.append(newSNRe)
                # Flatten the SNRendings
                SNRendings = flattenIntervals(SNRendings)
                # Remove the portions of SNRendings overlapped by Tendings
                SNRendings = removeOverlaps(SNRendings, Tendings)
                # Get the SNR coverage attributable to SNRs of this length
                snrCoverage = 0
                for covS, covE in SNRendings:
                    snrCoverage += np.sum(
                        bam.count_coverage(
                            contig = geneFeat.seqid,
                            start = covS, stop = covE,
                            quality_threshold = 0,
                            read_callback = lambda r:
                                r.is_reverse == (geneFeat.strand == '-')))
                prop = snrCoverage / (geneCoverage - endCoverage)
                totalSNRcovProp += prop
                lenToCovSNRs[length] = (prop, snrs)
        # Add this as an interesting gene only if the total proportion of
        #  non-canonical coverage attributable to SNRs exceeds the threshold
        #  Note that the SNR coverage between lengths may overlap, so this is
        #  only an estimate
        if totalSNRcovProp >= covThreshold:
            nonCanCovGenes.append(NoncGene(geneFeat, geneBPcov, transcripts,
                                           lenToCovSNRs, Tendings))
    bam.close()
    savePKL(out_NonCanCovGenes, nonCanCovGenes)
    
    return nonCanCovGenes


def getBaselineData(out_transBaselineData, out_db, bamfile, includeIntrons,
                    weightedCov = True):
    """Function to obtain the baseline reference data (covered transcripts and
    total TP/TN) as an input for faster calculations of sensitivity &
    specificity at each given end length. Note that here, only expressed (i.e.,
    with any coverage) transcripts are considered.
    
    Parameters
    ----------
    out_transBaselineData : (str)
        Path to the saved output, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    bamfile : (str)
        Path to the (filtered) bam file.
    includeIntrons : (bool)
        Whether (positive or negative) coverage of introns should also be
        considered.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Alternatively, the sensitivity/specificity
        relationship is purely binary (flat) and the covered bases (1) have the
        same weight as non-covered ones (0). The default is True.
    
    Returns
    -------
    covTransByStrdRef : (dict)
        { strd : { ref : [ Transcript ] } } of transcripts with non-0 exonic
        coverage.
    Pos : (int)
        Total coverage across all flattened exons (if weightedCov, otherwise
        only number of exon bps with any coverage) from covered transcripts.
    Neg : (int)
        Total number of exon bps with no coverage from covered transcripts.
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_transBaselineData):
        covTransByStrdRef, Pos, Neg = loadPKL(
            out_transBaselineData)
        return covTransByStrdRef, Pos, Neg
    
    # Otherwise initialize the variables    
    logger.info('Getting baseline data across covered transcripts, '
                f'{"including" if includeIntrons else "excluding"} introns, '
                f'in {bamfile}...')
    covTransByStrdRef = defaultdict(lambda: defaultdict(list))
    Pos = 0
    Neg = 0
    
    db = gffutils.FeatureDB(out_db, keep_order = True)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    references = bam.references
    
    # Get all exon-covered transcripts and their flattened exons by strd & ref
    for strd in (True, False):
        for refname in references:
            logger.info('Processing transcripts on reference '
                        f'{refname}{"+" if strd else "-"}...')
            # Initiate the to-be-flattened list of exons/transcripts for this
            #  strd/reference
            covPieces = []
            # Process this transcript only if it has ANY (exonic) coverage
            for trans in db.features_of_type(
                    featuretype = 'transcript',
                    limit = (refname, 1, bam.get_reference_length(refname)),
                    strand = '+' if strd else '-'):
                # Initiate the coverage measurement for this transcript
                tCov = 0
                # If introns are to be considered as well, check coverage
                #  before extracting exons
                if includeIntrons:
                    trans0start = trans.start - 1       # Conversion to 0-based
                    trans0end = trans.end               # Conversion to 0-based
                    tCov += np.sum(
                        bam.count_coverage(contig = refname,
                                           start = trans0start,
                                           stop = trans0end,
                                           quality_threshold = 0,
                                           read_callback = lambda r:
                                               r.is_reverse != strd))
                    # Only transcripts with ANY coverage are considered
                    if tCov == 0:
                        continue
                    # Add this covered transcript in the ref-wide piece list
                    covPieces.append((trans0start, trans0end))
                # Extract exons from the transcript, last-to-first; 0-based
                exons = sorted([(e.start - 1, e.end) for e
                                in db.children(trans, featuretype = 'exon')],
                               reverse = strd)
                # If only exons are to be considered, check coverage now
                if not includeIntrons:
                    for start, end in exons:
                            tCov += np.sum(
                                bam.count_coverage(contig = refname,
                                                   start = start,
                                                   stop = end,
                                                   quality_threshold = 0,
                                                   read_callback = lambda r:
                                                       r.is_reverse != strd))
                    # Only transcripts with ANY exonic coverage are considered
                    if tCov == 0:
                        continue
                    trans0start = trans.start - 1       # Conversion to 0-based
                    trans0end = trans.end               # Conversion to 0-based
                    # Add exons of this exon-wise covered transcript to the
                    #  ref-wide piece list
                    covPieces.extend(exons)
                [geneID] = trans.attributes['gene_id']
                # Save this transcript to the list of covered transcripts
                covTransByStrdRef[strd][refname].append(
                    Transcript(geneID, trans0start, trans0end, exons))
            # Once complete for this strd/ref, flatten the ref-wide exons/trans
            #  pieces & get the non/coverage basline data
            covPieces = flattenIntervals(covPieces)
            for start, end in covPieces:
                pieceCov = np.sum(bam.count_coverage(contig = refname,
                                                     start = start,
                                                     stop = end,
                                                     quality_threshold = 0,
                                                     read_callback = lambda r:
                                                         r.is_reverse != strd),
                                  axis = 0)
                flatCov = np.count_nonzero(pieceCov)
                Pos += np.sum(pieceCov) if weightedCov else flatCov
                Neg += len(pieceCov) - flatCov
    
    bam.close()
    savePKL(out_transBaselineData, (covTransByStrdRef, Pos, Neg))
    
    return covTransByStrdRef, Pos, Neg
    

def getTransEnding(exons, endLength, strd):
    """Helper function to extract the exon-wise end pieces of a transcript,
    consisting of a list of exons, of a given total cumulative length.

    Parameters
    ----------
    exons : (list)
        The complete list of exons constituting a transcript: [ (start, end) ]
    endLength : (int)
        The total length of the exon-wise ending.
    strd : (bool)
        The strand on which the transcript is found.
        
    Returns
    -------
    transEndPieces : (list)
        [ (start, end) ]
    """
    
    # Initiate the flattened exon-wise end pieces for this transcript
    transEndPieces = []
    # Initiate how much is left to progress from the end
    remaining = int(endLength)
    # For each exon (ordered last-to-first wrt/ strand)
    for eStart, eEnd in sorted(exons, reverse = strd):
        # Get the exon length
        exLen = eEnd - eStart
        # If it is not longer than the remainder, add the entire exon &
        #  decrease the remainder
        if exLen < remaining:
            transEndPieces.append((eStart, eEnd))
            remaining -= exLen
        # Otherwise add only the remaining piece of the exon & stop
        else:
            transEndPieces.append((eEnd - remaining, eEnd) if strd
                                  else (eStart, eStart + remaining))
            break
    
    return transEndPieces


def getTransStartEnding(exons, endLength, strd):
    """Helper function to split a transcript represented by a list of exons
    into an exon-wise start and ending of a given cumulative ending length.
    
    Parameters
    ----------
    exons : (list)
        The complete list of exons constituting a transcript: [ (start, end) ]
    endLength : (int)
        The total length of the exon-wise ending to be escaped before starts
        are retained.
    strd : (bool)
        The strand on which the transcript is found.
        
    Returns
    -------
    transStartPieces : (list)
        [ (start, end) ]
    transEndPieces : (list)
        [ (start, end) ]
    """
    
    # Initialize the lists of pieces
    transStartPieces = []
    transEndPieces = []
    # Initiate how much is left to progress to escape the end
    remaining = int(endLength)
    # For each exon (ordered last-to-first wrt/ strand)
    for eStart, eEnd in sorted(exons, reverse = strd):
        # If end has not been escaped yet
        if remaining > 0:
            # Measure the size of the exon
            exLen = eEnd - eStart
            # If not longer than the remainder, add the entire exon to ends
            if exLen <= remaining:
                transEndPieces.append((eStart, eEnd))
            # Otherwise add the remainder to end pieces and the escaped part of
            #  the exon to starts
            else:
                transEndPieces.append((eEnd - remaining, eEnd) if strd
                                      else (eStart, eStart + remaining))
                transStartPieces.append((eStart, eEnd - remaining) if strd
                                        else (eStart + remaining, eEnd))
            # Either way, decrease the remainder by the exon length
            remaining -= exLen
        # Otherwise add the entire exon as is
        else:
            transStartPieces.append((eStart, eEnd))
    
    return transStartPieces, transEndPieces
    
    
def getTransEndSensSpec(endLength, bamfile, BLdata, includeIntrons,
                        getSensSpec = True, weightedCov = True):
    """Function to provide data for a ROC curve across various transcript end
    lengths.
    
    Parameters
    ----------
    endLength : (int)
        The length of exon-wise transcript ends to consider.
    bamfile : (str)
        Path to the (filtered) bam file.
    BLdata : (tuple)
        Baseline data needed for Sens & Spec calulations, such that
        ( { strd : { ref : [ Transcript ] } }, Pos, Neg )
    includeIntrons : (bool)
        Whether (positive or negative) coverage of introns should also be
        considered.
    getSensSpec : (bool)m optional
        Whether to calculate the Sens/Spec data. The default is True.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Alternatively, the sensitivity/specificity
        relationship is purely binary (flat) and the covered bases (1) have the
        same weight as non-covered ones (0). The default is True.

    Returns
    -------
    TP : (int)
        The true positive rate of RNAseq coverage in transcript starts (total
        coverage if weightedCov, otherwise only number of bps with coverage).
    FN : (int)
        The false negative rate of RNAseq coverage in transcript starts (total
        coverage if weightedCov, otherwise only number of bps with coverage).
        [Pos baseline for transStarts.]
    TN : (int)
        The true negative rate of RNAseq coverage of transcript starts (sum of
        bps with no coverage). [Neg baseline for transStarts.]
    FP : (int)
        The false positive rate of RNAseq coverage of transcript starts (sum
        of bps with no coverage).
    eachTransStartByStrdRef : (dict)
        { strd : { refName : [ ( (start, end), ... ) ] } } for start exon
        pieces for each covered transcript.
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
    # !!! Only consider transcripts with coverage (expressed)!
    #
    # In order to get accuracy, the following is needed: TP, FP, TN, FN
    # This is binarized wrt/ the expected coverage: >= / <
    # This information also can be used to calculate Sens and Spec
    
    logger.info(f'Checking transcript ends of length {endLength:,d}...')
    
    # Initialize the baseline data
    covTransByStrdRef, Pos, Neg = BLdata
    
    bam = pysam.AlignmentFile(bamfile, 'rb')

    # POSITIVE ~ Inside exon-wise transcript ends
    TP = 0  # Total coverage in exon-wise ends
    FP = 0  # Number of bp with 0 coverage in exon-wise ends
    # NEGATIVE ~ Outside exon-wise transcript ends
    #  FN ~ Total coverage in outside exon-wise ends
    #  TN ~ Number of bp with 0 coverage outside exon-wise ends
    
    # Initialize a dictionary to save transcript starts for later use:
    eachTransStartByStrdRef = defaultdict(lambda: defaultdict(list))
    # Go over each strand and reference separately
    for strd, covTransByRef in covTransByStrdRef.items():
        for refname, trans in covTransByRef.items():
            # Initiate the flattened (exon-wise) transcript start/end pieces
            #  for this reference
            eachTransStart = []
            allTransEndPieces = []
            # Iterate through covered transcripts to extract start & end pieces
            #  of given length for each
            for Tran in trans:
                if includeIntrons:
                    # If the transcipt is no longer than end length, add the
                    #  entire transcript to ends and add no starts
                    if Tran.end - Tran.start <= endLength:
                        transEndPieces = [(Tran.start, Tran.end)]
                        transStartPieces = []
                    # Otherwise add both exon-wise and absolute transcript ends
                    else:
                        transEndPieces = getTransEnding(Tran.exons, endLength,
                                                        strd)                    
                        transEndPieces.append(
                            (Tran.end - endLength, Tran.end) if strd
                            else (Tran.start, Tran.start + endLength))
                        transEndPieces = flattenIntervals(transEndPieces)
                        transStartPieces = removeOverlaps(
                            [(Tran.start, Tran.end)], transEndPieces)
                else:
                    transStartPieces, transEndPieces = getTransStartEnding(
                        Tran.exons, endLength, strd)
                # Add the flattened ends and starts (if any) to the resp lists
                allTransEndPieces.extend(transEndPieces)
                if transStartPieces != []:
                    eachTransStart.append(tuple(transStartPieces))
            # Flatten the allTransEndPieces
            allTransEndPieces = flattenIntervals(allTransEndPieces)
            # Remove any overlaps of the starts with the ends
            for transStart in eachTransStart:
                eachTransStartByStrdRef[strd][refname].append(
                    tuple(removeOverlaps(list(transStart), allTransEndPieces)))
            # For each reference (on each strand), extract and add the exon
            #  coverage and non-covBP for the end pieces
            if getSensSpec:
                for pStart, pEnd in allTransEndPieces:
                    endsCov = np.sum(
                        bam.count_coverage(contig = refname,
                                           start = pStart,
                                           stop = pEnd,
                                           quality_threshold = 0,
                                           read_callback = lambda r:
                                               r.is_reverse != strd),
                            axis = 0)
                    flatCov = np.count_nonzero(endsCov)
                    TP += np.sum(endsCov) if weightedCov else flatCov
                    FP += len(endsCov) - flatCov
    
    bam.close()
    # Calculate the results
    TN = Neg - FP
    FN = Pos - TP
    sensitivity = TP / Pos
    specificity = TN / Neg
    if not 0 <= sensitivity <= 1:
        raise ValueError(f'Sensitivity is {sensitivity}.')
    if not 0 <= specificity <= 1:
        raise ValueError(f'Specificity is {specificity}.')
    
    logger.info(f'Processed transcript ends of length {endLength:,d}.')
    
    return TP, FN, TN, FP, eachTransStartByStrdRef


def keyVal(key, d, meth):
    """A helper function for finding a maximum ROC point according to the
    specified optimizing function. In practice, used as a partial function with
    'key' as the missing arg.
    
    Parameters
    ----------
    key : (int)
        The key of the dictionary.
    d : (dict)
        The ROC dictionary whose maximum point is being sought.
    meth : (str)
        The name of the optimization method.

    Returns
    -------
    val : (float)
        The value assigned to this dict key based on the method used.
    """
    
    TP, FN, TN, FP = d[key][:4]
    
    # Methods from Berrar, 2019 (DOI: 10.1016/B978-0-12-809633-8.20351-8)
    if meth == 'accuracy':
        val = (TP+TN) / (TP+FP+FN+TN)
    
    elif meth == 'Youden':
        val = TP/(TP+FN) + TN/(TN+FP) - 1
        
    elif meth == 'LR+':
        val = TP/(TP+FN) / (1 - TN/(TN+FP))
    
    elif meth == 'LR-':
        val = (1 - TP/(TP+FN)) / (TN/(TN+FP))
        
    elif meth == 'BACC':
        val = (TP/(TP+FN) + TN/(TN+FP)) / 2
        
    elif meth == 'F':
        val = 2 * (1/(1/(TP/(TP+FP)) + 1/(TP/(TP+FN))))
        
    elif meth == 'G':
        val = np.sqrt(TP/(TP+FP) * TP/(TP+FN))
        
    elif meth == 'MCC':
        val = ((TP*TN - FP*FN) /
               np.sqrt(float((TP+FN) * (TP+FP) * (TN+FP) * (TN+FN))))
    
    # Methods from Unal, 2017 (DOI: 10.1155/2017/3762651)
    elif meth == 'product':
        val = TP/(TP+FN) * TN/(TN+FP)
    
    elif meth == 'ER':
        val = -np.sqrt((1-TP/(TP+FN))**2 + (1-TN/(TN+FP))**2)
    
    elif meth == 'IU':
        fpr = [1 - d[k][2]/(d[k][2] + d[k][3]) for k in sorted(d)]
        tpr = [d[k][0]/(d[k][0] + d[k][1]) for k in sorted(d)]
        AUC = metrics.auc(x = [0] + fpr + [1], y = [0] + tpr + [1])
        val = -(abs(TP/(TP+FN) - AUC) + abs(TN/(TN+FP) - AUC))
        
    # From https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/
    #  Procedures/NCSS/One_ROC_Curve_and_Cutoff_Analysis.pdf
    elif meth == 'DOR':
        val = (TP/(TP+FN))/(1 - TN/(TN+FP)) / ((1 - TP/(TP+FN))/(TN/(TN+FP)))
    
    # Simplified IU, independent of AUC
    elif meth == 'diag':
        val = -abs(TP/(TP+FN) - TN/(TN+FP))
    
    else:
        raise ValueError('The method has to be one of the following:\n'
                         '"accuracy", "Youden", "LR+", "LR-", "BACC", "F", '
                         '"G", "MCC", "product", "ER", "IU", "DOR", "diag"')
    if np.isnan(val):
        val = 0
    
    return val


def getTransEndROC(out_TransEndROC, out_transBaselineData, out_db, bamfile,
                   endLenLo, endLenHi, tROC = {}, optMeth = 'Youden',
                   includeIntrons = False, endLenMax = 1000,
                   weightedCov = True):
    """A gradient descent-like wrapper function around getTransEndSensSpec() to
    manage result generation & collection. This function maximizes the function
    given by optMeth across exon-wise transcript endLengths. This function
    assumes that there is only one such local maximum. (!) Also note that while
    endLenLo and endLenHi serve to initiate the search, they are non-binding
    and this algorithm may look outside these boundaries (up to endLenMax).

    Parameters
    ----------
    out_TransEndROC : (str)
        Path to the saved output, if done before.
    out_transBaselineData : (str)
        Path to the saved baseline data, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    bamfile : (str)
        Path to the (filtered) bam file.
    endLenLo : (int)
        The minimal estimated length of exon-wise transcript ends to consider.
    endLenHi : (int)
        The maximum estimated length of exon-wise transcript ends to consider.
    tROC : (dict), optional
        { endLen : ( Sensitivity, Specificity, Accuracy ) }. The default is {}.
    optMeth : (str)
        What function is used for optimization; the options are 'Youden',
        'MCC', or 'product'. The default is 'Youden'.
    includeIntrons : (bool), optional
        Indicates whether to consider intron coverage. The default is False.
    endLenMax : (int), optional
        Indicates the maximum endLen to consider
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Alternatively, the sensitivity/specificity
        relationship is purely binary (flat) and the covered bases (1) have the
        same weight as non-covered ones (0). The default is True.

    Returns
    -------
    tROC : (dict)
        { endLen : ( TP, FN, TN, FP, eachTransStartByStrdRef ) }
    """
    
    # If the file already exists, simply load if it contains the final optLen
    if os.path.isfile(out_TransEndROC):
        tROC = loadPKL(out_TransEndROC)
        optLen = max(tROC, key = partial(keyVal, d = tROC, meth = optMeth))
        if (optLen-1) in tROC and ((optLen+1) in tROC or optLen == endLenMax):
            logger.info(f'The file {out_TransEndROC} already exists with the '
                        f'optimal end length by {optMeth} being {optLen}.')
            return tROC
    
    # Get the baseline data
    BLdata = getBaselineData(out_transBaselineData, out_db, bamfile,
                             includeIntrons, weightedCov = weightedCov)
    
    # Start with min, mid, and max and then apply the repetitive algorithm
    #  until 'covergence'
    for endLen in (endLenLo,
                   np.mean((endLenLo, endLenHi), dtype = int),
                   endLenHi):
        if endLen not in tROC:
            tROC[endLen] = getTransEndSensSpec(
                endLen, bamfile, BLdata, includeIntrons,
                weightedCov = weightedCov)
        
    # Always look at the midpoint between the endLen with the largest
    #  Sens*Spec product and an adjacent endLen, alternating above or below.
    #  If looking above the longest or below the shortest endLen, look as far
    #  as an existing adjacent value in the opposite direction.
    
    # Initiate the alternating indicator of looking above or below
    lookAbove = False
    # Initiate the indicators of having looked at both optLen+/-1
    checkedJustBelowAbove = [False, False]
    # Initiate the oldOptLen at a value that can't be true before the 1st loop
    oldOptLen = -1
    # Run the while loop until both optLen+/-1 have been checked and neither
    #  is better than optLen ~ the optimum has been found with max resolution
    while not all(checkedJustBelowAbove):
        checkedLens = sorted(tROC)
        # Get the current most optimal endLen using the J statistic or product
        optLen = max(tROC, key = partial(keyVal, d = tROC, meth = optMeth))
        logger.info('The current optimal end length by '
                    f'{optMeth} is {optLen}.')
        # Settings in the special case when the endLenMax is reached
        if optLen == endLenMax:
            lookAbove = False
            checkedJustBelowAbove = [False, True]
            oldOptLen = optLen
        # If the new opt is different from the old, reset the indicators & save
        #  the new "old"; otherwise just flip the direction to look further
        elif optLen == oldOptLen:
            lookAbove = not lookAbove
        else:
            checkedJustBelowAbove = [False, False]
            oldOptLen = optLen
        # Get the index of the optimal optLen
        i = checkedLens.index(optLen)
        # Get the index of an adjacent endLen
        j = i + (-1, 1)[lookAbove]
        
        # If the adjacent endLen does not exist in the given direction, look
        #  as far as an existing adjacent value in the opposite direction
        # Also, limit the length to at most endLenMax
        if j not in range(len(checkedLens)):
            lenToC = min(
                optLen + (optLen - checkedLens[i + (-1, 1)[not lookAbove]]),
                endLenMax)
        # Otherwise just look between the optimal & adjacent endLen
        else:
            lenToC = np.mean((optLen, checkedLens[j]), dtype = int)
        
        # If this length has already been checked, this must be either optLen
        #  itself (if lookAbove = True, which means that optLen + 1 must have
        #  been checked as well, as that must be the checkedLens[j])
        #  or optLen-1 (if lookAbove = False), so don't check it again & only
        #  mark that an immediately adjacent value has already been checked
        if lenToC in checkedLens:
            logger.info(f'End length of {(lenToC, checkedLens[j])[lookAbove]} '
                        'has already been checked.')
            checkedJustBelowAbove[lookAbove] = True
        else:
            tROC[lenToC] = getTransEndSensSpec(
                lenToC, bamfile, BLdata, includeIntrons,
                weightedCov = weightedCov)
    
    logger.info(f'The final optimal end length is {optLen}.')
    # Saving the results; note that this may rewrite a previously saved file
    #  with additional data (no data should be lost in this step)
    savePKL(out_TransEndROC, tROC)
    
    return tROC


def sortSNRsByLenStrdRef(lenToSNRs):
    """Helper function to pre-sort SNRs by length, strand, and reference.
    
    Parameters
    ----------
    lenToSNRs : (dict)
        { SNR length : [ SNR ] }

    Returns
    -------
    SNRsByLenStrdRef : (dict)
        { SNR length : { strand : { reference : [ SNR ] } } }
    """
    
    logger.info('Sorting SNRs by length, strand, and reference...')
    SNRsByLenStrdRef = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list)))
    for length, SNRs in lenToSNRs.items():
        for SNR in SNRs:
            SNRsByLenStrdRef[length][SNR.strand][SNR.record].append(SNR)
        
    return SNRsByLenStrdRef


def getSNREndROC(SNRsByLenStrdRef, tROC, out_SNREndROC, bamfile,
                 optMeth = 'Youden', sortedSNRs = True, weightedCov = True):
    """This alhorithm goes over ALL SNR lengths exactly once to determine the
    sensitivity and specificity of coverage captured by the SNR pieces. This is
    different from the transcript end algorithm, where there are many more
    possibilities and info on one end length is relatively independent of the
    others. In case of SNR lengths, we are interested in aggregating info about
    all SNRs above some threshold, so it makes sense to scan all of them in
    decreasing order exactly once. As an added benefit, the function *could*
    also coverge/terminate as soon as J statistic starts decreasing, as only
    one local maximum is assumed. (Not implemented.)
    
    Parameters
    ----------
    SNRsByLenStrdRef : (dict)
        { length : { strd : { ref : [ SNRs ] } } } or, alternatively when
        sortedSNRs = False: { length : [ SNRs ] }
    tROC : (dict)
        { endLen : ( TP, FN, TN, FP, eachTransStartByStrdRef ) }
    out_SNREndROC : (str)
        Location where the output of this function should be stored.
    bamfile : (str)
        Location of the bamfile to be scanned.
    optMeth : (str)
        What function is used for optimization; the options are 'Youden',
        'MCC', or 'product'. The default is 'Youden'.
    sortedSNRs : (bool)
        Indicates whether the SNRs are already sorted by length, strd, and ref.
        The default is True.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Alternatively, the sensitivity/specificity
        relationship is purely binary (flat) and the covered bases (1) have the
        same weight as non-covered ones (0). The default is True.

    Returns
    -------
    snrROC : (dict)
        { SNR length : ( TP, FN, TN, FP ) }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_SNREndROC):
        snrROC = loadPKL(out_SNREndROC)
        return snrROC
    
    # Sort SNRs by len, strd & ref     
    if not sortedSNRs:
        SNRsByLenStrdRef = sortSNRsByLenStrdRef(SNRsByLenStrdRef)
    
    snrROC = {}

    # Get the optimal length using the J statistic (or product) & announce
    optEndLen = max(tROC, key = partial(keyVal, d = tROC, meth = optMeth))
    logger.info(f'The optimal coverage distance length is {optEndLen}. '
                'Getting an ROC across SNR lengths...')
    # Extract the transROC data for this length
    eachTransStartByStrdRef = tROC[optEndLen][4]
    # Flatten to get transStartsByStrdRef
    transStartsByStrdRef = defaultdict(lambda: defaultdict(list))
    for strd, eachTransStartByRef in eachTransStartByStrdRef.items():
        for refName, eachTransStarts in eachTransStartByRef.items():
            # First, unpack the tuple of tuples; then flatten at once
            transStarts = []
            for eachTransStart in eachTransStarts:
                transStarts.extend(eachTransStart)
            transStartsByStrdRef[strd][refName] = flattenIntervals(transStarts)
    
    # The SNR length cutoff will describe the minimal SNR length included in
    #  non-canonical coverage. For SNRs at each length, flatten the SNRends and
    #  then find where they overlap with the transStarts & get the coverage
    SNRpiecesByStrdRef = defaultdict(lambda: defaultdict(list))
    Pos = tROC[optEndLen][1]  # FN from the transEnd coverage measurement
    Neg = tROC[optEndLen][2]  # TN from the transEnd coverage measurement
    # Initiate the values to be adjusted for each SNR length
    # "POSITIVE" ~ Inside exon-wise SNR pieces; T ~ coverage; F ~ zero cov bps
    TP = 0  # Total coveragee in exon-wise pieces
    FP = 0  # BPs with 0 coverage in exon-wise pieces

    bam = pysam.AlignmentFile(bamfile, 'rb')
    # Start with the longest length whose data will be included in all
    #  subsequent calculations. Always extract only the newly added pieces and
    #  add the coverage data to the already existing for speed.
    for length in sorted(SNRsByLenStrdRef, reverse = True):
        logger.info(f'Checking coverage for SNRs of length {length}+...')
        # Work by strand and reference
        for strd, transStartsByRef in transStartsByStrdRef.items():
            for refName, transStarts in transStartsByRef.items():
                refLen = bam.get_reference_length(refName)
                newSNRpieces = []
                # Get & flatten the SNRs' pieces at this length/strd/ref
                if strd:
                    for SNR in SNRsByLenStrdRef[length][strd][refName]:
                        start = max(0, SNR.start - optEndLen)
                        end = SNR.start
                        newSNRpieces.append((start, end))
                else:
                    for SNR in SNRsByLenStrdRef[length][strd][refName]:
                        start = SNR.end
                        end = min(refLen, SNR.end + optEndLen)
                        newSNRpieces.append((start, end))
                # Flatten the newSNRpieces
                newSNRpieces = flattenIntervals(newSNRpieces)
                # Retain only where these overlap with the transStarts
                newSNRpieces = getOverlaps(newSNRpieces, transStarts)
                # Remove already covered areas so only new ones are measured
                newSNRpieces = removeOverlaps(
                    newSNRpieces, SNRpiecesByStrdRef[strd][refName])
                # Measure the coverage of the new pieces only
                for pStart, pEnd in newSNRpieces:
                    pieceCov = np.sum(
                        bam.count_coverage(contig = refName,
                                           start = pStart,   # already 0-based
                                           stop = pEnd,
                                           quality_threshold = 0,
                                           read_callback = lambda r:
                                               r.is_reverse != strd),
                        axis = 0)
                    # Adjust the TP & TN accordingly
                    flatCov = np.count_nonzero(pieceCov)
                    TP += np.sum(pieceCov) if weightedCov else flatCov
                    FP += len(pieceCov) - flatCov
                # Add the newSNRpieces to those from previous SNR lengths
                SNRpiecesByStrdRef[strd][refName].extend(newSNRpieces)
                # Flatten the SNRs (should contain no overlaps, only adjacency)
                SNRpiecesByStrdRef[strd][refName] = flattenIntervals(
                    SNRpiecesByStrdRef[strd][refName])
        # Save the results
        TN = Neg - FP
        FN = Pos - TP
        sensitivity = TP / Pos
        specificity = TN / Neg        
        
        if not 0 <= sensitivity <= 1:
            raise ValueError(f'Sensitivity is {sensitivity}.')
        if not 0 <= specificity <= 1:
            raise ValueError(f'Specificity is {specificity}.')

        snrROC[length] = TP, FN, TN, FP
    
    bam.close()
    # Announce the optimal SNR length
    logger.info('The optimal SNR length is {}.'.format(
        max(snrROC, key = partial(keyVal, d = snrROC, meth = optMeth))))
    savePKL(out_SNREndROC, snrROC)

    return snrROC


def getSNRcovByTrans(SNRsByLenStrdRef, tROC, out_snrROC, bamfile,
                     optMeth = 'Youden', sortedSNRs = True,
                     weightedCov = True):
    """This function gets Sensitivity & Specificity of coverage accounted for
    by SNRs in expressed transcript starts (defined by the tROC), in which they
    occur. This is an alternative to getSNREndROC().
    
    Parameters
    ----------
    SNRsByLenStrdRef : (dict)
        { length : { strd : { ref : [ SNRs ] } } }
    tROC : (dict)
        { endLen : ( Sens, Spec, eachTransStartByStrdRef, TP, TN ) }
    out_snrROC : (str)
        Location of this function's saved output.
    bamfile : (str)
        Location of the bamfile to be scanned.
    optMeth : (str)
        What function is used for optimization; the options are 'Youden',
        'MCC', or 'product'. The default is 'Youden'.
    sortedSNRs : (bool)
        Indicates whether the SNRs are already sorted by length, strd, and ref.
        The default is True.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Alternatively, the sensitivity/specificity
        relationship is purely binary (flat) and the covered bases (1) have the
        same weight as non-covered ones (0). The default is True.
        
    Returns
    -------
    snrROC : (dict)
        { SNR length : ( Sensitivity, Specificity ) }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_snrROC):
        snrROC = loadPKL(out_snrROC)
        return snrROC
    
    # Sort SNRs by len, strd & ref     
    if not sortedSNRs:
        SNRsByLenStrdRef = sortSNRsByLenStrdRef(SNRsByLenStrdRef)
    
    # Derive the best transcript end length, which is also the SNR coverage len
    covLen = max(tROC, key = partial(keyVal, d = tROC, meth = optMeth))
    logger.info(f'The optimal coverage distance length is {covLen}.'
                'Getting ROC across SNR lengths...')
    # Extract the relevant transcripts starts
    eachTransStartByStrdRef = tROC[covLen][2]
    
    # Initialize the SNRpieces found on covered exons by strd & ref - coverage
    #  of these will consitute TP & FP measurements
    SNRpiecesByStrdRef = defaultdict(lambda: defaultdict(list))
    # Initialize the associated measurements that'll be added to over SNR lens
    TP = 0
    FP = 0
    # Initialize the running dict of flattened transStarts whose transcripts
    #  have been found to contain SNRpieces - coverage of these will constitute
    #  the baseline/denominator for Sens/Spec (Pos, Neg)
    flatTransStartsByStrdRef = defaultdict(lambda: defaultdict(list))
    # Initialize the associated measurements that'll be added to over SNR lens
    Pos = 0
    Neg = 0
    # Initialize the ROC results dict
    snrROC = {}
    # Load the bamfile
    bam = pysam.AlignmentFile(bamfile, 'rb')
    # Go over SNRs by length, from longest to shortest
    for length in sorted(SNRsByLenStrdRef, reverse = True):
        logger.info('Getting Sensitivity & Specificity for SNRs of length'
                    f' {length}...')
        # For each length, go over each strd & ref of trans starts
        for strd in eachTransStartByStrdRef.keys():
            for refName in eachTransStartByStrdRef[strd].keys():
                refLen = bam.get_reference_length(refName)
                # Initialize the list of SNR pieces *specific* for this SNR len
                newSNRpieces = []
                # Extract the SNR pieces for this strd/ref
                for SNR in SNRsByLenStrdRef[length][strd][refName]:
                    if strd:
                        start = max(0, SNR.start - covLen)
                        end = SNR.start
                    else:
                        start = SNR.end
                        end = min(refLen, SNR.end + covLen)
                    newSNRpieces.append((start, end))
                # Determine which portions of the flattened SNR pieces are NEW
                #  for this strd/ref, as only those will contribute to
                #  - adding new transcript starts (new Pos/Neg info) and
                #  - new TP/FP coverage information
                newSNRpieces = removeOverlaps(
                    flattenIntervals(newSNRpieces),
                    SNRpiecesByStrdRef[strd][refName])
                # These SNR pieces, specific for this SNR length, can overlap
                #  either the already included (flat list) expressed tStarts,
                #  or some of the new ones (each list), or both (!)
                
                # Retain overlaps with previously added (flat list) tStarts:
                #  positive selection of SNR pieces based on old transStarts
                ovSNRpieces = getOverlaps(
                    newSNRpieces, flatTransStartsByStrdRef[strd][refName])
                
                # Initialize the list of tStarts specific for this SNR len
                newTstarts = []
                # Initialize the set of tStarts that have already been added to
                #  the flattened list -> to be removed from "discovery"
                toRemove = set()
                # For *each* transcript start that had not been added yet:
                #  positive selection of new transStarts and SNR pieces based
                #  on mutual overlap (note: each tran is a tuple of pieces)
                for tran in eachTransStartByStrdRef[strd][refName]:
                    # Get the overlaps with SNR pieces specific for this length
                    SNRsInNewTransStart = getOverlaps(newSNRpieces, list(tran))
                    # If there are any,
                    if SNRsInNewTransStart != []:
                        # Add the overlap to the overlapped SNR pieces
                        ovSNRpieces.extend(SNRsInNewTransStart)
                        # Add the transcript start exons to the tStarts aggr
                        newTstarts.extend(tran)
                        # Add the transcript start to the removal list
                        toRemove.add(tran)
                # Remove the transcript starts from the original dict to
                #  increase efficiency (not to try to "rediscover" them again
                #  in the future, since they are already added)
                eachTransStartByStrdRef[strd][refName] = [
                    tS for tS in eachTransStartByStrdRef[strd][refName]
                    if tS not in toRemove]
                
                # Flatten the ovSNRpieces
                ovSNRpieces = flattenIntervals(ovSNRpieces)
                # Measure the contribution of the ovSNRpieces specific for this
                #  SNR length to the total TP/FP
                for ovStart, ovEnd in ovSNRpieces:
                    pieceCov = np.sum(
                        bam.count_coverage(contig = refName,
                                           start = ovStart,   # already 0-based
                                           stop = ovEnd,
                                           quality_threshold = 0,
                                           read_callback = lambda r:
                                               r.is_reverse != strd),
                        axis = 0)
                    # Adjust the TP & FP accordingly
                    flatCov = np.count_nonzero(pieceCov)
                    TP += np.sum(pieceCov) if weightedCov else flatCov
                    FP += len(pieceCov) - flatCov
                # Itegrate the ovSNRpieces into the aggregate dict & flatten
                SNRpiecesByStrdRef[strd][refName].extend(ovSNRpieces)
                SNRpiecesByStrdRef[strd][refName] = flattenIntervals(
                    SNRpiecesByStrdRef[strd][refName])
                    
                # Flatten the newTstarts
                newTstarts = flattenIntervals(newTstarts)
                # Remove the overlaps with the previously added tStarts
                newTstarts = removeOverlaps(
                    newTstarts, flatTransStartsByStrdRef[strd][refName])
                # Measure how much only the newTstarts, specific to this SNR
                #  length, contribute to the total Pos/Neg
                for tStart, tEnd in newTstarts:
                    tCov = np.sum(
                        bam.count_coverage(contig = refName,
                                           start = tStart,   # already 0-based
                                           stop = tEnd,
                                           quality_threshold = 0,
                                           read_callback = lambda r:
                                               r.is_reverse != strd),
                        axis = 0)
                    # Adjust the Pos & Neg accordingly
                    flatCov = np.count_nonzero(tCov)
                    Pos += np.sum(tCov) if weightedCov else flatCov
                    Neg += len(tCov) - flatCov
                # Itegrate the newTstarts into the aggregate dict & flatten
                flatTransStartsByStrdRef[strd][refName].extend(newTstarts)
                flatTransStartsByStrdRef[strd][refName] = flattenIntervals(
                    flatTransStartsByStrdRef[strd][refName])
                    
        # Calculate the Sens & Spec at this minimal SNR length and save
        if Pos != 0 and Neg != 0:
            TN = Neg - FP
            FN = Pos - TP
            Sens = TP / Pos
            Spec = TN / Neg
            if not 0 <= Sens <= 1:
                raise ValueError("Sensitivity is {}.".format(Sens))
            if not 0 <= Spec <= 1:
                raise ValueError("Specificity is {}.".format(Spec))
            snrROC[length] = TP, FN, TN, FP
    
    bam.close()
    savePKL(out_snrROC, snrROC)
    
    return snrROC


def getStatsByGene(covLen, minSNRlen, lenToSNRs, out_geneStats, out_db,
                   out_transBaselineData, bamfile, includeIntrons = False,
                   weightedCov = True):
    """Function to get several statistics in order to calculate correlation
    between gene coverage & SNR content. Specifically, it goes over all genes
    and for each gene, it measures its RNA-seq coverage. If > 0, it also gets
    gene length, starts coverage, # of SNRs (of given minimal length) and the
    associated total potential coverage area (that upstream of these SNRs
    overlapping exons). This data can be used to plot total coverage / starts
    coverage vs. # of SNRs / SNR area for each gene, normalized by gene length,
    if needed.

    Parameters
    ----------    
    covLen : (int)
        The optimal length of coverage with respect to the potential priming
        location.
    minSNRlen : (int)
        Minimal length of SNRs to be included.
    lenToSNRs : (dict)
        { SNR length : [ SNR ] }
    out_geneStats : (str)
        Path to the saved output of this function.
    out_db : (str)
        Path to the saved GTF/GFF database.
    out_transBaselineData : (str)
        Path to the saved baseline data.
    bamfile : (str)
        Location of the bamfile to be scanned.
    includeIntrons : (bool), optional
        Indicates whether to consider intron coverage. The default is False.
    weightedCov : (bool), optional
        Determines whether the covered bases are weighted by coverage amount
        (i.e., coverage is summed). Alternatively, the sensitivity/specificity
        relationship is purely binary (flat) and the covered bases (1) have the
        same weight as non-covered ones (0). The default is True.

    Returns
    -------
    statsByGene : (dict)
        { gene : (exon coverage, exon length,
                  starts coverage, starts length,
                  SNR exon area, SNR exon number,
                  SNR start area, SNR start number) }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_geneStats):
        statsByGene = loadPKL(out_geneStats)
        return statsByGene
    
    # Get a dictionary of SNRs by gene & length
    SNRsByGeneLen = getSNRsByGeneLen(lenToSNRs, concordant = True)
    # Get the list of covered transcripts by strd & ref
    BLdata = getBaselineData(out_transBaselineData, out_db, bamfile,
                             includeIntrons, weightedCov = weightedCov)
    covTransByStrdRef = BLdata[0]
    # Sort the covered transcripts by geneID
    covTransByGene = defaultdict(list)
    for covTransByRef in covTransByStrdRef.values():
        for covTrans in covTransByRef.values():
            for covTran in covTrans:
                covTransByGene[covTran.geneID].append(covTran)
    # Initialize a dictionary for the results
    statsByGene = {}
    # Connect the gff database
    db = gffutils.FeatureDB(out_db, keep_order = True)
    # Count the total number of genes
    totalGenes = len(covTransByGene)
    # Load the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    # Initialize the progress tracker variables
    progressTracker = 0
    prog = 0
    # Go over all genes
    for gene in db.features_of_type(featuretype = 'gene'):
        geneID = gene.id
        if geneID not in covTransByGene:
            continue
        progressTracker += 1
        newProg = round(progressTracker / totalGenes, 3)
        if newProg != prog:
            prog = newProg
            logger.info(f'Getting gene coverage statistics... ({prog:.1%})')
            clear_output(wait = True)
        strd = gene.strand == '+'
        refName = gene.seqid
        refLen = bam.get_reference_length(refName)
        # Extract all exons from this gene's covered transcripts
        geneExons = [(eS, eE) for Trans in covTransByGene[geneID]
                     for eS, eE in Trans.exons]
        # Flatten these exons into non-overlapping intervals
        geneExons = flattenIntervals(geneExons)
        # Measure the total exon-wise length & RNA-seq coverage of this gene
        eLen = 0
        eCov = 0
        for eStart, eEnd in geneExons:
            eLen += eEnd - eStart
            eCov += np.sum(bam.count_coverage(contig = refName,
                                              start = eStart,
                                              stop = eEnd,
                                              quality_threshold = 0,
                                              read_callback = lambda r:
                                                  r.is_reverse != strd))
        # Process this gene further only if it is expressed
        if eCov != 0:
            # Remove the exon-wise transcript endings for this gene by going
            #  over each transcript to get its ending, then merging & removing
            # Initiate the list of endings
            transEnds = []
            for Trans in covTransByGene[geneID]:
                # Extract the endings for this transcript and add
                transEnds.extend(getTransEnding(Trans.exons, covLen, strd))
            # Flatten the endings
            transEnds = flattenIntervals(transEnds)
            # Remove these endings from the gene exons
            geneStartsOnly = removeOverlaps(geneExons, transEnds)
            # Get the starts length & coverage
            sLen = 0
            sCov = 0
            for sStart, sEnd in geneStartsOnly:
                sLen += sEnd - sStart
                sCov += np.sum(bam.count_coverage(contig = refName,
                                                  start = sStart,
                                                  stop = sEnd,
                                                  quality_threshold = 0,
                                                  read_callback = lambda r:
                                                      r.is_reverse != strd))
            # Count & get the exon-overlapping pieces for all SNRs of given
            #  minimal length overlapping this gene; collect the data for all
            #  exons of the gene, or for starts only
            SNRnum = 0
            SNRpieces = []
            SNRstartNum = 0
            SNRstartPieces = []
            for length, SNRs in SNRsByGeneLen[geneID].items():
                if length >= minSNRlen:
                    for SNR in SNRs:
                        if strd:
                            start = max(0, SNR.start - covLen)
                            end = SNR.start
                        else:
                            start = SNR.end
                            end = min(refLen, SNR.end + covLen)
                        # Count this SNR if its piece overlaps with gene exons
                        exonPieceOverlaps = getOverlaps([(start, end)],
                                                        geneExons)
                        if exonPieceOverlaps != []:
                            SNRpieces.extend(exonPieceOverlaps)
                            SNRnum += 1
                        # See if this SNR's piece overlaps with starts only
                        startPieceOverlaps = getOverlaps([(start, end)],
                                                         geneStartsOnly)
                        if startPieceOverlaps != []:
                            SNRstartPieces.extend(startPieceOverlaps)
                            SNRstartNum += 1
            # Flatten the SNRpieces
            SNRpieces = flattenIntervals(SNRpieces)
            # Calculate the total length of overlap between SNR pieces & exons
            SNRlen = sum(pEnd - pStart for pStart, pEnd in SNRpieces)
            # Flatten the SNR start pieces
            SNRstartPieces = flattenIntervals(SNRstartPieces)
            # Calculate the total length of overlap between SNR pieces & starts
            SNRstartLen = sum(pEnd - pStart for pStart, pEnd in SNRstartPieces)
            # Save the results for this gene in the following order:
            #  (exon coverage, exon length, starts coverage, SNR area, SNR #)
            statsByGene[geneID] = (eCov, eLen, sCov, sLen, SNRlen, SNRnum,
                                    SNRstartLen, SNRstartNum)
    bam.close()
    savePKL(out_geneStats, statsByGene)
    
    return statsByGene
