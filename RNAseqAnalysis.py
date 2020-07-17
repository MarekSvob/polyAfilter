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
from IPython.display import clear_output
from SNRanalysis import getFlatFeatsByTypeStrdRef, getSNRsByGeneLen, \
    flattenIntervals, getOverlaps, removeOverlaps
from SNRdetection import savePKL, loadPKL

expCovByConcFeat = collections.defaultdict(lambda:collections.defaultdict(int))

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
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db,
                                          (featType, ))
    
    # Connect to the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # Initiate the varuables needed
    # This needs to start off as a set to avoid double-adding reads aligned to
    #  two features of the same type
    filt_reads = set()
    
    # Get only reads on the same/other strand (depends if concordant or not)
    #  as the feature in question
    for strd, featsByStrdRef in flatFeats[featType].items():
        # Cycle over the feats under the appropriate key
        for ref, feats in featsByStrdRef.items():
            print('Processing {} strand of reference {}'.format(
                '+' if strd else '-', ref))
            for start, end in feats:
                for read in bam.fetch(contig = ref, start = start, end = end):
                    # Make sure the read is on the correct strand before adding:
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
    
    print('A filtered bam file has been created.')
    print('!!! Before continuing, index it in terminal using command '\
          '"samtools index {}"'.format(filt_bamfile))


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
    
    global expCovByConcFeat
    
    # Try to retrieve the value if already calculated in this session
    if expCovByConcFeat[concordant][baselineFeat] != 0:
        expected = expCovByConcFeat[concordant][baselineFeat]
        return expected
    
    # In most cases, this should just load a file
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db,
                                          (baselineFeat, ))
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')

    totalCoverage = 0
    totalLength = 0

    # Get the total length of the feats and the total coverage over them
    for strd, featsByRef in flatFeats[baselineFeat].items():
        print('Scanning {}{}s for total coverage and length...'.format(
            '+' if strd else '-', baselineFeat))
        for ref, feats in featsByRef.items():
            for start, end in feats:
                totalLength += end - start
                totalCoverage += np.sum(
                    bam.count_coverage(
                        contig = ref, start = start, stop = end,
                        quality_threshold = 0,
                        read_callback = lambda r: r.is_reverse != strd),
                    dtype = 'int')
    bam.close()
    
    # Normalize the coverage by #SNRs, total coverage, and total length
    expected = totalCoverage / totalLength
    
    # Save the value for later use in this session
    expCovByConcFeat[concordant][baselineFeat] = expected
    
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
    outlierPeaks = collections.defaultdict(list)
    totalCount = 0
    filteredOut = 0
    zeroCovs = 0
    # Go over the SNRs by length, starting with the longest
    for length in sorted(lenToSNRs.keys(), reverse = True):
        if length in SNRlengths:
            # For each length, initialize the count & array for window coverage
            print('Going over {} coverage of SNR{}s...'.format(
                'concordant' if concordant else 'discordant', length))
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
    print('Filtered out {:.2%} outliers and {:.2%} SNRs had no '\
          'coverage.'.format(filteredOut / totalCount, zeroCovs / totalCount))
    
    return covByLen, outlierPeaks


def getCovPerTran(out_TranCov, out_db, out_strandedFeats, exonic_bamfile,
                  window = 2000):
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
    expCov = getExpectedCoverage(out_db, out_strandedFeats, exonic_bamfile,
                                 concordant = True,
                                 baselineFeat = 'transcript')
    
    # In most cases, this should just load a file
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db,
                                          ('transcript', ))
    
    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(exonic_bamfile, 'rb')
    
    # Initiate the variables needed
    transCount = 0
    coverage = np.zeros((window), 'L')
    
    for strd, featsByRef in flatFeats['transcript'].items():
        print('Going over coverage of {}transcripts'.format(
            '+' if strd else '-'))
        for ref, feats in featsByRef.items():
            refLen = bam.get_reference_length(ref)
            for feat in feats:
                if strd:
                    mid = feat[1]
                else:
                    mid = feat[0]
                start = int(mid - window / 2)
                stop = int(mid + window / 2)
                # Include correctins for the start & end if the window
                #  falls out of the reference size range
                if start < 0:
                    corrStart = 0
                else:
                    corrStart = start
                if stop > refLen:
                    corrStop = refLen
                else:
                    corrStop = stop
                # Get the coverage summed over A/T/C/G; count only concordant
                refCoverage = np.sum(
                    bam.count_coverage(contig = ref,
                                       start = corrStart,
                                       stop = corrStop,
                                       quality_threshold = 0,
                                       read_callback = lambda r:
                                           r.is_reverse != strd),
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
                # Add the SNR
                coverage += refCoverage
    # Normalize the coverage only once at the end
    normCov = coverage / (transCount * expCov)
    
    bam.close()
    savePKL(out_TranCov, normCov)
    
    return normCov


def getNonCanCovGenes(out_NonCanCovGenes, lenToSNRs, out_db, bamfile,
                      lastBP = 250, covThreshold = 0.05, minLen = 0):
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
    geneLenToSNRs = getSNRsByGeneLen(lenToSNRs, concordant = True)
    db = gffutils.FeatureDB(out_db, keep_order = True)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    progressTracker = 0
    prog = 0
    # For each gene that has an SNR detected
    for geneID, glenToSNRs in geneLenToSNRs.items():
        progressTracker += 1
        newProg = round(progressTracker / len(geneLenToSNRs), 3)
        if newProg != prog:
            prog = newProg
            print('Looking for genes with non-canonical'\
                  ' coverage... ({:.1%})'.format(prog))
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
        #  add up the coverage of exon-wise last X bp
        for trans in db.features_of_type(
                featuretype = 'transcript',
                limit = (geneFeat.seqid, gene0start, gene0end),
                strand = geneFeat.strand):
            # Extract the exons as (start, end); 1- => 0-based
            exons = sorted([(e.start - 1, e.end) for e
                            in db.children(trans, featuretype = 'exon')],
                           reverse = strd) # Order from last to 1st wrt/ strand
            trans0start = trans.start - 1       # Conversion to 0-based
            trans0end = trans.end               # Conversion to 0-based
            [geneID] = trans.attributes['gene_id']
            # Save all transcripts (lite) on the same strand as the gene
            transcripts.append(Transcript(geneID, trans.seqid, trans0start,
                                          trans0end, exons))
            # Determine the transcripts's exon-wise end
            remaining = lastBP
            covStart = None
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
            #  beyond the gene end (in case it is from another gene)
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


def getBaselineData(out_transBaselineData, out_db, bamfile):
    """Function to obtain the baseline reference data (Transcripts & exons) as
    an input for calculation of sensitivity, specificity, and accuracy. Note
    that here, only expressed (i.e., with coverage) transcripts are considered.
    
    Parameters
    ----------
    out_transBaselineData : (str)
        Path to the saved output, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    bamfile : (str)
        Path to the (filtered) bam file.

    Returns
    -------
    covTransByStrdRef : (dict)
        { strd : { ref : [ Transcript ] } } of transcripts with non-0 exonic
        coverage.
    covExonsByStrdRef : (dict)
        { strd : { ref : [ exon ] } } flattened from transcripts above.
    Pos : (int)
        Total coverage across all flattened exons from covered transcripts.
    Neg : (int)
        Total number of exon bps with no coverage from covered transcripts.
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_transBaselineData):
        covTransByStrdRef, covExonsByStrdRef, Pos, Neg = loadPKL(
            out_transBaselineData)
        return covTransByStrdRef, covExonsByStrdRef, Pos, Neg
    
    # Otherwise initialize the variables    
    print('Getting baseline data for ROC across '\
          'exon-wise covered transcripts...')
    covTransByStrdRef = collections.defaultdict(
        lambda: collections.defaultdict(list)
        )
    covExonsByStrdRef = collections.defaultdict(
        lambda: collections.defaultdict(list)
        )
    Pos = 0
    Neg = 0
    
    db = gffutils.FeatureDB(out_db, keep_order = True)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    references = bam.references
    
    # Get all exon-covered transcripts and their flattened exons by strd & ref
    for strd in (True, False):
        for refname in references:
            # Process this transcript only if it has ANY exonic coverage
            for trans in db.features_of_type(
                    featuretype = 'transcript',
                    limit = (refname, 1, bam.get_reference_length(refname)),
                    strand = '+' if strd else '-'):
                # Extract exons from the transcript, last-to-first
                exons = sorted([(e.start - 1, e.end) for e
                                in db.children(trans, featuretype = 'exon')],
                               reverse = strd)
                eCov = 0
                for start, end in exons:
                        eCov += np.sum(
                            bam.count_coverage(contig = refname,
                                               start = start,
                                               stop = end,
                                               quality_threshold = 0,
                                               read_callback = lambda r:
                                                   r.is_reverse != strd))
                # Only transcripts with ANY exonic coverage are considered here
                if eCov == 0:
                    continue
                trans0start = trans.start - 1       # Conversion to 0-based
                trans0end = trans.end               # Conversion to 0-based
                [geneID] = trans.attributes['gene_id']
                # Save this transcript to the list of covered transcripts
                covTransByStrdRef[strd][refname].append(
                    Transcript(geneID, trans0start, trans0end, exons))
                # Add these exons in the ref-wide exon list
                covExonsByStrdRef[strd][refname].extend(exons)
            # Once complete for this strd/ref, flatten the ref-wide exon list
            covExonsByStrdRef[strd][refname] = flattenIntervals(
                covExonsByStrdRef[strd][refname])
            # Get the non/coverage of these exons to speed up future processing
            for start, end in covExonsByStrdRef[strd][refname]:
                exonicCov = np.sum(
                    bam.count_coverage(contig = refname,
                                       start = start,
                                       stop = end,
                                       quality_threshold = 0,
                                       read_callback = lambda r:
                                           r.is_reverse != strd),
                    axis = 0)
                Pos += np.sum(exonicCov)
                Neg += np.count_nonzero(exonicCov == 0)
    
    savePKL(out_transBaselineData,
            (covTransByStrdRef, covExonsByStrdRef, Pos, Neg))
    
    return covTransByStrdRef, covExonsByStrdRef, Pos, Neg
    

def getTransEnding(exons, endLength, strd):
    """Helper function to extract the exon-wise end pieces of a transcript,
    consisting of a list of exons, of a given total length.

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


def getTransStart(exons, endLength, strd):
    """Helper function to extract the exon-wise start pieces of a transcript,
    consisting of a list of exons, excluding an ending of a given total length.
    Note that this function could also be achieved by the combination of 
    getTransEnding() and then removeOverlaps() between the complete list of
    exons and transEndPieces, but this custom function should be a bit faster.
    
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
    """
    
    # Initialize the list of start pieces
    transStartPieces = []
    # Initiate how much is left to progress to escape the end
    remaining = int(endLength)
    # For each exon (ordered last-to-first wrt/ strand)
    for eStart, eEnd in sorted(exons, reverse = strd):
        # If end has not been escaped yet
        if remaining > 0:
            # Measure the size of the exon
            exLen = eEnd - eStart
            # If longer than the remainder, add the escaped part of the exon
            if exLen > remaining:
                transStartPieces.append((eStart, eEnd - remaining) if strd
                                        else (eStart + remaining, eEnd))
            # Either way, decrease the remainder by the exon length
            remaining -= exLen
        # Otherwise add the entire exon as is
        else:
            transStartPieces.append((eStart, eEnd))
    
    return transStartPieces
    
    
def getTransEndSensSpec(endLength, bamfile, BLdata):
    """Function to provide data for a ROC curve for various transcript end
    lengths.
    
    Parameters
    ----------
    endLength : (int)
        The length of exon-wise transcript ends to consider.
    bamfile : (str)
        Path to the (filtered) bam file.
    BLdata : (tuple)
        Baseline data needed for Sens & Spec calulations, such that
        ( { strd : { ref : [ Transcript ] } }, { strd : { ref : [ exon ] } },
         Pos, Neg )

    Returns
    -------
    sensitivity : (float)
        Proportion of total exonic coverage captured by the transcript ends of
        a given length: TP / (TP + FN)
    specificity : (float)
        Proportion of non-covered exon bps not overlapped by transcript ends:
        TN / (TN + FP)
    transStartsByStrdRef : (dict)
        { strd : { refName : [ (start, end) ] } } for exon pieces not covered
        by transcript ends.
    FN : (int)
        The total RNA-seq coverage in exon-wise transcript starts.
    TN : (int)
        The number of BPs with no coverage in exon-wise transcript starts.
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
    
    print('Checking transcript ends of length {}...'.format(endLength))
    
    # Initialize the baseline data
    covTransByStrdRef, covExonsByStrdRef, Pos, Neg = BLdata
    
    bam = pysam.AlignmentFile(bamfile, 'rb')

    # POSITIVE ~ Inside exon-wise transcript ends
    TP = 0  # Total coverage in exon-wise ends
    FP = 0  # Number of bp with 0 coverage in exon-wise ends
    # NEGATIVE ~ Outside exon-wise transcript ends
    #  FN ~ Total coverage in outside exon-wise ends
    #  TN ~ Number of bp with 0 coverage outside exon-wise ends
    
    # Initialize a dictionary to save refStarts:
    transStartsByStrdRef = collections.defaultdict(dict)
    # Go over each strand and reference separately
    for strd, covTransByRef in covTransByStrdRef.items():
        for refname, trans in covTransByRef.items():
            # Initiate the flattened exon-wise transcript end pieces for
            #  this reference
            transEndPieces = []
            # Iterate through previously found transcripts to extract &
            #  save exon-wise endpieces of given length for each
            for Tran in trans:
                transEndPieces.extend(getTransEnding(Tran.exons,
                                                     endLength,
                                                     strd))            
            # Flatten the transEndPieces
            transEndPieces = flattenIntervals(transEndPieces)
            # For each reference (on each strand), extract and add the exon
            #  coverage and non-covBP for the exon-wise end pieces
            for pStart, pEnd in transEndPieces:
                endsCov = np.sum(bam.count_coverage(contig = refname,
                                                    start = pStart,
                                                    stop = pEnd,
                                                    quality_threshold = 0,
                                                    read_callback = lambda r:
                                                        r.is_reverse != strd),
                                 axis = 0)
                TP += np.sum(endsCov)
                FP += np.count_nonzero(endsCov == 0)
            # Remove the portions of flattened exons overlapped by
            #  endpieces and save as starts
            transStartsByStrdRef[strd][refname] = removeOverlaps(
                covExonsByStrdRef[strd][refname],
                transEndPieces)
            
    bam.close()
    # Calculate the results
    FN = Pos - TP
    TN = Neg - FP
    sensitivity = TP / Pos
    specificity = TN / Neg
    
    return sensitivity, specificity, transStartsByStrdRef, FN, TN


def getTransEndROC(out_TransEndROC, out_transBaselineData, out_db, bamfile,
                   endLenMin, endLenMax, tROC = {}, product = False):
    """A gradient descent-like wrapper funciton around getTransEndSensSpec() to
    manage result generation & collection. This function maximizes Youden's J
    statistic (Sensitivity + Specificity - 1) across exon-wise transcript
    endLengths. This function assumes that there is only one such local
    maximum. (!) Also note that while endLenMin and endLenMax serve to initiate
    the search, they are non-binding and this algorithm may look outside these
    boundaries.

    Parameters
    ----------
    out_transROC : (str)
        Path to the saved output, if done before.
    out_transBaselineData : (str)
        Path to the saved baseline data, if done before.
    out_db : (str)
        Path to the saved GTF/GFF database.
    bamfile : (str)
        Path to the (filtered) bam file.
    endLenMin : (int)
        The minimal estimated length of exon-wise transcript ends to consider.
    endLenMax : (int)
        The maximum estimated length of exon-wise transcript ends to consider.
    tROC : (dict), optional
        { endLen : ( Sensitivity, Specificity, Accuracy ) }. The default is {}.
    product : (bool), optional
        Indicates whether the optimization is guided by maximizing the product
        of Sens*Spec; otherwise the J statistic (Sens+Spec-1) is used. The
        default is False.

    Returns
    -------
    tROC : (dict)
        { endLen : ( Sensitivity, Specificity, Accuracy ) }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_TransEndROC):
        tROC = loadPKL(out_TransEndROC)
        return tROC
    
    # Get the baseline data
    BLdata = getBaselineData(out_transBaselineData, out_db, bamfile)
    
    # Start with min, mid, and max and then apply the repetitive algorithm
    #  until 'covergence'
    for endLen in (endLenMin,
                   np.mean((endLenMin, endLenMax), dtype = int),
                   endLenMax):
        if endLen not in tROC:
            tROC[endLen] = getTransEndSensSpec(endLen, bamfile, BLdata)
        
    # Always look at the midpoint between the endLen with the largest
    #  Sens*Spec product and an adjacent endLen, alternating above or below.
    #  If looking above the longest or below the shortest endLen, look as far
    #  as an existing adjacent value in the opposite direction.
    
    # Initiate the alternating indicator of looking above or below
    lookAbove = False
    # Initiate the indicators of having looked at both optLen+/-1
    checkedJustAboveBelow = [False, False]
    # Initiate the oldOptLen at a value that can't be true before the 1st loop
    oldOptLen = -1
    # Run the while loop until both optLen+/-1 have been checked and neither
    #  is better than optLen ~ the optimum has been found with max resolution
    while not all(checkedJustAboveBelow):
        checkedLens = sorted(tROC)
        # Get the current most optimal endLen using the J statistic or product
        optLen = max(tROC, key = lambda l: (tROC[l][0] * tROC[l][1] if product
                                            else tROC[l][0] + tROC[l][1] - 1))
        print('The current optimal end length is {}.'.format(optLen))
        # If the new opt is different from the old, reset the indicators & save
        #  the new "old"; otherwise just flip the direction to look further
        if optLen == oldOptLen:
            lookAbove = not lookAbove 
        else:
            checkedJustAboveBelow = [False, False]
            oldOptLen = optLen
        # Get the index of the 
        i = checkedLens.index(optLen)
        # Get the index of the adjacent endLen
        j = i + (-1, 1)[lookAbove]
        
        # If the adjacent endLen does not exist in the given direction, look
        #  as far as an existing adjacent value in the opposite direction
        if j in (-1, len(checkedLens)):
            lenToC = optLen + (optLen
                               - checkedLens[i + (-1, 1)[not lookAbove]])
        # Otherwise just look between the optimal & adjacent endLen
        else:
            lenToC = np.mean((optLen, checkedLens[j]), dtype = int)
        
        # If this length has already been checked, this must be either optLen
        #  itself (if lookAbove = True, which means that optLen+1 must have
        #  been checked as well, as that must be the checkedLens[j])
        #  or optLen-1 (if lookAbove = False), so don't check it again & only
        #  mark that an immediately adjacent value has already been checked
        if lenToC in checkedLens:
            print('End length of {} has already been checked.'.format(
                (lenToC, checkedLens[j])[lookAbove]))
            checkedJustAboveBelow[lookAbove] = True
        else:
            tROC[lenToC] = getTransEndSensSpec(lenToC, bamfile, BLdata)
    
    print('The optimal end length is {}.'.format(optLen))
    savePKL(out_TransEndROC, tROC)
    
    return tROC


def getSNREndROC(tROC, lenToSNRs, out_SNREndROC, bamfile, product = False):
    """This alhorithm goes over ALL SNR lengths exactly once to determine the
    sensitivity and specificity of coverage captured by the SNR pieces. This is
    different from the transcript end algorithm, where there are many more
    possibilities and info on one end length is relatively independent of the
    others. In case of SNR lengths, we are interested in aggregating info about
    all SNRs above some threshold, so it makes sense to scan all of them in
    decreasing order exactly once. As an added benefit, the function *could*
    also coverge/terminate as soon as J statistic starts decreasing, as only
    one local maximum is assumed. (Not implemented.)
    Note that this function hass been replaced by getSNRcovByGene().
    
    Parameters
    ----------
    tROC : (dict)
        { endlength : (Sens, Spec, transStartsByStrdRef, FN, TN) }
    lenToSNRs : (dict)
        { SNRlength : [ SNR ] }
    out_SNREndROC : (str)
        Location where the output of this function should be stored.
    bamfile : (str)
        Location of the bamfile to be scanned.

    Returns
    -------
    snrROC : (dict)
        { covLength : (sensitivity, specificity) }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_SNREndROC):
        snrROC = loadPKL(out_SNREndROC)
        return snrROC
    
    snrROC = {}
    
    # Get the optimal length using the J statistic & announce
    optEndLen = max(tROC, key = lambda l: (tROC[l][0] * tROC[l][1] if product
                                           else tROC[l][0] + tROC[l][1] - 1))
    print('The optimal coverage distance length is {}. Getting ROC across ' \
          'SNR lengths...'.format(optEndLen))
    # Extract the transROC data for this length
    transStartsByStrdRef = tROC[optEndLen][2]
    # The SNR length cutoff will describe the minimal SNR length included in
    #  non-canonical coverage. For SNRs at each length, flatten the SNRends and
    #  then find where they overlap with the transStarts & get the coverage
    SNRpiecesByStrdRef = collections.defaultdict(
        lambda: collections.defaultdict(list))
    Pos = tROC[optEndLen][3]
    Neg = tROC[optEndLen][4]
    # "POSITIVE" ~ Inside exon-wise SNR pieces
    TP = 0  # Total coveragee in exon-wise pieces
    # "NEGATIVE" ~ Outside exon-wise SNR pieces
    TN = Neg  # BPs with 0 coverage outside exon-wise pieces
    
    bam = pysam.AlignmentFile(bamfile, 'rb')
    # Start with the longest length whose data will be included in all
    #  subsequent calculations. Always extract only the newly added pieces and
    #  add the coverage data (TP/FP/TN/FN) to the already existing for speed.
    for length in sorted(lenToSNRs, reverse = True):
        print('Checking coverage for SNRs of length {}+...'.format(length))
        # Work by strand and reference
        for strd, transStartsByRef in transStartsByStrdRef.items():
            for refName, transStarts in transStartsByRef.items():
                newSNRpieces = []
                # Get & flatten the SNRs' pieces at this length/strd/ref
                for SNR in lenToSNRs[length]:
                    # Filter the SNRs by strd and ref
                    if SNR.strand == strd and SNR.record == refName:
                        start = SNR.start - optEndLen if strd else SNR.end
                        end = SNR.start if strd else SNR.end + optEndLen
                        newSNRpieces.append((start, end))
                # Flatten the newSNRpieces
                newSNRpieces = flattenIntervals(newSNRpieces)
                # Retain only where these overlap with the transStarts (exons)
                newSNRpieces = getOverlaps(newSNRpieces, transStarts)
                # Remove already covered areas so only new ones are measured
                newSNRpieces = removeOverlaps(
                    newSNRpieces, SNRpiecesByStrdRef[strd][refName])                
                # Measure the coverage of the new pieces
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
                    TP += np.sum(pieceCov)
                    TN -= np.count_nonzero(pieceCov == 0)
                # Add the newSNRpieces to those from previous SNR lengths
                SNRpiecesByStrdRef[strd][refName].extend(newSNRpieces)
                # Flatten the SNRs (should contain no overlaps, only adjacency)
                SNRpiecesByStrdRef[strd][refName] = flattenIntervals(
                    SNRpiecesByStrdRef[strd][refName])
        # Save the results
        sensitivity = TP / Pos
        specificity = TN / Neg
        
        snrROC[length] = sensitivity, specificity
    
    bam.close()
    # Announce the optimal SNR length using the J statistic
    print('The optimal SNR length is {}.'.format(max(
        snrROC, key = lambda l: (snrROC[l][0] * snrROC[l][1] if product
                                 else snrROC[l][0] + snrROC[l][1] - 1))))
    savePKL(out_SNREndROC, snrROC)
    
    return snrROC


def getSNRcovByGene(covLen, lenToSNRs, out_snrROC, out_transBaselineData,
                    out_db, bamfile):
    """This function gets Sensitivity & Specificity of coverage accounted for
    by SNRs (in exons only, but NOT exon-wise) in genes (expressed transcripts)
    in which they occur.
    
    Parameters
    ----------
    covLen : (int)
        The optimal length of coverage with respect to the potential priming
        location.
    BLdata : (tuple)
        (covTransByStrdRef, covExonsByStrdRef, Pos, Neg)
    lenToSNRs : (dict)
        { SNR length : [ SNR ] }
    bamfile : (str)
        Location of the bamfile to be scanned.

    Returns
    -------
    snrROC : (dict)
        { SNR length : (Sensitivity, Specificity) }
    """
    
    # If the file already exists, simply load
    if os.path.isfile(out_snrROC):
        snrROC = loadPKL(out_snrROC)
        return snrROC
    
    # Get the list of full transcripts by strd & ref
    BLdata = getBaselineData(out_transBaselineData, out_db, bamfile)
    covTransByStrdRef = BLdata[0]
    
    print('Extracting the exons up to the last {} bp of each expressed' \
          ' transcript...'.format(covLen))
    # Precalculate starts (one sorted list of tuples) for each expressed
    #  transcript based on the optimal endLen & save by strd & ref
    eachTransStartByStrdRef = collections.defaultdict(
        lambda: collections.defaultdict(list))
    for strd, covTransByRef in covTransByStrdRef.items():
        for ref, trans in covTransByRef.items():
            for Tran in trans:
                # Initialize the list of start exons
                startPieces = []
                # Initiate how much is left to progress to escape the end
                remaining = int(covLen)
                # For each exon (ordered last-to-first wrt/ strand)
                for eStart, eEnd in sorted(Tran.exons, reverse = strd):
                    # If end has not been escaped yet
                    if remaining > 0:
                        # Measure the size of the exon
                        exLen = eEnd - eStart
                        # If longer than the remainder, add the escaped part of
                        #  the exon
                        if exLen > remaining:
                            startPieces.append(
                                (eStart, eEnd - remaining) if strd
                                else (eStart + remaining, eEnd))
                        # Either way, decrease the remainder by the exon length
                        remaining -= exLen
                    # Otherwise add the entire exon as is
                    else:
                        startPieces.append((eStart, eEnd))
                # Append (not extend!) each exon list as a tuple
                eachTransStartByStrdRef[strd][ref].append(tuple(startPieces))
    
    # Initialize the SNRpieces found on covered exons by strd & ref - coverage
    #  of these will consitute TP & FP measurements
    SNRpiecesByStrdRef = collections.defaultdict(
        lambda: collections.defaultdict(list))
    # Initialize the associated measurements that'll be added to over SNR lens
    TP = 0
    FP = 0
    # Initialize the flattened transStarts whose transcripts had SNRpieces -
    #  coverage of these will constitute the baseline/denominator for Sens/Spec
    flatTransStartsByStrdRef = collections.defaultdict(
        lambda: collections.defaultdict(list))
    # Initialize the associated measurements that'll be added to over SNR lens
    Pos = 0
    Neg = 0
    # Initialize the ROC results dict
    snrROC = {}
    # Load the bamfile
    bam = pysam.AlignmentFile(bamfile, 'rb')
    # Go over SNRs by length, from longest to shortest
    for length in sorted(lenToSNRs, reverse = True):
        print('Getting Sensitivity & Specificity for SNRs of length' \
              ' {}...'.format(length))
        # For each length, go over each strd & ref of trans starts
        for strd in eachTransStartByStrdRef.keys():
            for refName in eachTransStartByStrdRef[strd].keys():
                # Initialize the list of SNR pieces *specific* for this SNR len
                newSNRpieces = []
                # Filter & extract the SNR pieces for each strd/ref
                for SNR in lenToSNRs[length]:
                    # Filter the SNRs by strd and ref
                    if SNR.strand == strd and SNR.record == refName:
                        start = SNR.start - covLen if strd else SNR.end
                        end = SNR.start if strd else SNR.end + covLen
                        newSNRpieces.append((start, end))
                # Determine which portions of the flattened SNR pieces are NEW
                #  for this strd/ref, as only those will contribute to
                #  - adding new transcript starts (exons) and
                #  - new TP/FP coverage information
                newSNRpieces = removeOverlaps(
                    flattenIntervals(newSNRpieces),
                    SNRpiecesByStrdRef[strd][refName])
                # These SNR pieces, specific for this SNR length, can overlap
                #  either the already included (flat list) trans starts, or
                #  some of the new ones (each list), or both (!)
                
                # Retain overlaps with previously added (flat list) exons:
                #  positive selection of SNR pieces based on old transStarts
                ovSNRpieces = getOverlaps(
                    newSNRpieces, flatTransStartsByStrdRef[strd][refName])
                
                # Initialize the list of tStart exons specific for this SNR len
                newTstarts = []
                # Initialize the list of tStarts whose exons were added
                toRemove = set()
                # For *each* transcript start that had not been added yet:
                #  positive selection of new transStarts and  SNR pieces based
                #  on mutual overlap (note: tran is a tuple of exons)
                for tran in eachTransStartByStrdRef[strd][refName]:
                    # Get the overlaps with SNR pieces specific for this length
                    SNRsInNewTransStart = getOverlaps(newSNRpieces, list(tran))
                    # If there are any,
                    if SNRsInNewTransStart != []:
                        # Add the overlap to the overlapped SNRpieces
                        ovSNRpieces.extend(SNRsInNewTransStart)
                        # Add the transcript start exons to the tStarts aggr
                        newTstarts.extend(tran)
                        # Add the transcript start to the removal list
                        toRemove |= {tran}
                # Remove the transcript starts from the original dict to
                #  increase efficiency (not to try to "rediscover" them again
                #  in the future, since they are already added)
                eachTransStartByStrdRef[strd][refName] = tuple(
                    tS for tS in eachTransStartByStrdRef[strd][refName]
                    if tS not in toRemove)
                
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
                    TP += np.sum(pieceCov)
                    FP += np.count_nonzero(pieceCov == 0)
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
                    Pos += np.sum(tCov)
                    Neg += np.count_nonzero(tCov == 0)
                # Itegrate the newTstarts into the aggregate dict & flatten
                flatTransStartsByStrdRef[strd][refName].extend(newTstarts)
                flatTransStartsByStrdRef[strd][refName] = flattenIntervals(
                    flatTransStartsByStrdRef[strd][refName])
                    
        # Calculate the Sens & Spec at this minimal SNR length and save
        if Pos != 0 and Neg != 0:
            TN = Neg - FP
            Sens = TP / Pos
            Spec = TN / Neg
            if not 0 <= Sens <= 1:
                raise Exception("Sensitivity is {}.".format(Sens))
            if not 0 <= Spec <= 1:
                raise Exception("Specificity is {}.".format(Spec))
            snrROC[length] = Sens, Spec
    
    bam.close()
    savePKL(out_snrROC, snrROC)
    
    return snrROC


def getStatsByGene(covLen, minSNRlen, lenToSNRs, out_geneStats, out_db,
                   out_transBaselineData, bamfile):
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
    BLdata = getBaselineData(out_transBaselineData, out_db, bamfile)
    covTransByStrdRef = BLdata[0]
    # Sort the covered transcripts by geneID
    covTransByGene = collections.defaultdict(list)
    for covTransByRef in covTransByStrdRef.values():
        for covTrans in covTransByRef.values():
            for covTran in covTrans:
                covTransByGene[covTran.geneID].append(covTran)
    # Initialize a dictionary for the results
    statsByGene = {}
    # Connect the gff database
    db = gffutils.FeatureDB(out_db, keep_order = True)
    # Count the total number of genes
    print('Counting the number of genes in the genome...')
    totalGenes = sum(1 for g in db.features_of_type(featuretype = 'gene'))
    # Load the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    # Initialize the progress tracker variables
    progressTracker = 0
    prog = 0
    # Go over all genes
    for gene in db.features_of_type(featuretype = 'gene'):
        progressTracker += 1
        newProg = round(progressTracker / totalGenes, 3)
        if newProg != prog:
            prog = newProg
            print('Getting gene coverage statistics... ({:.1%})'.format(prog))
            clear_output(wait = True)
        strd = gene.strand == '+'
        refName = gene.seqid
        # Extract all exons from this gene's covered transcripts
        geneExons = [(eS, eE) for Trans in covTransByGene[gene.id]
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
            for Trans in covTransByGene[gene.id]:
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
            for length, SNRs in SNRsByGeneLen[gene.id].items():
                if length >= minSNRlen:
                    for SNR in SNRs:
                        start = SNR.start - covLen if strd else SNR.end
                        end = SNR.start if strd else SNR.end + covLen
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
            statsByGene[gene.id] = (eCov, eLen, sCov, sLen, SNRlen, SNRnum,
                                    SNRstartLen, SNRstartNum)
    bam.close()
    savePKL(out_geneStats, statsByGene)
    
    return statsByGene