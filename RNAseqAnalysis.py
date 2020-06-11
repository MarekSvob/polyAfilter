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

expCovConc = 0
expCovDisc = 0


class NoncGene:
    """An object that stores information about a gene identified to be
    non-canonically covered.
    """    
    
    def __init__(self, geneF, cov, prop, trans, SNRs):
        self.geneFeat = geneF   # gffutils.Feature, gff db gene in quesion
        self.cov = cov          # np.array, the bp-wise coverage of this gene
        self.prop = prop        # float, the non-canonical proportion of cov
        self.trans = trans      # list, [ gffutils.Feature ] of overlapping Ts
        self.SNRs = SNRs        # dict, { SNRlength : [ SNR ] }
        
    def __str__(self):
        # Displays the base, number of mismatches, and location
        return '{} NoncGene @ {} strand'.format(
            self.geneFeat.ID,
            self.geneFeat.strand
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
    for strd in (True, False):
        print(
            'Scanning {}{}s for total length...'.format(
                '+' if strd else '-',
                baselineFeat
                )
            )
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


def getNonCanCovGenes(
        out_NonCanCovGenes,
        lenToSNRs,
        out_db,
        bamfile,
        lastBP = 250,
        covThreshold = 0.05,
        minLen = 0
        ):
    """Function to scan each gene with SNRs and determine whether its RNA-seq
    coverage can be accounted for by canonical reads in exon-wise ends.

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
        The minimal proportion of non-canonical coverage to include the gene.
        The default is 0.05.
    minLen : int, optional
        Minimal SNR length that needs to be present for a gene to be
        considered. The default is 0.

    Returns
    -------
    nonCanCovGenes : (list)
        Genes found to have a significant proportion of coverage outside of
        the canonically expected area.
    """
    # For each gene, what proportion of (concordant) coverage is outside of
    #  the last 250 bp of all exon-wise transcripts?
    
    # If the file already exists, simply load
    if os.path.isfile(out_NonCanCovGenes):
        nonCanCovGenes = loadPKL(out_NonCanCovGenes)
        return nonCanCovGenes
    
    nonCanCovGenes = []
    
    # Sort SNRs by gene
    geneLenToSNRs = getGeneLenToSNRs(lenToSNRs, concordant = True)
    db = gffutils.FeatureDB(out_db, keep_order = True)
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # For each gene that has an SNR detected
    for geneID, lenToSNRs in geneLenToSNRs.items():
        # If no SNRs of sufficient length are found, skip
        if all(length < minLen for length in lenToSNRs.keys()):
            continue
        # Get the total bp-wise coverage for the gene
        geneFeat = db[geneID]
        geneBPcov = np.sum(
            bam.count_coverage(
                contig = geneFeat.seqid,
                start = geneFeat.start,
                stop = geneFeat.end,
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
        endings = []
        endCoverage = 0
        # Go over ALL transcripts overlapping this gene (even from other genes)
        #  add up the coverage of exon-wise last X bp
        for trans in db.features_of_type(
                featuretype = 'transcript',
                limit = (geneFeat.seqid, geneFeat.start, geneFeat.end),
                strand = geneFeat.strand
                ):
            # Extract the exons as (start, end)
            exons = sorted(
                [
                    (e.start, e.end) for e in db.children(
                        trans,
                        featuretype = 'exon'
                        )
                    ]
                )
            # Determine the transcripts's exon-wise end
            remaining = lastBP
            covStart = None
            strd = geneFeat.strand == '+'
            i = -strd
            f = 1 if strd else -1   # multiplication factor
            while covStart is None:
                exLen = exons[i][1] - exons[i][0] + 1
                if exLen < remaining:
                    remaining -= exLen
                    # If we have gone over all exons, take the beginning of
                    #  the transcript; otherwise keep going
                    if len(exons) <= abs(i - 1 * f + strd) :
                        covStart = trans.start if strd else trans.end
                    else:
                        i -= 1 * f
                else:
                    covStart = exons[i][strd] - remaining * f
            # Only add the transcript if the transcripts's exon-wise end is
            #  not beyond the gene end
            if (strd and covStart <= geneFeat.end) \
                or (not strd and covStart >= geneFeat.start):
                    # Save this transcript feature in the list
                    transcripts.append(trans)
                    # Add an ending whose coverage is to be assessed
                    if strd:
                        endings.append(
                            (covStart, min(trans.end, geneFeat.end))
                            )
                    else:
                        endings.append(
                            (max(trans.start, geneFeat.start), covStart)
                            )
                    
        # Once done going over the transcripts, add up the coverage of endings
        for (covS, covE) in endings:
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
        # Add this as an interesting gene only if the unaccounted-for coverage
        #  (not in endings) exceeds the threshold
        prop = (geneCoverage - endCoverage) / geneCoverage
        [biotype] = geneFeat.attributes['gene_biotype']
        if prop >= covThreshold:
            nonCanCovGenes.append(
                NoncGene(geneFeat, geneBPcov, prop, transcripts, lenToSNRs)
                )

    bam.close()
    savePKL(out_NonCanCovGenes, nonCanCovGenes)    
    
    return nonCanCovGenes