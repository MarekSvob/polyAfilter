#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 18:07:19 2020

@author: marek
"""

import os
import collections
import gffutils
import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO

from SNRdetection import loadPKL, savePKL, getGenomeLength, compDict, capsDict

# Default exclusivePairs
pairs = [
    ('gene', 'Intergenic'),
    ('transcript', 'Regulatory'),
    ('exon', 'Intron')
    ]


def getBaseComp(out_bases, fasta = None, showGC = True):
    """Function to obtain the base composition of a genome by scanning. Once
    this is done once for a given fasta file, the result is saved and from then
    on, always retreived. Optionally, the GC% is shown as a fact-check on
    correct processing.

    Parameters
    ----------
    fasta : (str), optional
        File path of the fasta reference file. The default is None.
    showGC : (bool), optional
        Switch for showing the GC% content of this genome. The default is True.

    Returns
    -------
    bases : (defaultdict)
        { base : count }
    """
    
    # Load the bases dictionary, if available
    if os.path.isfile(out_bases):
        bases = loadPKL(out_bases)
    # Otherwise, scan the genome to create it
    else:
        bases = collections.defaultdict(int)
    
        print('Scanning the genome for base composition...')
        # Make sure that the fasta file is open only temporarily
        with open(fasta, 'r') as genome:
            # Go over each base in each record in the genome and count each
            for record in SeqIO.parse(genome, "fasta"):
                for base in record.seq:
                    bases[base] += 1
        # Save the bases
        savePKL(out_bases, bases)
        # Optionally, show the GC% content of this genome
        if showGC:
            gc = bases['C'] + bases['c'] + bases['G'] + bases['g']
            sumKnownBases = bases['A'] + bases['a'] + bases['T'] + bases['t'] \
                + bases['C'] + bases['c'] + bases['G'] + bases['g']
            print(
                'Scanning finished, G/C content: {:.2%}'.format(
                    round(gc / sumKnownBases, 2)
                    )
                )
    
    return bases


def countsCheck(base, lengthToSNRcounts, lenToFeats, out_bases, fasta = None):
    """Check whether the sum of the bases in question obtained from SNR counts
    and scanning the genome are congruent.
    Note: 'bases' has to be a (defaultdict) because if non-caps bases are not
    defined, a simple (dict) would return an error.

    Parameters
    ----------
    base : (str)
        The base in question.
    lengthToSNRcounts : (dict)
        { SNRlength : SNRcount }
    out_bases : (str)
        File path of the pre-calculated bases file.
    fasta : (str), optional
        File path of the fasta reference file. The default is None.

    Returns
    -------
    None.
    """
    
    # Sum the number of bases contained in the SNR counts
    SNRbases = 0
    for k,v in lengthToSNRcounts.items():
        SNRbases += int(k) * int(v)
    
    # Get the base composition of the genome
    bases = getBaseComp(out_bases, fasta = fasta)
    # Get the number of bases (including complementary and non-caps) from the
    #  genome scan.
    scanned = bases[capsDict[base][0]] \
        + bases[capsDict[base][1]] \
            + bases[capsDict[compDict[base]][0]] \
                + bases[capsDict[compDict[base]][1]]
    
    featbases = sum([c*l for l,cs in lenToFeats.items() for c in cs.values()])
    
    if SNRbases == scanned == featbases:
        print(
            'The total number of {}/{} bases ({:,}) '\
                'checks out across SNRs, feature sets, and genome!'.format(
                base,
                compDict[base],
                SNRbases
                )
            )
    else:
        print(
            '{}/{} bases in SNRs: {:,}; from genome scanning {:,}; '\
                'from feature counting: {:,}'.format(
                base,
                compDict[base],
                SNRbases,
                scanned,
                featbases
                )
            )
    
    


def SNRcountTable(base, lengthToSNRcounts, out_bases, fasta = None):
    """Obtain a df of SNR Length, Observed counts, and Observed/Expected
    ratios.
    Note: 'bases' has to be a (defaultdict) because if non-caps bases are not
    defined, a simple (dict) would return an error.

    Parameters
    ----------
    base : (str)
        The base in question.
    lengthToSNRcounts : (dict)
        { SNRlength : SNRcount }
    out_bases : (str)
        File path of the pre-calculated bases file.
    fasta : (str)
        File path of the fasta reference file. The default is None.

    Returns
    -------
    df : (pandas.dataframe)
        The table with calculated values.
    """
    
    # Get the base composition of the genome
    bases = getBaseComp(out_bases, fasta = fasta)
    
    sumKnownBases = bases['A'] + bases['a'] + bases['T'] + bases['t'] \
        + bases['C'] + bases['c'] + bases['G'] + bases['g']
    
    # Calculate the frequency of the base in question and its complement
    pBase = (bases[capsDict[base][0]] \
             + bases[capsDict[base][1]]) / sumKnownBases
    pComp = (bases[capsDict[compDict[base]][0]] \
             + bases[capsDict[compDict[base]][1]]) / sumKnownBases
    
    # Obtain the table data as a list of 3 respective lists: SNRlength,
    #  Observed absolute # of SNRs, and O/E ratio
    polyLen = []
    obs = []
    oe = []
    for k,v in lengthToSNRcounts.items():
        polyLen.append(int(k))
        obs.append(int(v))
        oe.append(
            (int(v)/int(sumKnownBases)) \
                / ( (1-pBase) * pBase**int(k) * (1-pBase) \
                   + (1-pComp) * pComp**int(k) * (1-pComp) )
            )
    table_data = [polyLen,obs,oe]
    # Construct the df with the desired column names & types
    df = pd.DataFrame(data = pd.DataFrame(data = table_data).T)
    df.columns = ('SNR Length','Observed', 'O/E')
    df['Observed'] = df['Observed'].astype(int)
    df['SNR Length'] = df['SNR Length'].astype(int)
    df.sort_values(by = ['SNR Length'], inplace = True)
    
    return df


def SNRfeatureSets(lenToFeats):
    """Function that compiles a bool df of Features vs. FeatureSets.
    
    Parameters
    ----------
    lengthToSNRs : (dict)
        The dict of SNRs detected, sorted by length.

    Returns
    -------
    df : (dataframe)
        A dataframe of features vs. featureSets, where the presence of the
        former in the latter is expressed as a bool.
    """
    
    # Get the set of (frozen)sets of features represented among the SNRs
    featureSets = list({s for c in lenToFeats.values() for s in c.keys()})
    # Flatten this set of features represented
    features = list({f for featureSet in featureSets for f in featureSet})
    
    # Create a { feature : [ bool ] } dict to describe representation of each
    #  feature in each featureset and use this to create a df
    h_data = collections.defaultdict(list)
    for featureSet in featureSets:
        for feature in features:
            h_data[feature].append(feature in featureSet)
    df = pd.DataFrame(data = pd.DataFrame(data = h_data).T)
    
    # Sort the rows (features) and columns (featureSets) from most to least
    rows = list(df.index)
    rows.sort(key = lambda x: df.sum(axis = 1)[x], reverse = True)
    cols = list(df.columns)
    cols.sort(key = lambda x: df.sum(axis = 0)[x], reverse = True)
    df = df.loc[rows, cols]

    return df


def getColNames(df, exclusivePairs = pairs, other = 'Exon'):
    """Function to create a label to each df column (featureSet) according to
    the exclusivePairs - i.e., if feature is absent, label is assigned.
    
    Parameters
    ----------
    df : (dataframe)
        A bool dataframe describing the presence of Features in FeatureSets.
    exclusivePairs : (list), optional
        A list of tuples with features and the respective labels assigned to
        mark their absence, such that [ (feature, label) ]. The default is
        [('gene','Intergenic'),('transcript','Regulatory'),('exon','Intron')].
    other : (str), optional
        The label given if all features are present. The default is 'Exon'.

    Returns
    -------
    named_cols : (list)
        The list of df labels.
    """
    
    named_cols = []
    # For each column (featureSet), check each feautre in order, one-by-one.
    for col in df.columns:
        for pair in exclusivePairs:
            # If the feature is absent, assign the label & break the inner
            #  loop, continuing with the next column.
            if not df.loc[pair[0], col]:
                named_cols.append(pair[1])
                break
        #  If no labels have been assigned, assign the remaining 'other' label.
        else:
            named_cols.append(other)
            
    return named_cols


def SNRlabelProps(lenToFeats, exclusivePairs = pairs, other = 'Exon'):
    """Function that prepares a df of proportions of labeled SNRs for each
    length.    

    Parameters
    ----------
    lenToFeats : (dict)
        { length : { { feature } : count } }
    exclusivePairs : (list), optional
        A list of tuples with features and the respective labels assigned to
        mark their absence, such that [ (feature, label) ]. The default is
        [('gene','Intergenic'),('transcript','Regulatory'),('exon','Intron')].
    other : (str), optional
        The label given if all features are present. The default is 'Exon'.

    Returns
    -------
    df_p : (DataFrame)
        Table of proportions of labeled SNRs for each length.
    """
    
    # Initiate the data matrix of zeros, SNRlength x label
    SNRlabels_data = np.zeros(
        (max(lenToFeats.keys()) + 1, len(exclusivePairs) + 1),
        dtype = int
        )
    
    # For each length, add up the feats corresponding to the respective labels
    for length, featureCounts in lenToFeats.items():
        for featureSet, count in featureCounts.items():
            for i in range(len(exclusivePairs)):
                if exclusivePairs[i][0] not in featureSet:
                    index = i
                    break
            else:
                index = len(exclusivePairs)
            SNRlabels_data[length, index] += count

    # Derive the labels from the input
    colNames = [p[1] for p in exclusivePairs]
    colNames.append(other)
    df = pd.DataFrame(data = SNRlabels_data, columns = colNames)
    # Remove all the lengths that have less than 10 SNRs in them
    df.drop(
        labels = [i for i in df.index if df.sum(axis = 1)[i] in range(10)],
        axis = 0,
        inplace = True
        )
    # Get proportions, instead of absolute values, for each SNR length
    props = df.values / df.values.sum(axis = 1)[:, np.newaxis]
    df_p = pd.DataFrame(data = props, columns = df.columns, index = df.index)
    
    return df_p


def getStrandedFeats(out_strandedFeats, out_db, featsOfInterest):
    """Function to scan the GTF file and flatten the feats, preserving their
    strandedness.    

    Parameters
    ----------
    out_strandedFeats : (str)
        Path to saved output.
    out_db : (str)
        Path to the database.
    featsOfInterest : (iter)
        An iterable (list or tuple) containing the features (str) of interest.

    Returns
    -------
    flatStrandedFeats : (dict)
        { featureType+/- : [ (ref, start, end) ] }

    """
    
    # If the flat feats have been created before, load, otherwise process
    if os.path.isfile(out_strandedFeats):
        flatStrandedFeats = loadPKL(out_strandedFeats)
        return flatStrandedFeats
    
    # Connect the db
    db_conn = gffutils.FeatureDB(out_db, keep_order = True)
    
    # Initiate the variables needed
    strands = ('+', '-')
    flatStrandedFeats = {}
    
    for featType in featsOfInterest:
        # Go over each strand separately
        for strd in strands:
            print('Flattening {}{}s'.format(featType, strd))
            featList = []
            # Iterate through ALL features of this type, for each strand.
            for feat in db_conn.all_features(
                    featuretype = featType,
                    strand = strd
                    ):
                # Save in a tuple: (ref, start, end)
                nFeat = (feat.seqid, feat.start, feat.end)
                # Go over the previous ones, check for overlap on the same
                #  ref and save all the overlaps found in a list.
                overlaps = []
                for oFeat in featList:
                    # If the ref is the same and (
                    #  (the new Feat's beginning is within the old Feat)
                    #  or (the new Feat's end is within the oldFeat)
                    #  or (the newFeat spans the oldFeat)
                    #  ), add old Feat to the overlaps list
                    if nFeat[0] == oFeat[0] and (
                        (nFeat[1] >= oFeat[1] and nFeat[1] <= oFeat[2]) \
                        or (nFeat[2] >= oFeat[1] and nFeat[2] <= oFeat[2]) \
                        or (nFeat[1] < oFeat[1] and nFeat[2] > oFeat[2])
                        ):
                        overlaps.append(oFeat)
                #  If overlaps have been found, merge with the new one
                if overlaps != []:
                    # Initialize the start & end of the merged feature
                    #  using the new Feat (not in the overlaps list)
                    ref = nFeat[0]; start = nFeat[1]; end = nFeat[2]
                    # Go over all the overlaps & update the start & end of
                    #  the merged feature as necessary
                    for ft in overlaps:
                        if ft[1] < start:
                            start = ft[1]
                        if ft[2] > end:
                            end = ft[2]
                        # When done considering this feature, remove it
                        #  from the master list
                        featList.remove(ft)
                    # When done, add the new merged feature to the list
                    featList.append((ref, start, end))
                # If no overlap has been found, add new as a new one
                else:
                    featList.append(nFeat)
            # Save the feats in the master list to be saved
            flatStrandedFeats['{}{}'.format(featType, strd)] = sorted(featList)
    # Save the flattened feats to speed up the process in the future
    savePKL(out_strandedFeats, flatStrandedFeats)
    
    return flatStrandedFeats


def normalizeLabels(
        df,
        out_strandedFeats,
        out_db = None,
        fasta = None,
        exclusivePairs = pairs
        ):
    """Function to calculate what proportion of the genome is covered by each
    of the respective labels and subsequenly normalize the df by these
    proportions. If the flattened feats have not been processed previously,
    they will be processed and saved (can take up to 5 hrs).
    
    Parameters
    ----------
    df : (dataframe)
        The df of measured SNR labels to be normalized.
    out_strandedFeats : (str)
        The location of the PKL file containing the dict of flattened feats.
    fasta : (str), optional
        The location of the FASTA reference sequence file. The default is None.
    out_db : (str), optional
        The location of the reference annotation database file. The default is
        None.
    exclusivePairs : (list), optional
        A list of tuples with features and the respective labels assigned to
        mark their absence, such that [ (feature, label) ]. The default is
        [('gene','Intergenic'),('transcript','Regulatory'),('exon','Intron')].

    Returns
    -------
    df_norm : (dataframe)
        Normalized dataframe
    """
    
    # Extract the features of interest
    featsOfInterest = [p[0] for p in exclusivePairs]
    # Initiate the counter of total feature length
    regions = collections.defaultdict(int)
    
    # First, if not done yet, flatten the relevant features, using the
    #  simplified peak-merging strategy from Biostats Final Project
    # Note that GFF features are 1-based, closed intervals
    flatStrandedFeats = getStrandedFeats(
        out_strandedFeats,
        out_db,
        featsOfInterest
        )

    # Add up the cumulative feature length, summing up both strands
    for k,feats in flatStrandedFeats.items():
        # 1-based, closed intervals [start, end]
        regions[k[:-1]] += sum([feat[2] - feat[1] + 1 for feat in feats])
    
    # Get the genome length for this fasta file
    genLen = getGenomeLength(fasta)

    # Now calculate the proportion of genome covered by each label - e.g.,
    #  Intergenic = genome - genes
    #  Regulatory = genes - transcripts
    #  Intron = transcripts - exons
    #  Exon = exons

    # The first label in exclusivePairs depends on the length of the genome
    labelProps = np.array(
        [(genLen*2 - regions[exclusivePairs[0][0]]) / (genLen*2)]
        )
    # Add values for the subsequent labels in exclusivePairs
    for i in range(1,len(exclusivePairs)):
        labelProps = np.concatenate((labelProps, np.array([
            (regions[exclusivePairs[i-1][0]] - regions[exclusivePairs[i][0]]) \
                / (genLen*2)
            ])))
    # Add the last one, which is independent of other labels
    labelProps = np.concatenate(
        (labelProps, np.array([regions[exclusivePairs[-1][0]] / (genLen*2)]))
        )
    
    # Normalize the df using these proportions
    df_norm = pd.DataFrame(
        data = df.values / labelProps[np.newaxis, :],
        columns = df.columns,
        index = df.index
        )
    
    return df_norm


def filterReads(
        bamfile,
        filt_bamfile,
        flatStrandedFeats,
        featType = 'exon',
        concordant = True
        ):
    """Function to create a bam file from a subset of reads that are aligned
    to a specific feature type on the same or opposite strand.

    Parameters
    ----------
    bamfile : (str)
        The path to the original bamfile.
    filt_bamfile : (str)
        The path to the filtered bamfile.
    flatStrandedFeats : (dict)
        { featureType+/- : [ (ref, start, end) ] } where [start, end] is a
        1-based closed interval.
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
        return    
    
    # Connect to the bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # Initiate the varuables needed
    # This needs to start off as a set to avoid double-adding reads aligned to
    #  two features of the same type
    filt_reads = set()
    current_ref = ''
    strds = ('+', '-')
    
    # Get only reads on the same strand as the feature in question
    for strd in strds:
        # Cycle over the feats under the appropriate key
        for feat in sorted(flatStrandedFeats['{}{}'.format(featType, strd)]):
            if current_ref != feat[0]:
                current_ref = feat[0]
                print(
                    'Processing {} strand of reference {}'.format(
                        strd,
                        current_ref
                        )
                    )
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
    

def getCoveragePerSNR(
        lenToSNRs,
        bamfile,
        SNRlengths,
        window = 2000,
        concordant = True,
        SNRfeat = 'transcript'
        ):
    """Obtain RNA-seq coverage in the vicinity of SNRs of given length.
    
    Parameters
    ----------
    lengthToSNRs : (dict)
        { length : [SNR] }
    bamfile : (str)
        Path to the (pre-filtered) bamfile.
    SNRlengths : (iter)
        Any iterable (list, tuple, etc.) containing the SNR lengths of interest.
    window : (int), optional
        Total size of the window around the SNR covered. The default is 2000.
    concordant : (bool), optional
        Switch between concordant & discordant coverage. The default is True.
    SNRfeat : (str), optional
        The feature by which SNRs are selected. The default is 'transcript'.

    Returns
    -------
    coverageByLen : (dict)
        { SNRlength : ( SNRcount, np.array ) } 
    """

    # Attach a pre-filtered bam file
    bam = pysam.AlignmentFile(bamfile, 'rb')
    
    # Initialize the dict to collect all the data
    coverageByLen = {}
    # Go over the SNRs by length, starting with the longest
    for length in sorted(lenToSNRs.keys(), reverse = True):
        if length in SNRlengths:
            # For each length, initialize the count & array for the window coverage
            SNRcount = 0
            coverage = np.zeros((window), 'L')
            print('Going over SNR{}s...'.format(length))
            # Note that SNRs are 0-based
            for SNR in lenToSNRs[length]:
                # Only consider SNRs in transcripts, as those are the
                #  sources of internal (concordant or discordant) priming
                if (concordant and SNRfeat in SNR.concFeats) \
                    or (not concordant and SNRfeat in SNR.discFeats):
                    # Count the SNR
                    SNRcount += 1
                    # Get the mid point of the SNR wrt/ the reference
                    mid = (SNR.end + SNR.start) // 2
                    start = int(mid - window/2)
                    stop = int(mid + window/2)
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
                    refCoverage = np.sum(
                        bam.count_coverage(
                            contig = SNR.record,
                            start = corrStart,
                            stop = corrStop,
                            quality_threshold = 0,
                            read_callback = lambda r: concordant == \
                                (r.is_reverse != SNR.strand)
                            ),
                        axis = 0
                        )
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
                    #  orientation wrt/ the transcript direction
                    if SNR.strand == concordant:
                        coverage += refCoverage
                    else:
                        coverage += refCoverage[::-1]
            if SNRcount != 0:
                coverageByLen[length] = (SNRcount, coverage)
   
    bam.close()
    
    return coverageByLen


def getExpectedCoverage(bamfile, flatStrandedFeats, SNRfeat = 'transcript'):
    """Get the per-base expected RNA-seq coverage across the given featuretype.
    
    Parameters
    ----------
    bamfile : (str)
        Path to the (pre-filtered) bamfile.
    flatStrandedFeats : (dict)
        { featureType : [ (ref, start, end) ] }
    SNRfeat : (str), optional
        The featureType from which SNRs are drawn. The default is 'transcript'.

    Returns
    -------
    expected : (float)
        The expected per-base coverage from the RNA-seq data.
    """
    
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
        print('Scanning {}{}s for total length...'.format(SNRfeat, strd))
        for feat in flatStrandedFeats['{}{}'.format(SNRfeat, strd)]:
            # Note that pysam is 0-based; my vars are 1-based (like the gtf)
            totalLength += feat[2] - feat[1] + 1
    
    # Normalize the coverage by #SNRs, total coverage, and total length
    expected = totalCoverage/totalLength
    
    return expected