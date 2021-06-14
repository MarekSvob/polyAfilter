#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 18:07:19 2020

@author: marek
"""

import os
import collections
import gffutils
import logging
import pandas as pd
import numpy as np
from Bio import SeqIO

from SNRdetection import loadPKL, savePKL, getGenomeLength, compDict, capsDict


logging.basicConfig(level = logging.INFO,
                    format = '%(asctime)s - %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Default exclusivePairs
pairs = [('gene', 'Intergenic'),
         ('transcript', 'Regulatory'),
         ('exon', 'Intron')]


def getBaseComp(out_bases, fasta = None, showGC = True):
    """Function to obtain the base composition of a genome by scanning. Once
    this is done once for a given fasta file, the result is saved and from then
    on, always retreived. Optionally, the GC% is shown as a fact-check on
    correct processing.

    Parameters
    ----------
    out_bases : (str)
        File path of the pre-calculated bases file.
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
    else:
        if fasta is None:
            raise ValueError('A path to either a fasta file or a valid '
                             'pre-calculated base composition file needs to be'
                             ' provided.')
        # Otherwise, scan the genome to create it
        bases = collections.defaultdict(int)
    
        logger.info('Scanning the genome for base composition...')
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
        logger.info('Scanning finished, G/C content: '
                    f'{round(gc / sumKnownBases, 2):.2%}.')
    
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
    lenToFeats : (dict)
        { SNRlength : { { feature } : count } }
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
    scanned = (sum(bases[cap] for cap in capsDict[base])
               + sum(bases[cap] for cap in capsDict[compDict[base]]))
    # scanned = (bases[capsDict[base][0]]
    #            + bases[capsDict[base][1]]
    #            + bases[capsDict[compDict[base]][0]]
    #            + bases[capsDict[compDict[base]][1]])
    
    featbases = sum([c*l for l,cs in lenToFeats.items() for c in cs.values()])
    
    if SNRbases == scanned == featbases:
        logger.info(f'The total number of {base}/{compDict[base]} bases '
                    f'({SNRbases:,d}) checks out across SNRs, feature sets, '
                    'and genome!')
    else:
        logger.info(f'{base}/{compDict[base]} bases in SNRs: {SNRbases:,d}; '
                    f'from genome scanning {scanned:,d}; '
                    f'from feature counting: {featbases:,d}.')


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
    
    sumKnownBases = (bases['A'] + bases['a'] + bases['T'] + bases['t']
                     + bases['C'] + bases['c'] + bases['G'] + bases['g'])
    
    # Calculate the frequency of the base in question and its complement
    pBase = sum(bases[cap] for cap in capsDict[base]) / sumKnownBases
    # pBase = (bases[capsDict[base][0]]     # capsDict changed to dict(sets)
    #          + bases[capsDict[base][1]]) / sumKnownBases
    pComp = sum(bases[cap] for cap in capsDict[compDict[base]]) / sumKnownBases
    # pComp = (bases[capsDict[compDict[base]][0]]
    #          + bases[capsDict[compDict[base]][1]]) / sumKnownBases
    
    # Obtain the table data as a list of 3 respective lists: SNRlength,
    #  Observed absolute # of SNRs, and O/E ratio
    polyLen = []
    obs = []
    oe = []
    for k,v in lengthToSNRcounts.items():
        polyLen.append(int(k))
        obs.append(int(v))
        oe.append((int(v)/int(sumKnownBases))
                  / ((1-pBase) * pBase**int(k) * (1-pBase)
                     + (1-pComp) * pComp**int(k) * (1-pComp)))
    table_data = [polyLen, obs, oe]
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
    lenToFeats : (dict)
        { length : { { feature } : count } }

    Returns
    -------
    df : (dataframe)
        A dataframe of features vs. featureSets, where the presence of the
        former in the latter is expressed as a bool.
    """
    
    # Get the list of (frozen)sets of features represented among the SNRs
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
        dtype = int)
    
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
    df.drop(labels = [i for i in df.index if df.sum(axis = 1)[i] in range(10)],
            axis = 0,
            inplace = True)
    # Get proportions, instead of absolute values, for each SNR length
    props = df.values / df.values.sum(axis = 1)[:, np.newaxis]
    df_p = pd.DataFrame(data = props, columns = df.columns, index = df.index)
    
    return df_p


def flattenIntervals(intervals):
    """Improved helper function to flatten a list of overlapping intervals in
    linear time. Note that this requires 0-based coordinates.

    Parameters
    ----------
    intervals : (list)
        [ (start, end) ]

    Returns
    -------
    flat : (list)
        [ (start, end) ]
    """
    
    # Initialize the list of merged intervals
    flat = []
    # If the list is empty, return
    if intervals == []:
        return flat
    # Sort the intervals by start
    intervals.sort(key = lambda interval: interval[0])
    # Initialize the first 'current' interval
    currentStart, currentEnd = intervals[0]
    # Make sure that this is indeed a valid interval
    if currentStart >= currentEnd:
        raise Exception(f"Interval's start [{currentStart}] >= end "
                        f"[{currentEnd}]!")
    for iStart, iEnd in intervals[1:]:
        # Make sure that this is indeed a valid interval
        if iStart >= iEnd:
            raise Exception(f"Interval's start [{iStart}] >= end [{iEnd}]!")
        # If the new interval overlaps the existing, merge by updating the end.
        # Note that when iStart = currentEnd, this is technically an adjacency
        #  (not overlap) but even then a merge is desirable, nevertheless.
        # Also note that because the intervals have been sorted,
        #  currentStart > iStart can never be the case
        if iStart <= currentEnd:
            currentEnd = max(iEnd, currentEnd)
        # Otherwise add the previous (no longer overlapped) merged interval
        #  into the list and create a new 'current' interval.
        else:
            flat.append((currentStart, currentEnd))
            currentStart = iStart
            currentEnd = iEnd
    # Don't forget to add the last interval!
    flat.append((currentStart, currentEnd))
    
    # Check that the intervals in each list are indeed not overlapping
    nonOverlapCheck(flat)
        
    return flat


def getOverlaps(aIntervals, bIntervals):
    """Helper function to get the mutual overlaps between two lists of
    intervals in linear time. Note that this function assumes flattened
    (internally non-overlapping) 0-based intervals.
    
    Parameters
    ----------
    aIntervals : (list)
        [ (start, end) ]
    bIntervals : (list)
        [ (start, end) ]

    Returns
    -------
    overlaps : (list)
        [ (start, end) ]
    """
    
    # Sort the two lists to be in the ascending order
    aIntervals.sort(key = lambda interval: interval[0])
    bIntervals.sort(key = lambda interval: interval[0])
    # Check that the sorted intervals in each list are indeed not overlapping
    for intervals in (aIntervals, bIntervals):
        nonOverlapCheck(intervals)
    # Initialize the output and bInt index, and extract the first bInt info
    overlaps = []
    # If bIntervals is empty, return
    if bIntervals == []:
        return overlaps
    bI = 0
    bStart, bEnd = bIntervals[bI]
    # Once bI reaches lenOfB, return overlaps
    lenOfB = len(bIntervals)
    # Test each possible pair of bIntervals & bIntervals in linear time
    for aStart, aEnd in aIntervals:
        # For this given aInt, keep going over the bIntervals as long as they
        #  do not occur after the aInt, including when bStart = aEnd (adjacent)
        #  once the bInterval occurs after the aInt, go to next aInt
        while bStart < aEnd:
            # If the bInterval is entirely before the aInt, go to the next one
            #  without saving; when bEnd = aStart, they are only adjacent
            if bEnd <= aStart:
                bI += 1
                if bI < lenOfB:
                    bStart, bEnd = bIntervals[bI]
                else:
                    return overlaps
            # Othewise see which of the following overlaps is found & save
            # If the bInterval start is not inside the aInt
            elif bStart <= aStart:
                # If bEnd is inside the aInt, just add up to bEnd
                if bEnd < aEnd:
                    overlaps.append((aStart, bEnd))
                    bI += 1
                    if bI < lenOfB:
                        bStart, bEnd = bIntervals[bI]
                    else:
                        return overlaps
                # Otherwise the bInterval must overhang the entire aInt, so add
                #  the entire aInt, BUT DO NOT move onto the next bInterval
                #  until next aInt is also checked, as it also can be
                #  overlapped by this bInterval!
                else:
                    overlaps.append((aStart, aEnd))
                    break
            # Otherwise the bStart must be inside the aInt; if the bEnd is also
            #  inside the aInt, add the entire bInterval
            elif bEnd < aEnd:
                overlaps.append((bStart, bEnd))
                bI += 1
                if bI < lenOfB:
                    bStart, bEnd = bIntervals[bI]
                else:
                    return overlaps
            # The only option left is that the bInterval end is outside, while
            #  the bInterval start is inside, so add the bInterval up to the
            #  aEnd BUT again, do not move onto the following bInterval yet!
            else:
                overlaps.append((bStart, aEnd))
                break
                
    return overlaps


def removeOverlaps(aIntervals, bIntervals):
    """Helper function to remove the mutual overlaps of two lists of intervals
    from the first one of them (aIntervals) in linear time. This function uses
    similar general logic as getOverlaps(). Note that this function assumes
    two lists of flattened (internally non-overlapping) 0-based intervals.
    
    Parameters
    ----------
    aIntervals : (list)
        [ (start, end) ]
    bIntervals : (list)
        [ (start, end) ]

    Returns
    -------
    overlaps : (list)
        [ (start, end) ]
    """
    
    # Sort the two lists to be in the ascending order
    aIntervals.sort(key = lambda interval: interval[0])
    bIntervals.sort(key = lambda interval: interval[0])
    # Check that the sorted intervals in each list are indeed not overlapping
    for intervals in (aIntervals, bIntervals):
        nonOverlapCheck(intervals)
    # Initialize the output, bInt index, and max bInt index
    nonOverlaps = []
    bI = 0
    maxbI = len(bIntervals)
    # Test each possible pair of aIntervals & bIntervals in linear time
    for aStart, aEnd in aIntervals:
        # Initialize the non-overlap start of the aInterval to be modified
        # Note: when the end of the aInterval is modified, it is added right
        #  away and next aInterval is evaluated without changing the bInterval
        noStart = aStart
        # Once bI reaches maxbI, overlaps are no longer possible and remainder
        #  of all aIntervals are added (this also means that when bIntervals is
        #  an empty list, all aIntervals are returned as is)
        while bI in range(maxbI):
            # Initialize the overlapping piece start & end
            bStart, bEnd = bIntervals[bI]
            # For this given aInt, keep going over the bIntervals as long as
            #  they do not occur after the aInt, (including when bStart = aEnd
            #  (adjacent)). Once the bInterval occurs after the aInt, add the
            #  non-overlapped portion of aInt and go to next w/o increasing bI
            if bStart >= aEnd:
                nonOverlaps.append((noStart, aEnd))
                break
            # If the bInterval is entirely before the aInt, go to the next one
            #  without modifying; when bEnd = noStart, they are only adjacent
            elif bEnd <= noStart:
                bI += 1
            # Otherwise see which of the following overlaps is found & modify
            # If the bInterval start is before the aInt
            elif bStart <= noStart:
                # If the bEnd is inside the aInt, remove up to bEnd and move on
                #  to the next bInt
                if bEnd < aEnd:
                    noStart = bEnd
                    bI += 1
                # Otherwise the bInterval must overhang the entire aInt, so
                #  skip adding the rest of aInt, without moving onto the next
                #  bInterval, as the next aInt may be overlapped by this bInt
                else:
                    break
            # Otherwise the bStart must be inside the aInt; if the bEnd is also
            #  inside the aInt, add the left aInt overhang and update the start
            elif bEnd < aEnd:
                nonOverlaps.append((noStart, bStart))
                noStart = bEnd
                bI += 1
            # The only option left is that the bEnd is outside, while the
            #  bStart is inside, so add the aInt up to the bStart BUT again,
            #  do not move onto the following bInterval just yet!
            else:
                nonOverlaps.append((noStart, bStart))
                break
        # Once a bInt outside the aInterval is reached, add the non-overlapped
        #  portion of the current aInterval and go to the next aInterval
        else:
            nonOverlaps.append((noStart, aEnd))
            
    return nonOverlaps


def nonOverlapCheck(intervals):
    """For a given sorted list of intervals, check that they are not
    overlapping.
    
    Parameters
    ----------
    intervals : (list)
        [ (start, end) ]

    Returns
    -------
    None.
    """
    oldEnd = 0
    for newStart, newEnd in intervals:
        if newStart < oldEnd:
            raise Exception(f'The intervals overlap! [newStart = {newStart}, '
                            f'newEnd = {newEnd}, oldEnd = {oldEnd}]')
        # Note that if newStart == oldEnd, they are adjacent
        else:
            oldEnd = newEnd
        


def getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db, featsOfInterest):
    """Function to scan the GTF file and flatten the feats, preserving their
    strandedness. Converts 1-based to 0-based feats.   

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
    flatFeats : (dict)
        { featuretype : { strand : { ref : [ ( start, end ) ] } } }
        where the T/F bool indicates +/- strand, respectively.
    """
    
    # If the flat feats have been created before, load
    if os.path.isfile(out_strandedFeats):
        flatFeats = loadPKL(out_strandedFeats)
        # If all feat types needed have been processed, simply return;
        #  otherwise new feature(s) will be added below
        if all(featType in flatFeats.keys() for featType in featsOfInterest):
            return flatFeats
    # Only initiate a completely new dict if none has been made before
    else:
        flatFeats = collections.defaultdict(
            lambda: collections.defaultdict(
                lambda: collections.defaultdict()))
        
    # Connect the db
    db_conn = gffutils.FeatureDB(out_db, keep_order = True)
    
    # Go over each strand separately
    for featType in featsOfInterest:
        for strd in (True, False):
            logger.info(f'Flattening {"+" if strd else "-"}{featType}s...')
            featsByRef = collections.defaultdict(list)
            # Iterate through ALL features of this type, for each strand.
            for feat in db_conn.all_features(featuretype = featType,
                                             strand = '+' if strd else '-'):
                # Add this new feat as a tuple by ref: (start, end), 0-based
                featsByRef[feat.seqid].append((feat.start - 1, feat.end))
            # Once done for this Type & strd, flatten for each ref & save
            for ref, feats in featsByRef.items():
                flatFeats[featType][strd][ref] = flattenIntervals(feats)

    # Save the flattened feats to speed up the process in the future (note
    #  that this may overwrite a previously existing file if a new featType
    #  was processed)
    savePKL(out_strandedFeats, flatFeats)
    
    return flatFeats
    

def normalizeLabels(df, out_strandedFeats, fasta, out_db = None,
                    exclusivePairs = pairs):
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
    fasta : (str)
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
    
    # First, if not done yet, flatten the relevant features
    flatFeats = getFlatFeatsByTypeStrdRef(out_strandedFeats, out_db,
                                          featsOfInterest)

    # Add up the cumulative feature length, summing up over strands & refs
    for featType, featsByStrdRef in flatFeats.items():
        # Go over strands
        for featsByRef in featsByStrdRef.values():
            # Go over refs
            for feats in featsByRef.values():
                regions[featType] += sum([feat[1] - feat[0] for feat in feats])
    
    # Get the genome length for this fasta file
    genLen = getGenomeLength(fasta)

    # Now calculate the proportion of genome covered by each label - e.g.,
    #  Intergenic = genome - genes
    #  Regulatory = genes - transcripts
    #  Intron = transcripts - exons
    #  Exon = exons

    # The first label in exclusivePairs depends on the length of the genome
    labelProps = np.array([(genLen*2 - regions[exclusivePairs[0][0]])
                           / (genLen*2)])
    # Add values for the subsequent labels in exclusivePairs
    for i in range(1,len(exclusivePairs)):
        labelProps = np.concatenate((labelProps, np.array(
            [(regions[exclusivePairs[i-1][0]] - regions[exclusivePairs[i][0]])
             / (genLen*2)])))
    # Add the last one, which is independent of other labels
    labelProps = np.concatenate((labelProps, np.array(
        [regions[exclusivePairs[-1][0]] / (genLen*2)])))
    
    # Normalize the df using these proportions
    df_norm = pd.DataFrame(data = df.values / labelProps[np.newaxis, :],
                           columns = df.columns,
                           index = df.index)
    
    return df_norm


def getSNRsByGeneLen(SNRsByLenStrdRef, concordant = True, sortedSNRs = True):
    """Function to sort SNRs by all genes that each of them maps to.

    Parameters
    ----------
    SNRsByLenStrdRef : (dict)
        { length : { strd : { ref : [ SNRs ] } } } or, alternatively when
        sortedSNRs = False: { length : [ SNRs ] }
    concordant : (bool), optional
        Indicates whether the sorting is by conrordant or discordant genes.
        The default is True.
    sortedSNRs : (bool), optional
        Indicates whether the SNRs are already sorted by length, strd, and ref.
        The default is True.

    Returns
    -------
    SNRsByGeneLen : (dict)
        { gene : { SNRlength : [ SNR ] } }
    """
    
    # Initialize a nested defaultdict
    SNRsByGeneLen = collections.defaultdict(
        lambda: collections.defaultdict(list))
    
    if sortedSNRs:
        for length, SNRsByStrdRef  in SNRsByLenStrdRef.items():
            for SNRsByRef in SNRsByStrdRef.values():
                for SNRs in SNRsByRef.values():
                    for snr in SNRs:
                        genesOfInterest = list(snr.concGenes.keys()) \
                            if concordant else list(snr.discGenes.keys())
                        for gene in genesOfInterest:
                            SNRsByGeneLen[gene][length].append(snr)
    else:
        # Go over all SNRs and sort them by concordant or discordant genes
        for length, SNRs in SNRsByLenStrdRef.items():
            for snr in SNRs:
                genesOfInterest = list(snr.concGenes.keys()) if concordant \
                    else list(snr.discGenes.keys())
                for gene in genesOfInterest:
                    SNRsByGeneLen[gene][length].append(snr)
    
    return SNRsByGeneLen



