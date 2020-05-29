#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 18:07:19 2020

@author: marek
"""

import os
import collections
import gffutils
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


def countsCheck(base, lengthToSNRcounts, out_bases, fasta = None):
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
    
    if SNRbases == scanned:
        print(
            'The total number of {}/{} bases ({:,}) checks out!'.format(
                base,
                compDict[base],
                SNRbases
                )
            )
    else:
        print(
            '{}/{} bases in SNRs: {:,}; from genome scanning {:,}'.format(
                base,
                compDict[base],
                SNRbases,
                scanned
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


def SNRfeatureSets(lengthToSNRs):
    """Function that compiles a bool df of Features vs. FeatureSets.
    
    Parameters
    ----------
    lengthToSNRs : (dict)
        The dict of SNRs detected, sorted by length.

    Returns
    -------
    df : (dataframe)
        A dataframe of features vs. featureSets.
    """
    
    # Get the set of (frozen)sets of features represented among the SNRs
    featureSets = set()
    for key,vals in lengthToSNRs.items():
        for SNR in vals:
            featureSets |= {frozenset(SNR.feats)}
    # Flatten this set of features represented
    features = {feat for featureSet in featureSets for feat in featureSet}
    
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


def SNRlabelProps(lengthToSNRs, exclusivePairs = pairs, other = 'Exon'):
    """Function that prepares a df of proportions of labeled SNRs for each
    length.    

    Parameters
    ----------
    lengthToSNRs : (dict)
        The dict of SNRs detected, sorted by length.
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
        (max(lengthToSNRs.keys()) + 1, len(exclusivePairs) + 1),
        dtype = int
        )
    # For each SNR, add +1 for each label @ the appropriate length
    for key, vals in lengthToSNRs.items():
        for SNR in vals:
            for i in range(len(exclusivePairs)):
                if exclusivePairs[i][0] not in SNR.feats:
                    index = i
                    break
            else:
                index = len(exclusivePairs)
            SNRlabels_data[key, index] += 1
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
    

def normalizeLabels(
        df,
        out_feats,
        out_db = None,
        fasta = None,
        exclusivePairs = pairs,
        other = 'Exon'
        ):
    """Function to calculate what proportion of the genome is covered by each
    of the respective labels and subsequenly normalize the df by these
    proportions. If the flattened feats have not been processed previously,
    they will be processed and saved (can take up to 5 hrs).
    
    Parameters
    ----------
    df : (dataframe)
        The df of measured SNR labels to be normalized.
    out_feats : (str)
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
    
    regions = collections.defaultdict(int)
    # First, if not done yet, flatten the relevant features, using the
    #  simplified peak-merging strategy from Biostats Final Project
    # Note that GFF features are 1-based, closed intervals
    if os.path.isfile(out_feats):
        flatFeats = loadPKL(out_feats)
        for k,feats in flatFeats.items():
            # 1-based, closed intervals [start, end]
            regions[k[:-1]] += sum([feat[2] - feat[1] + 1 for feat in feats])
    else:
        # Connect the db
        db_conn = gffutils.FeatureDB(out_db, keep_order = True)
        
        # Initiate the variables needed
        strands = ('+', '-')
        featsOfInterest = [p[0] for p in exclusivePairs]
        flatFeats = {}
        
        for featType in featsOfInterest:
            # Go over each strand separately
            for strd in strands:
                print('Going over {}{}s'.format(strd, featType))
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
                # Once all features have been flattened, add up their lengths
                for feat in featList:
                    # 1-based, closed intervals [start, end]
                    regions[featType] += (feat[2] - feat[1] + 1)
                # Save the feats in the master list to be saved
                flatFeats['{}{}'.format(featType, strd)] = sorted(featList)
        # Save the flattened feats to speed up the process in the future
        savePKL(out_feats, flatFeats)
    
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
    # Add the last one, independent of other labels
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