#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 18:07:19 2020

@author: marek
"""

import collections
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Dict of complementary bases
compDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
# Handling of non-caps in the sequence
capsDict = {
    'A':('a', 'A'),
    'T':('t', 'T'),
    'C':('c', 'C'),
    'G':('g', 'G'),
    'N':'N'
    }
# Default exclusivePairs
pairs = [
    ('gene', 'Intergenic'),
    ('transcript', 'Regulatory'),
    ('exon', 'Intron')
    ]
   

def SNRcountOEplot(
        dfs,
        names,
        fontSize = 20,
        figSize = (15, 15),
        xlim = None,
        ylim = None
        ):
    """Plotter of the SNR dfs provided,

    Parameters
    ----------
    dfs : (list)
        [ df from SNRcountTable ]
    names : (list)
        [ series name ]
    fontSize : (int), optional
        Size of the plot's font. The default is 20.
    figSize : (tuple), optional
        (width, height) size in inches. The default is (15,15).
    xlim : (tuple)
        (left, right) limits of the x-axis. The defaults is None.
    ylim : (tuple)
        (bottom, top) limits of the y-axis. The defaults is None.

    Returns
    -------
    None.
    """
    
    # Set the figure parameters
    plt.rcParams.update({'font.size': fontSize})
    plt.figure(figsize = figSize)
    # Set the properties of each series
    for i in range(len(dfs)):
        plt.plot(
            dfs[i]['SNR Length'],
            np.log10(dfs[i]['O/E']),
            marker = 'o',
            linestyle = '-',
            label = names[i]
            )
    # Set the axes
    plt.xlabel('poly(A) length')
    plt.ylabel('Log10 O/E')
    plt.xlim(xlim)
    plt.ylim(ylim)
    # Add a gray line @0 (where O/E = 1)
    plt.axhline(y = 0, color = 'gray', linestyle = '--')
    # Place the legend
    plt.legend(loc = 'upper left')
    plt.show()


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
    
