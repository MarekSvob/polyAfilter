#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 18:07:19 2020

@author: marek
"""

import collections
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from Bio import SeqIO
from IPython.display import display


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


def sumKnownBases(bases):
    """Helper function to sum all bases across each known base, regardless of
    caps.
    Note: 'bases' has to be a (defaultdict) because if non-caps bases are not
    defined, a simple (dict) would return an error.
    
    Parameters
    ----------
    bases : (defaultdict)
        { base : count }

    Returns
    -------
    (int)
        Sum of all known bases in the genome.
    """
    
    return bases['A'] + bases['a'] + bases['T'] + bases['t'] \
        + bases['C'] + bases['c'] + bases['G'] + bases['g']


def getBaseComp(fasta, showGC = True):
    """Function to obtain the base composition of a genome by scanning.
    Optionally, the GC% is shown as a fact-check on correct processing.

    Parameters
    ----------
    fasta : (str)
        File path of the fasta reference file.
    showGC : (bool), optional
        Switch for showing the GC% content of this genome. The default is True.

    Returns
    -------
    bases : (defaultdict)
        { base : count }
    """

    bases = collections.defaultdict(int)

    print('Scanning the genome for base composition...')
    # Make sure that the fasta file is open only temporarily
    with open(fasta, 'r') as genome:
        # Go over each base in each record in the genome and count each base
        for record in SeqIO.parse(genome, "fasta"):
            for base in record.seq:
                bases[base] += 1
    # Optionally, show the GC% content of this genome
    if showGC:
        gc = bases['C'] + bases['c'] + bases['G'] + bases['g']
        print(
            'Scanning finished, G/C content: {:.2%}'.format(
                round(gc / sumKnownBases(bases), 2)
                )
            )
    
    return bases


def countsCheck(base, lengthToSNRcounts, bases):
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
    bases : defaultdict
        { base : count }

    Returns
    -------
    None.
    """
    
    # Sum the number of bases contained in the SNR counts
    SNRbases = 0
    for k,v in lengthToSNRcounts.items():
        SNRbases += int(k) * int(v)
    
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

def SNRcountTable(base, lengthToSNRcounts, bases, show = True):
    """Obtain a df of SNR Length, Observed counts, and Observed/Expected
    ratios.
    Note: 'bases' has to be a (defaultdict) because if non-caps bases are not
    defined, a simple (dict) would return an error.

    Parameters
    ----------
    base : (str)
        The base in question
    lengthToSNRcounts : (dict)
        { SNRlength : SNRcount }
    bases : (defaultdict)
        { base : count }
    show : (bool), optional
        Switch for showing the df as a table. The default is True.

    Returns
    -------
    df : (pandas.dataframe)
        The table with calculated values.
    """
    
    # Calculate the frequency of the base in question and its complement
    pBase = (bases[capsDict[base][0]] \
             + bases[capsDict[base][1]]) / sumKnownBases(bases)
    pComp = (bases[capsDict[compDict[base]][0]] \
             + bases[capsDict[compDict[base]][1]]) / sumKnownBases(bases)
    
    # Obtain the table data as a list of 3 respective lists: SNRlength,
    #  Observed absolute # of SNRs, and O/E ratio
    polyLen = []
    obs = []
    oe = []
    for k,v in lengthToSNRcounts.items():
        polyLen.append(int(k))
        obs.append(int(v))
        oe.append(
            (int(v)/int(sumKnownBases(bases))) \
                / ( (1-pBase) * pBase**int(k) * (1-pBase) \
                   + (1-pComp) * pComp**int(k) * (1-pComp) )
            )
    table_data = [polyLen,obs,oe]
    # Construct the df with the desired column names & types
    df = pd.DataFrame(data = pd.DataFrame(data = table_data).T)
    df.columns = ('SNR Length','Observed', 'O/E')
    df['Observed'] = df['Observed'].astype(int)
    df['SNR Length'] = df['SNR Length'].astype(int)
    df.sort_values(by = ['SNR Length'])
    # Optionally, show the dataframe as a styled table
    if show:
        styler = df.style.format('{:.2E}', subset = ['O/E']).hide_index()
        display(styler)
    
    return df
    

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
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize = (15,15))
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

