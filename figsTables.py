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
    
    return bases['A'] + bases['a'] + bases['T'] + bases['t'] \
        + bases['C'] + bases['c'] + bases['G'] + bases['g']

def getBaseComp(fasta):
    
    bases = collections.defaultdict(int)

    print('Scanning the genome for base composition...')
    with open(fasta, 'r') as genome:
        for record in SeqIO.parse(genome, "fasta"):
            for base in record.seq:
                bases[base] += 1
    
    gc = bases['C'] + bases['c'] + bases['G'] + bases['g']
    
    print(
        'Scanning finished, G/C content: {}%'.format(
            round(gc*100 / sumKnownBases(bases), 2)
            )
        )
    
    return bases


def countsCheck(base, lengthToSNRcounts, bases):
    
    #Check that the number of A/T bases is the same gotten by either method
    SNRbases = 0
    
    for k,v in lengthToSNRcounts.items():
        SNRbases += int(k) * int(v)
    
    scanned = bases[capsDict[base][0]] \
        + bases[capsDict[base][1]] \
            + bases[capsDict[compDict[base]][0]] \
                + bases[capsDict[compDict[base]][0]]
    
    if SNRbases == scanned:
        print(
            'The total number of {}/{} bases checks out!'.format(
                base,
                compDict[base]
                )
            )
    else:
        print(
            '{}/{} bases in SNRs: {}; from genome scanning {}'.format(
                base,
                compDict[base],
                SNRbases,
                scanned
                )
            )

def SNRcountTable(base, lengthToSNRcounts, bases):
    
    pBase = (bases[capsDict[base][0]] \
             + bases[capsDict[base][1]]) / sumKnownBases(bases)
    pComp = (bases[capsDict[compDict[base]][0]] \
             + bases[capsDict[compDict[base]][1]]) / sumKnownBases(bases)

    table_data = []
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
    df = pd.DataFrame(data = pd.DataFrame(data = table_data).T)
    df.columns = ('SNR Length','Observed', 'O/E')
    df['Observed'] = df['Observed'].astype(int)
    df['SNR Length'] = df['SNR Length'].astype(int)
    styler = df.style.format('{:.2E}', subset = ['O/E']).hide_index()
    display(styler)
    
    return df
    

def SNRcountPlot(df):
    
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize=(15,15))
    plt.plot(df['SNR Length'], np.log10(df['O/E']), '-o')
    plt.xlabel('poly(A) length')
    plt.ylabel('Log10 O/E')
    plt.axhline(y = 0, color = 'gray', linestyle='--')
    plt.show()

