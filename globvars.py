#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 16:07:57 2020

@author: marek

File containing the global variables for SNR analysis.
"""

import collections


# Global result, tracking, and help variables to be initiated
resultSNRs = collections.defaultdict(list)
resultSNRcounts = collections.defaultdict(int)
totalPieces = 0
processed = 0
genomeLength = 0
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
# A dictionary of letters to avoid for splitting, given the base sought
avoidSplits = {
    'A':('a', 'A', 't', 'T'),
    'T':('a', 'A', 't', 'T'),
    'C':('c', 'C', 'g', 'G'),
    'G':('c', 'C', 'g', 'G'),
    }
# Default exclusivePairs
pairs = [
    ('gene', 'Intergenic'),
    ('transcript', 'Regulatory'),
    ('exon', 'Intron')
    ]