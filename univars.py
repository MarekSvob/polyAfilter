#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 16:07:57 2020

@author: marek

File containing the global variables for SNR analysis.
"""

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
