#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 08:42:02 2020

@author: marek
"""

import time
import collections
from IPython.display import clear_output

# Global variable to be imported
resultDict = collections.defaultdict(list)
total = 0
processed = 0

class SNR:
    """An object that stores information about a Single Nucleotide Repeat
    sequence found.
    Notes:
        - Information wrt/ the reference genome.
        - Coordinates are 0-based.
    """
    def __init__(self, base, seq, chrom, start, end, mism, strand, feat):
        self.base = base                # str, 'A', 'T', 'C' or 'G'
        self.seq = seq                  # str, the sequence of the SNR
        self.chrom = chrom              # str, chromosome where SNR is found
        self.start = start              # int, start of the Site
        self.end = end                  # int, end of the Site
        self.mism = mism                # int, number of mismatches 
        self.strand = strand        # bool: T (+) or F (-)
        self.feat = feat                # set, GFFdb features for the
                                        #  respective strand, such as region,
                                        #  gene, exon, mRNA, CDS, etc.

    def __str__(self):
        # Displays the base, length, number of mismatches, chromosome, start,
        #  end wrt/ reference genome
        return 'poly({}): {} @ ch{}:{}-{}'.format(
            self.base,
            self.seq,
            self.chrom,
            self.start,
            self.end
            )



def splitter(base, refseq, minlength = 1000000):
    """
    Parameters
    ----------
    base : (str)
        Base that is being sought for by findSNRs.
    refseq : (str)
        Reference sequence to be divided into splits.
    minlength : (int), optional
        The minimum length for a split. The default is 1000000.

    Returns
    -------
    splits : (list)
        [ (first, last) ].
    """
    
    # Initiate the list and the first split
    splits = []; first = 0; last = minlength
    
    # Make a dictionary of letters to avoid for splitting
    avoidSplits = {
        'A':('a', 'A', 't', 'T'),
        'T':('a', 'A', 't', 'T'),
        'C':('c', 'C', 'g', 'G'),
        'G':('c', 'C', 'g', 'G'),
        }
    avoid = avoidSplits[base]
    
    # Generate valid splits until the last split exceeds length of the sequence
    while last < len(refseq):
        if refseq[last] in avoid:
            last += 1
        else:
            splits.append((first, last))
            first = last
            last += minlength
    
    # Add the last, sub-minlength split    
    splits.append((first, len(refseq)))
    
    return splits
    

def findSNRs(base, refname, refseq, extent):
    """Function to scan a given reference sequence (range) & find unique SNRs
    of given base type and length with given maximum number of mismatches
    allowed, such that the SNR contains mostly the primary bases and starts &
    ends with the primary base
    
    Parameters
    ----------
    base : (str)
        DNA base constituting SNRs to be searched for: "A", "T", "C", "G", "N"
    refname : (str)
        The name of the reference sequence to be searched.
    refseq : (str)
        Reference sequence or its part to be searched.
    extent : (tuple), optional
        The range of the sequence to search.

    Returns
    -------
    lengthToSNRs : (dict)
        { length : [ SNRs ] }
    lengthToSNRcounts : (dict)
        { length : SNRcount }
    """
    
    # Define the complement base to be sought
    compDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    cBase = compDict[base]
    # Handle small letters in the sequence
    capsDict = {
        'A':('a', 'A'),
        'T':('t', 'T'),
        'C':('c', 'C'),
        'G':('g', 'G'),
        'N':'N'
        }
    
    # Initiate the dicts for SNRs
    lengthToSNRs = collections.defaultdict(list)
    #lengthToSNRcounts = collections.defaultdict(int)
    
    # Initialize trackers
    first = 0; truEnd = len(refseq)
    
    
    # Keep scanning the sequence until the true end is reached
    while first < truEnd:
        # Set the last base in the range currently tested
        last = first + 1
        # Test the base & its complement at the beginning of each range
        for b in (base, cBase):
            # If found, keep testing the subsequent bases for this base
            nonCaps = capsDict[b]
            if refseq[first] in nonCaps:
                # As long as last hasn't reached the end and the same base
                #  is found on the next position, increase the "last"
                #  index of the range confirmed
                while last < truEnd and refseq[last] in nonCaps:
                    last += 1
                # Save the SNR into the dict by its length & break
                lengthToSNRs[last - first].append(
                    SNR(
                        base,
                        refseq[first:last],
                        refname,
                        first + extent[0],
                        last + extent[0],
                        0,
                        b == base,
                        None
                        )
                    )
                break
        # Move the first base to that which was tested last
        first = last
        
    # Announce finishing this sequence
    #print('Finished sequence {}[{}:{}]'.format(refname, extent[0], extent[1]))
    
    return lengthToSNRs


def getFeatures(lengthToSNRs, db):
    """Function that takes a list generated by findSNRs, looks up the set of
    features for each SNR, and sorts the SNRs into a dictionary by their length
    
    Parameters
    ----------
    SNRs : (list)
        List of Single Nucleotide Repeats (SNRs).
    db : (gffutils.FeatureDB)
        Database of genomic features.

    Returns
    -------
    lengthToSNRs : (defaultdict)
        A dictionary of SNRs sorted by their length.
    """
    
    # Pre-calculate the length of the list
    length = 0
    for SNRs in lengthToSNRs.values(): length += len(SNRs)
    length *= 100
        
    # Initiate current progress
    current = 0; i = 0
    
    for SNRs in lengthToSNRs.values():
        for SNR in SNRs:
            # Track progres through the SNRs
            i += 1
            # Calculate the current progress
            progress = round(i/length, 2)
            # If the new progress is different from the current
            if progress != current:
                # Announce progress through the list
                print('Processing SNRs: {}%'.format(progress))
                # Clear output but wait until the next one appears
                clear_output(wait = True)
                # Save the progress as current
                current = progress
            # Get the set of features from the database
            SNR.feat = {
                ft.featuretype for ft in db.region(
                    region = (SNR.chrom, SNR.start, SNR.end),
                    strand = '+' if SNR.strand else '-'
                    )
                }
    
    return lengthToSNRs


def howLongSince(t_start):
    """Function that prints out the amount of time elapsed since t_start.
    
    Parameters
    ----------
    t_start : (time.time())
        The timestamp since which to count.

    Returns
    -------
    None.
    """
    
    # Save the current time
    t_end = time.time()
    t_diff = t_end - t_start
    print(
        "{}h:{}m:{}s elapsed)".format(
            int(t_diff//(60*60)),
            int((t_diff%(60*60))//60),
            int(t_diff%60)
            )
        )

    
def collectResult(result):
    """Callback function for pooled function mapping.
    
    Parameters
    ----------
    result : (dict)
        Result yielded by multiprocessing.pool.map_async()

    Returns
    -------
    None.
    """
    global processed
    
    # Go over the results in the dict and add them to the result dict
    for length, SNRs in result.items(): resultDict[length].extend(SNRs)
    # Add to the global variable
    processed += 1
    # Announce progress
    print('Processed split: {}/{}'.format(processed, total))
    # Clear output but wait until the next one appears
    clear_output(wait = True)