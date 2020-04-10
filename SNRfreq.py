#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 08:42:02 2020

@author: marek
"""

import os
import time
import collections
import pickle
import gffutils
from Bio import SeqIO
from IPython.display import clear_output
from multiprocessing import Pool


# Global variables to be initiated
resultSNRs = collections.defaultdict(list)
resultSNRcounts = collections.defaultdict(int)
totalPieces = 0
processed = 0
# Dict of complementary letters to consider both
compDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
# Handling of small letters in the sequence
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
        self.strand = strand            # bool: T (+) or F (-)
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


def splitter(base, refseq, minlength = 5000000):
    """Function to find splits in a given DNA sequence such that contigs of
    the base in question are not interrupted
    
    Parameters
    ----------
    base : (str)
        Base that is being sought for by findSNRs.
    refseq : (str)
        Reference sequence to be divided into splits.
    minlength : (int), optional
        The minimum length for a split. The default is 5000000.

    Returns
    -------
    splits : (list)
        [ (first, last) ].
    """
    
    # Initiate the list and the first split
    splits = []; first = 0; last = minlength
    
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

def getSplits(base, fasta, minlength = 5000000, mincont = 5):
    """Function that cuts up the genome into smaller pieces for parallel
    processing.

    Parameters
    ----------
    base : (str)
        Base that is being sought for by findSNRs.
    fasta : (str)
        Location of the fasta file.
    minlength : TYPE, optional
        The minimum length for a splits. The default is 5000000.
    mincont : TYPE, optional
        The minimal length of the SNRs saved by findSNRs. The default is 5.

    Returns
    -------
    pieces : (list)
        [ pieces ]
    """
    
    global totalPieces
    
    # Initiate a dictionary for splitters
    pieces = []
    # Keep genome open only as long as necessary
    with open(fasta, "r") as genome:
        for ch in SeqIO.parse(genome, 'fasta'):
            splitters = splitter(base, str(ch.seq), minlength)
            totalPieces += len(splitters)
            for ext in splitters:
                pieces.append(
                    (base, ch.id, str(ch.seq)[ext[0]:ext[1]], ext, mincont)
                    )
    
    return pieces

    
def findSNRs(base, refname, refseq, extent, minlen = 5):
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
    minlen : (int), optional
        The minimal length of the SNRs saved. (Note that all SNRs will be stil
        counted.) The default is 5.

    Returns
    -------
    lengthToSNRs : (dict)
        { length : [ SNRs ] }
    lengthToSNRcounts : (dict)
        { length : SNRcount }
    """
    
    # Define the complement base to be sought
    cBase = compDict[base]
    
    # Initiate the dicts for SNRs
    lengthToSNRs = collections.defaultdict(list)
    lengthToSNRcounts = collections.defaultdict(int)
    
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
                # Calculate the length of the SNR
                truLen = last - first
                # Count the SNR by its length
                lengthToSNRcounts[truLen] += 1
                # If the size is sufficient, save the SNR into the dict by
                #  its length & break
                if truLen >= minlen:
                    lengthToSNRs[truLen].append(
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
    
    return lengthToSNRs, lengthToSNRcounts

def assignFeatures(lengthToSNRs, gff, db_out):
    """Function that takes a dict generated by findSNRs and looks up & saves
    the set of features for each SNR.
    
    Parameters
    ----------
    lengthToSNRs : (dict)
        { length : [ SNRs ] }
    gff : (str)
        Location of the gff annotation file.
    db_out : (str)
        Location of the database file.

    Returns
    -------
    None.
    """
    
    # Connect to database; if it does not exist, create it from the gff first
    if os.path.isfile(db_out):
        db = gffutils.FeatureDB(db_out, keep_order = True)
    else:
        db = gffutils.create_db(
            gff,
            dbfn = db_out,
            force = True,
            merge_strategy = 'create_unique',
            id_spec = ['ID', 'Name'],
            verbose = True
            )
    
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
            progress = round(i/length, 3)
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
    """Callback function for pooled function mapping. Takes only one input var.
    
    Parameters
    ----------
    result : (dict)
        Result yielded by multiprocessing.pool.map_async()

    Returns
    -------
    None.
    """
    global processed
    
    SNRs, counts = result
    # Go over the result dicts and add them to the respective dict
    for length, SNRs in SNRs.items(): resultSNRs[length].extend(SNRs)
    for length, count in counts.items(): resultSNRcounts[length] += count
    # Add to the global variable
    processed += 1
    # Announce progress
    print('Processed split: {}/{}'.format(processed, totalPieces))
    # Clear output but wait until the next one appears
    clear_output(wait = True)
    
def saveSNRcsv(loc, base, lengthToSNRcounts):
    """Saving the SNR count dictionary as a csv file.

    Parameters
    ----------
    loc : (str)
        Location of where the csv should be saved
    lengthToSNRcounts : (dict)
        { length : SNRcount }

    Returns
    -------
    None.
    """
    with open(loc, 'w') as f:
        # Write the header
        line = 'poly{}/{} Length,Count\n'.format(base, compDict[base])
        f.write(line)
        for key in sorted(lengthToSNRcounts.keys()):
            f.write('{},{}\n'.format(key, lengthToSNRcounts[key]))
            
def saveSNRpkl(loc, lengthToSNRs):
    """Saving the SNR dictionary as a pkl file for later use.

    Parameters
    ----------
    loc : (str)
        Location of where the csv should be saved
    lengthToSNRs : (dict)
        { length : [ SNR ] }

    Returns
    -------
    None.
    """
    with open(loc, 'wb') as f:
        pickle.dump(lengthToSNRs, f)

def loadSNRs(genome, base):
    """Loading the SNR dictionary as a pkl file for further processing.

    Parameters
    ----------
    genome : (str)
        DESCRIPTION.
    base : (str)
        Identifies the base originally sought for.

    Returns
    -------
    lengthToSNRs : (dict)
        { length : [ SNR ] }
    """
    with open('{}_{}/{}.pkl'.format(genome, base, compDict[base]), 'rb') as f:
        lengthToSNRs = pickle.load(f)
    
    return lengthToSNRs

def processPool(pieces, cpus):
    """Parallel processing of genome scanning.

    Parameters
    ----------
    pieces : (list)
        [ pieces ]
    cpus : (int)
        Number of cpus to be used in parallel.

    Returns
    -------
    None.
    """
    pool = Pool(processes = cpus)
    for p in pieces:
        pool.apply_async(findSNRs, args = p, callback = collectResult)
    pool.close()
    pool.join()