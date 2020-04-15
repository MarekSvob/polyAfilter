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
import tempfile
import shutil
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
    def __init__(self, base, seq, rec, start, end, mism, strd, feats, genes):
        self.base = base        # str, 'A', 'T', 'C' or 'G'
        self.seq = seq          # str, the sequence of the SNR
        self.rec = rec          # str, genomic record where SNR is found
        self.start = start      # int, start of the Site
        self.end = end          # int, end of the Site
        self.mism = mism        # int, number of mismatches 
        self.strd = strd        # bool: T (+) or F (-)
        self.feats = feats      # set, GFFdb features for the respective strd,
                                #  such as region, gene, exon, mRNA, CDS, etc.
        self.genes = genes      # set, genes onto which SNR maps (or empty)

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


def splitter(base, refseq, minlength = 20000000):
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
    # Determine which bases to avoid using a global dictionary
    avoid = avoidSplits[base]
    
    # Generate valid splits until the last split exceeds length of the sequence
    while last < len(refseq):
        # If the following base is to be avoided, move on
        if refseq[last] in avoid:
            last += 1
        # Otherwise save the split & start a new one
        else:
            splits.append((first, last))
            first = last
            last += minlength
    
    # Add the last, sub-minlength split    
    splits.append((first, len(refseq)))
    
    return splits


def getPieces(
        base,
        fasta,
        db_out,
        temp = '.',
        minlength = 20000000,
        mincont = 5
        ):
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
    
    # Announce that the function is currently active.
    print('Obtaining splits of minimum length {}...'.format(minlength))
    # Initiate a dictionary for splitters
    pieces = []
    # Keep genome open only as long as necessary
    with open(fasta, "r") as genome:
        for ch in SeqIO.parse(genome, 'fasta'):
            # Obtain the splits for this record
            splits = splitter(base, str(ch.seq), minlength)
            for split in splits:
                # Create pieces such that each can serve as an input to
                #  findSNRs() & add them to the list of pieces
                pieces.append(
                    (
                        base,
                        ch.id,
                        str(ch.seq)[split[0]:split[1]],
                        split,
                        db_out,
                        temp,
                        mincont
                        )
                    )
    # Edit the global variable to reflect the total number of pieces
    totalPieces = len(pieces)
    
    print('Sorting {} pieces by length...'.format(totalPieces))
    # Sort the pieces by sequence length, starting with the largest
    pieces.sort(key = lambda x: len(x[2]), reverse = True)
    # Now sort the pieces in an alternating manner:
    # Initiate the new list and starting indices
    optim_pieces = []
    left = 0
    right = len(pieces) - 1
    # Only add more as long as the optim_pieces list is shorter than pieces
    while len(optim_pieces) < len(pieces):
        # Add the left one first & increase index
        optim_pieces.append(pieces[left])
        left += 1
        if len(optim_pieces) < len(pieces):
            # Then add the other one & increase index
            optim_pieces.append(pieces[right])
            right -= 1
    
    return optim_pieces

    
def findSNRs(base, refname, refseq, split, db_out, temp = '.', mincont = 5):
    """Function to scan a given reference sequence (range) & find unique SNRs
    of given base type and length with given maximum number of mismatches
    allowed, such that the SNR contains mostly the primary bases and starts &
    ends with the primary base
    
    Parameters
    ----------
    base : (str)
        DNA base constituting SNRs to be searched for: "A", "T", "C", "G", "N"
    refname : (str)
        The name of the reference sequence searched.
    refseq : (str)
        Reference sequence or its part to be searched.
    split : (tuple), optional
        The range of the refseq that this specific split corresponds to.
    db_out : (str)
        Address of the master copy of the database.
    temp : (str), optional
        Address of temporary space to create working copies of the database.
    mincont : (int), optional
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
    
    # Make a temporary copy of the db file
    temp_f = tempfile.NamedTemporaryFile(suffix = '.db', dir = temp)
    shutil.copy(db_out, temp_f.name)
    
    # Open the db connection (for each process separately)
    db_conn = gffutils.FeatureDB(temp_f.name, keep_order = True)
    
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
                if truLen >= mincont:
                    feats = set()
                    genes = set()
                    for ft in db_conn.region(
                            region = (refname, split[0]+first, split[0]+last),
                            strand = '+' if b == base else '-'
                            ):
                        feat = ft.featuretype
                        feats.update(feat)
                        if feat in ('mRNA', 'transcript'):
                            genes.update(set(ft.attributes['gene']))
                    lengthToSNRs[truLen].append(
                        SNR(
                            base,
                            refseq[first:last],
                            refname,
                            split[0]+first,
                            split[0]+last,
                            0,
                            b == base,
                            feats,
                            genes
                            )
                        )
                break
        # Move the first base to that which was tested last
        first = last
        
    # Announce finishing this sequence
    #print('Finished sequence {}[{}:{}]'.format(refname, extent[0], extent[1]))
    
    return lengthToSNRs, lengthToSNRcounts


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
        "{}h:{}m:{}s elapsed".format(
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
    
    SNRdict, countDict = result
    # Go over the result dicts and add them to the respective dict
    for length, SNRs in SNRdict.items(): resultSNRs[length].extend(SNRs)
    for length, count in countDict.items(): resultSNRcounts[length] += count
    # Add to the global variable
    processed += 1
    # Clear output but wait until the next one appears
    clear_output(wait = True)
    # Announce progress
    print('Processed split: {}/{}'.format(processed, totalPieces))
    
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

def processPool(
        base,
        fasta,
        gff,
        db_out,
        temp = '.',
        minlength = 20000000,
        mincont = 5,
        cpus = 50,
        maxtasks = None
        ):
    """Parallel processing of genome scanning.

    Parameters
    ----------
    pieces : (list)
        [ pieces ]
    cpus : (int), optional
        Number of processes to be used in parallel. The default is 20.
    maxtasks : (int), optional
        Number of tasks each process does before it is removed & replaced. The
        default is None.

    Returns
    -------
    None.
    """
    # Create a database from the gff if it does not exist yet
    print('Checking database...')
    if not os.path.isfile(db_out):
        gffutils.create_db(
            gff,
            dbfn = db_out,
            force = True,
            merge_strategy = 'create_unique',
            id_spec = ['ID', 'Name'],
            verbose = True
            )
    # Create the list of pieces to be processed in parallel
    pieces = getPieces(
        base,
        fasta,
        db_out,
        temp = temp,
        minlength = minlength,
        mincont = mincont
        )
    print('Looking for SNRs in the splits across {} parallel processes...'.format(cpus))
    # Create a pool of workers with given # of processes & max tasks for each
    #  before it is replaced
    pool = Pool(processes = cpus, maxtasksperchild = maxtasks)
    # Queue up the processes
    for p in pieces:
        pool.apply_async(findSNRs, args = p, callback = collectResult)
    # Close the pool
    pool.close()
    # Join the processes
    pool.join()