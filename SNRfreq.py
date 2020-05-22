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
from multiprocessing import Pool
from random import randint


# Global result, tracking, and help variables to be initiated
resultSNRs = collections.defaultdict(list)
resultSNRcounts = collections.defaultdict(int)
totalPieces = 0
processed = 0
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


class SNR:
    """An object that stores information about a Single Nucleotide Repeat
    sequence found.
    Notes:
        - Information wrt/ the reference genome.
        - Coordinates are 0-based.
    """    
    
    def __init__(self, base, record, start, end, mism, strand, feats, genes):
        self.base = base        # str, 'A', 'T', 'C' or 'G'
        self.record = record    # str, record (chromosome) where SNR is found
        self.start = start      # int, start of the SNR locus
        self.end = end          # int, end of the SNR locus
        self.mism = mism        # int, number of mismatches 
        self.strand = strand    # bool: T (+) or F (-)
        self.feats = feats      # set, GFFdb features on the respective strand,
                                #  such as region, gene, exon, transcript, etc.
        self.genes = genes      # dict, { gene : bool } genes (transcripts)
                                #  onto which SNR maps & whether it is an exon

    def __str__(self):
        # Displays the base, number of mismatches, and location
        return 'poly({})[-{}] @ {}:{:,}-{:,}'.format(
            self.base,
            self.mism,
            self.record,
            self.start,
            self.end,
            )


def getPieces(base, fasta, cpus, cFR):
    """Function that organizes the genome into pieces of lengths following a
    uniform distribution.

    Parameters
    ----------
    base : (str)
        The DNA base to comprise the SNRs.
    fasta : (str)
        Reference to the FASTA file with the reference sequence.
    cpus : (int)
        The number of concurrent processes to be used in downstream analysis.
    cFR : (tuple)
        The range of denominators to be used in deciding on the length of
        pieces, i.e.: (A,B) -> (genomeLen//(A*cpus), genomeLen//(B*cpus)).

    Returns
    -------
    allPieces : (list)
        List that contains all the pieces (lists), each piece containing
        multiple slices (tuples), such that [ [ (recordID, start, seq) ] ].
    """
    
    # Grab a global variable to be changed in the course of this function
    global totalPieces
    
    def addSlice():
        """Helper function to manage adding of individual slices to the
        current piece and piece to allPieces, if necessary.
        """
        # Nonlocal ('global' in the nearest scope) variables to be changed
        nonlocal piece, pieceLen, allPieces, unit
        
        # Add a slice to the current piece
        piece.append((ch.id, first, str(ch.seq)[first:last]))
        # Increase the size of the piece
        pieceLen += last - first
        # If the piece has exceeded the target length, add the piece to
        #  the master list, reset it, and set a new target length
        if pieceLen >= unit:
            allPieces.append(piece)
            piece = []
            pieceLen = 0
            unit = randint(*ran)
    
    print('Measuring the size of the genome...')
    # Set the bases to avoid for splitting
    avoid = avoidSplits[base]
    # Initiate the genome length to be measured
    genomeLength = 0
    # Keep genome open only as long as necessary, while going over each record
    #  to measure the total genome length
    with open(fasta, "r") as genome:
        for ch in SeqIO.parse(genome, 'fasta'):
            genomeLength += len(ch.seq)
    
    # Initiate the range of piece sizes as a reusable tuple (which is a range
    #  of pre-set franctions of the genome divided by the # of cpus) & announce
    ran = (genomeLength//(cFR[0]*cpus), genomeLength//(cFR[1]*cpus))
    print('The total genome length is {:,} bp'.format(genomeLength))
    print(
        'Each piece will contain between {:,} and {:,} bases...'.format(*ran)
        )
    # Initiate the parameters for the first piece & the master list itself
    unit = randint(*ran)    
    piece = []
    pieceLen = 0
    allPieces = []
    
    # Go over each record again and create the pieces by adding slices to them
    with open(fasta, "r") as genome:
        # For each record, initiate the first slice to be between 0 & up to
        #  whatever length remains to be added to the piece
        for ch in SeqIO.parse(genome, 'fasta'):     
            first = 0
            last = unit - pieceLen
            
            # Keep adding slices until 'last' has exceeded the record length
            while last < len(ch.seq):
                # If the current last base is one to be avoided, move on
                if ch.seq[last] in avoid:
                    last += 1
                # Otherwise add the piece & reset the frame of the next slice
                else:
                    addSlice() 
                    first = last
                    last += unit - pieceLen
            # Once the end of the record is reached, add the last piece and
            #  move on to the next record
            last = len(ch.seq)
            addSlice()
        # Once all the records have been scanned, add the last piece to the
        #  master list (unless just added, which is unlikely)
        if piece != []:
            allPieces.append(piece)
    
    # Save the total number of pieces to the respective global variable
    totalPieces = len(allPieces)
    # Sort the pieces in place by their total length, largest to smallest
    print(
        'Sorting {:,} pieces by the total number of bp...'.format(totalPieces)
        )
    allPieces.sort(key = lambda x: sum([len(i[2]) for i in x]), reverse = True)
    
    return allPieces

    
def findSNRs(base, piece, db_out, temp, mincont):
    """Function to scan a given reference sequence (range) & find unique SNRs
    of given base type and length with given maximum number of mismatches
    allowed, such that the SNR contains mostly the primary bases and starts &
    ends with the primary base
    
    Parameters
    ----------
    base : (str)
        DNA base constituting SNRs to be searched for: "A", "T", "C", "G", "N"
    piece : (list)
        A list of splits, such that: [ (record, start, sequence) ]
    db_out : (str)
        Address of the master copy of the database.
    temp : (str)
        Address of temporary space to create working copies of the database.
    mincont : (int)
        The minimal length of the SNRs saved. (Note that all SNRs will be still
        counted.)

    Returns
    -------
    lengthToSNRs : (dict)
        { length : [ SNRs ] }
    lengthToSNRcounts : (dict)
        { length : SNRcount }
    """

    # Define the complement base to be sought
    cBase = compDict[base]
    
    # Initiate the dicts for SNRs and their counts
    lengthToSNRs = collections.defaultdict(list)
    lengthToSNRcounts = collections.defaultdict(int)
    
    # Make a temporary copy of the db file for this thread to speed up the
    #  process (SQLite3 dbs do not truly support parallel inquiries)
    #  more info: https://github.com/mapbox/node-sqlite3/issues/408
    temp_f = tempfile.NamedTemporaryFile(suffix = '.db', dir = temp)
    shutil.copy(db_out, temp_f.name)
    
    # Open the db connection
    db_conn = gffutils.FeatureDB(temp_f.name, keep_order = True)
    
    # For each slice in the piece
    for s in piece:
        # Save the sequence to scan & initialize trackers
        refseq = s[2]
        first = 0
        end = len(refseq)
        # Keep scanning the sequence until the end is reached
        while first < end:
            # Initiate the last base of the potential new SNR
            last = first + 1
            # Test the base & its complement at the beginning of each new SNR
            for b in (base, cBase):
                # If found, keep testing the subsequent bases for this base in
                #  either caps or non-caps
                eitherCaps = capsDict[b]
                if refseq[first] in eitherCaps:
                    # As long as last hasn't reached the end of the record and
                    #  the same base is found on the subsequent position,
                    #  keep increasing the "last" index of the range confirmed
                    while last < end and refseq[last] in eitherCaps:
                        last += 1
                    # Calculate the length of the SNR when done
                    length = last - first
                    # Count the SNR by its length
                    lengthToSNRcounts[length] += 1
                    # If the size is sufficient, save the SNR into the dict by
                    #  its length with its attributes
                    if length >= mincont:
                        # Initialize the set of features & dict of genes
                        feats = set()
                        genes = collections.defaultdict(bool)
                        # For each feature from the db spanned by the SNR
                        for ft in db_conn.region(
                                region = (s[0], s[1] + first, s[1] + last + 1),
                                strand = '+' if b == base else '-'
                                ):
                            # Note that the region() query is 1-based (a,b)
                            #  open interval, as documented here:
                            #  https://github.com/daler/gffutils/issues/129
                            # Save the type of the feature & add to the set
                            feat = ft.featuretype
                            feats.update({feat})
                            # If this is a transcript, make sure that the gene
                            #  (a list of 1) it represents is in the dict
                            if feat in ('transcript', 'mRNA'):
                                # If gene_id is N/A, get 'gene' (DM genome)
                                [g] = ft.attributes['gene_id'] \
                                    if 'gene_id' in ft.attributes.keys() \
                                        else ft.attributes['gene']
                                genes[g]
                                # The gene is either added to the dict with
                                #  'False' or if it already exists, its bool
                                #  is not over-written
                            # If this is an exon, make sure that the gene it
                            #  represents is in the dict with val 'True'
                            elif feat == 'exon':
                                # If gene_id is N/A, get 'gene' (DM genome)
                                [g] = ft.attributes['gene_id'] \
                                    if 'gene_id' in ft.attributes.keys() \
                                        else ft.attributes['gene']
                                genes[g] = True
                        # Now add the new SNR to the appropriate dict
                        lengthToSNRs[length].append(
                            SNR(
                                base,
                                s[0],
                                s[1] + first,
                                s[1] + last,
                                0,
                                b == base,
                                feats,
                                genes
                                )
                            )
                    # For each SNR, break the for-loop so that if 'base' was
                    #  found in the first position of the SNR, 'cBase' is not
                    #  tested for in the same position again
                    break
            # Move over the frame for the new potential SNR
            first = last
    
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
    
    # Save how much has elapsed
    t_diff = time.time() - t_start
    print(
        "{}h:{}m:{}s elapsed".format(
            int(t_diff//(60*60)),
            int((t_diff%(60*60))//60),
            int(t_diff%60)
            )
        )

def collectResult(result):
    """Callback function for pooled function mapping. Can take only one input
    variable.
    
    Parameters
    ----------
    result : (dict)
        Result yielded by multiprocessing.pool.map_async(); in this case, it
        is a tuple such that: (SNRdict, countDict)

    Returns
    -------
    None.
    """
    
    # Get the global variables to be changed
    global processed, resultSNRs, resultSNRcounts
    
    # Unpack the result
    SNRdict, countDict = result
    # Go over the result dicts and add them to the respective global dicts
    for length, SNRs in SNRdict.items(): resultSNRs[length].extend(SNRs)
    for length, count in countDict.items(): resultSNRcounts[length] += count
    # Count this piece as processed
    processed += 1
    # Announce progress
    print('Processed split: {:,}/{:,}'.format(processed, totalPieces))
    
def saveSNRcsv(loc, lengthToSNRcounts):
    """Saving the SNR count dictionary as a csv file.

    Parameters
    ----------
    loc : (str)
        Location of where the csv should be saved.
    base : (str)
        The DNA base that has been sought.
    lengthToSNRcounts : (dict)
        { length : SNRcount }

    Returns
    -------
    None.
    """

    with open(loc, 'w') as f:
        for key in sorted(lengthToSNRcounts.keys()):
            f.write('{},{}\n'.format(key, lengthToSNRcounts[key]))
    
    print('SNR counts saved to {}'.format(loc))
            
def savePKL(loc, var):
    """Saving any variable as a pkl file for later use.

    Parameters
    ----------
    loc : (str)
        Location of where the pkl should be saved.
    var : (any)
        Any file

    Returns
    -------
    None.
    """

    with open(loc, 'wb') as f:
        pickle.dump(var, f)
        
    print('File saved to {}'.format(loc))
        
def loadSNRcsv(loc):
    """Loading the SNR count dictionary as a csv file for further processing.

    Parameters
    ----------
    loc : (str)
        Location of where the csv had been saved.

    Returns
    -------
    lengthToSNRcounts : (dict)
        { length : SNRcount }
    """
    
    lengthToSNRcounts = {}
    
    with open(loc, 'r') as f:
        for line in f:
            k,v = line.rstrip('\n').split(',')
            lengthToSNRcounts[int(k)] = int(v)
    print('SNR counts loaded from {}'.format(loc))
            
    return lengthToSNRcounts
    

def loadPKL(loc):
    """Loading any variable from a pkl file for further processing.
    Note: The approprite class needs to have been loaded into the environment.

    Parameters
    ----------
    loc : (str)
        Location of where the pkl had been saved.

    Returns
    -------
    var : (any)
        Any class variable
    """

    with open(loc, 'rb') as f:
        var = pickle.load(f)
    print('File loaded from {}'.format(loc))
    
    return var

def processPoolSNRs(
        base,
        fasta,
        gff,
        db_out,
        temp = '.',
        mincont = 5,
        cpus = 70,
        maxtasks = None,
        cFR = (40, 5)
        ):
    """Wrapper function for parallel processing of genome scanning for SNRs.

    Parameters
    ----------
    base : (str)
        The base whose SNRs to scan for.
    fasta : (str)
        The location of the FASTA reference sequence file.
    gff : (str)
        The location of the GFF reference annotation file.
    db_out : (str)
        The location of the reference annotation database file. If none exists
        at that location, one will be created.
    temp : (str)
        The location to create temporary copies of the reference annotation
        database file for each parallel thread. Requires ample space.
    mincont : (int)
        The smallest number of bases in an SNR for it to be saved with
        annotations (saves space & processing time). Note that all SNRs will
        be still counted.
    cpus : (int), optional
        Number of processes to be used in parallel. The default is 70.
    maxtasks : (int), optional
        Number of tasks each process does before it is removed & replaced. The
        default is None.
    cFR : (tuple)
        A parameter that determines the range of sizes of the pieces to be
        processed across processes. The default is (40, 5).

    Returns
    -------
    None.
    """

    # Create a database from the gff if it does not exist yet
    print('Checking the database...')
    if not os.path.isfile(db_out):
        print('None had been created - creating a new database...')
        gffutils.create_db(
            gff,
            dbfn = db_out,
            force = True,
            merge_strategy = 'create_unique',
            id_spec = ['ID', 'Name'],
            verbose = True
            )
    print('The database is ready.')
    # Create the list of pieces to be processed in parallel
    allPieces = getPieces(base, fasta, cpus, cFR)
    print('Looking for SNRs across {} parallel processes...'.format(cpus))
    # Create a pool of workers with given the # of processes & max tasks for
    #  each one before it is restarted
    pool = Pool(processes = cpus, maxtasksperchild = maxtasks)
    # Queue up the processes in the order of the list. The results are sent to
    #  the collectResult() function and from there in turn to the global
    #  result variables.
    for piece in allPieces:
        pool.apply_async(
            func = findSNRs,
            args = (base, piece, db_out, temp, mincont),
            callback = collectResult
            )
    # Close the pool
    pool.close()
    # Join the processes
    pool.join()