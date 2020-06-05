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
resSNRs = collections.defaultdict(list)
resSNRcounts = collections.Counter()
resConcFeats = collections.defaultdict(collections.Counter)
resDiscFeats = collections.defaultdict(collections.Counter)
genomeLengths = collections.defaultdict(int)
totalPieces = 0
processed = 0
timeStart = None

# A dictionary of letters to avoid for splitting, given the base sought
avoidSplits = {
    'A':('a', 'A', 't', 'T'),
    'T':('a', 'A', 't', 'T'),
    'C':('c', 'C', 'g', 'G'),
    'G':('c', 'C', 'g', 'G'),
    }

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


class SNR:
    """An object that stores information about a Single Nucleotide Repeat
    sequence found.
    Notes:
        - Information wrt/ the reference genome.
        - Coordinates are 0-based.
    """    
    
    def __init__(
            self,
            base,
            record,
            start,
            end,
            mism,
            strand,
            cFeats,
            cGenes,
            dFeats,
            dGenes
            ):
        self.base = base        # str, 'A', 'T', 'C' or 'G'
        self.record = record    # str, record (chromosome) where SNR is found
        self.start = start      # int, start of the SNR locus
        self.end = end          # int, end of the SNR locus
        self.mism = mism        # int, number of mismatches 
        self.strand = strand    # bool: T (+) or F (-)
        self.concFeats = cFeats # set, GFFdb features on the respective strand,
        self.discFeats = dFeats #  such as region, gene, exon, transcript, etc.
        self.concGenes = cGenes # dict, { gene : bool } genes (transcripts)
        self.discGenes = dGenes #  onto which SNR maps on the respective strand
                                #  & whether it is in an exon wrt/ this gene

    def __str__(self):
        # Displays the base, number of mismatches, and location
        return 'poly({})[-{}] @ {}:{:,}-{:,}'.format(
            self.base if self.strand else compDict[self.base],
            self.mism,
            self.record,
            self.start,
            self.end
            )


def getGenomeLength(fasta):
    """Function to measure the size of the genome.

    Parameters
    ----------
    fasta : (str)
        Reference to the FASTA file with the reference sequence.

    Returns
    -------
    genomeLength : (int)
        The total length of the genome.
    """
    
    global genomeLengths
    
    # Check if the fasta in question has already been measured
    if genomeLengths[fasta] == 0:
        print('Measuring the size of the genome...')
        # Initiate the genome length to be measured
        genLen = 0
        # Keep genome open only as long as necessary, while going over each
        #  record to measure the total genome length.
        with open(fasta, 'r') as genome:
            for ch in SeqIO.parse(genome, 'fasta'):
                genLen += len(ch.seq)
        genomeLengths[fasta] = genLen
        print(
            'The total genome length is {:,} bp.'.format(genomeLengths[fasta])
            )
    
    return genomeLengths[fasta]


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
    
    # Grab global variable to be changed in the course of this function
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
    
    # Set the bases to avoid for splitting
    avoid = avoidSplits[base]
    # Get the genome length
    genLen = getGenomeLength(fasta)
    
    # Initiate the range of piece sizes as a reusable tuple (which is a range
    #  of pre-set franctions of the genome divided by the # of cpus) & announce
    ran = (genLen//(cFR[0]*cpus), genLen//(cFR[1]*cpus))
    print(
        'Each piece will contain between {:,} and {:,} bp...'.format(*ran)
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
        'Sorting {:,} pieces by the total number of bp...'.format(
            totalPieces
            )
        )
    allPieces.sort(key = lambda x: sum([len(i[2]) for i in x]), reverse = True)
    
    return allPieces

    
def findSNRs(base, piece, db_out, temp, minFeatLen, minSNRlen):
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
    minFeatLen : (int)
        The minimal length of the SNRs to save its features; limited to save
        processing time. (Note that all SNRs will be still counted.)
    minSNRlen : (int)
        The minimal length of the SNRs saved; limited to save storage space
        needed. (Note that all SNRs will be still counted.)

    Returns
    -------
    lengthToSNRs : (dict)
        { length : [ SNRs ] }
    lengthToSNRcounts : (counter)
        { length : SNRcount }
    lenToConcFeats : (dict)
        { length : { { feature } : count } } for concordant features.
    lenToDiscFeats : (dict)
        { length : { { feature } : count } } for discordant features.
    """

    # Define the complement base to be sought
    cBase = compDict[base]
    
    # Initiate the dicts for SNRs and their counts
    lenToSNRs = collections.defaultdict(list)
    lenToSNRcounts = collections.Counter()
    lenToConcFeats = collections.defaultdict(collections.Counter)
    lenToDiscFeats = collections.defaultdict(collections.Counter)
    
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
                    lenToSNRcounts.update({length : 1})
                    # If the size is sufficient, look up the features
                    if length >= minFeatLen:
                        # Set the order of strands to look up features, always
                        #  first concordant with the SNR, then discordant
                        strds = ('+', '-') if b == base else ('-', '+')
                        # Initialize the list [cFeats, cGenes, dFeats, dGenes]
                        elems = []
                        for strd in strds:
                            # Initialize the set of features & genes
                            feats = set()
                            genes = collections.defaultdict(bool)
                            # For each feature from the db spanned by the SNR
                            for ft in db_conn.region(
                                region = (s[0], s[1] + first, s[1] + last + 1),
                                strand = strd
                                ):
                            # Note that the region() query is a 1-based (a,b)
                            #  open interval, as documented here:
                            #  https://github.com/daler/gffutils/issues/129
                                # Save the feature type & add to the set
                                feat = ft.featuretype
                                feats.update({feat})
                                # Look up gene info only if eventually added
                                if length >= minSNRlen:
                                    # If this is a transcript, make sure that
                                    #  the gene (a list of 1) it represents is
                                    #  added to the dict
                                    if feat in ('transcript', 'mRNA'):
                                        # If gene_id is N/A, get 'gene' (DM)
                                        [g] = ft.attributes['gene_id'] if \
                                            'gene_id' in ft.attributes.keys() \
                                                else ft.attributes['gene']
                                        genes[g]
                                        # The gene is either added to the dict
                                        #  w/ 'False' or if it already exists,
                                        #  its bool is not over-written
                                    # If this is an exon, make sure that the
                                    #  gene it represents is in the dict with
                                    #  value 'True'
                                    elif feat == 'exon':
                                        # If gene_id is N/A, get 'gene' (DM)
                                        [g] = ft.attributes['gene_id'] if \
                                            'gene_id' in ft.attributes.keys() \
                                                else ft.attributes['gene']
                                        genes[g] = True                            
                            # Save the elements
                            elems.append(feats)
                            elems.append(genes)
                        # Count the feats in the appropriate counter
                        lenToConcFeats[length].update({frozenset(elems[0]):1})
                        lenToDiscFeats[length].update({frozenset(elems[2]):1})
                        # If the SNR is large enough, add the SNR to the dict
                        if length >= minSNRlen:
                            lenToSNRs[length].append(
                                SNR(
                                    base,
                                    s[0],
                                    s[1] + first,
                                    s[1] + last,
                                    0,
                                    b == base,
                                    *elems
                                    )
                                )
                        # For each SNR, break the for-loop so that if 'base'
                        #  was found in the first position of the SNR, 'cBase'
                        #  is not tested for in the same position again
                    break
            # Move over the frame for the new potential SNR
            first = last
    # Close (and hence also delete) the temporary file
    temp_f.close()
    
    return lenToSNRs, lenToSNRcounts, lenToConcFeats, lenToDiscFeats


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
    
    # Save how much has elapsed (in seconds)
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
        is a tuple such that:
        (lenToSNRs, lenToSNRcounts, lenToConcFeats, lenToDiscFeats)

    Returns
    -------
    None.
    """
    
    # Get the global variables to be changed
    global processed, resSNRs, resSNRcounts, resConcFeats, resDiscFeats
    
    # Unpack the result
    lenToSNRs, lenToSNRcounts, lenToConcFeats, lenToDiscFeats = result
    # Go over the result dicts and add them to the respective global dicts
    for length, SNRs in lenToSNRs.items(): resSNRs[length].extend(SNRs)
    resSNRcounts.update(lenToSNRcounts)
    for length, counts in lenToConcFeats.items():
        resConcFeats[length].update(counts)
    for length, counts in lenToDiscFeats.items():
        resDiscFeats[length].update(counts)
    # Count this piece as processed
    processed += 1
    # Announce progress
    t_since = time.time() - timeStart
    t_left = t_since/processed * (totalPieces - processed)
    print(
        'Processed {:,} / {:,} ({:.2%}) splits in {}h:{}m:{}s. '\
            'Time remaining: ~{}h:{}m:{}s'.format(
                processed,
                totalPieces,
                processed / totalPieces,
                int(t_since//(60*60)),
                int((t_since%(60*60))//60),
                int(t_since%60),
                int(t_left//(60*60)),
                int((t_left%(60*60))//60),
                int(t_left%60),
                )
        )
    
    
def saveSNRcsv(loc, lengthToSNRcounts):
    """Saving the SNR count dictionary as a csv file.

    Parameters
    ----------
    loc : (str)
        Location of where the csv should be saved.
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
        Location of the pkl file where the variable should be saved.
    var : (any)
        Any variable

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


def getSNRs(
        base,
        fasta,
        gff,
        out_db,
        out_snrs,
        out_csv,
        out_concf,
        out_discf,
        temp = '.',
        minFeatLen = 1,
        minSNRlen = 5,
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
    out_db : (str)
        The location of the reference annotation database file. If none exists
        at that location, one will be created.
    out_snrs : (str)
        The location of the previously created SNRs pkl file. If none exists
        at that location, one will be created.
    out_csv : (str)
        The location of the previously created SNRs count pkl. If none exists
        at that location, one will be created.
    out_concf : (str)
        The location of the saved concordant feats dictionary.
    out_discf : (str)
        The location of the saved discordant feats dictionary.
    temp : (str)
        The location to create temporary copies of the reference annotation
        database file for each parallel thread. Requires ample space.
    minFeatLen : (int)
        The smallest number of bases for an SNR to have its features saved
        (processing time). Note that all SNRs will be still counted.
    minSNRlen : (int)
        The smallest number of bases for an SNR to be saved (saves space). Note
        that all SNRs will be still counted.
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
    resSNRs : (dict)
        { length : [ SNRs ] }
    resSNRcounts : (counter)
        { length : SNRcount }
    resConcFeats : (dict)
        { length : { { feature } : count } } for concordant features.
    resDiscFeats : (dict)
        { length : { { feature } : count } } for discordant features.
    """
    
    global timeStart, resSNRs, resSNRcounts, resConcFeats, resDiscFeats
    
    # If available, just load the previously saved data
    if os.path.isfile(out_snrs) and os.path.isfile(out_csv) \
        and os.path.isfile(out_concf) and os.path.isfile(out_discf):
        resSNRs = loadPKL(out_snrs)
        resSNRcounts = loadSNRcsv(out_csv)
        resConcFeats = loadPKL(out_concf)
        resDiscFeats = loadPKL(out_discf)
        return resSNRs, resSNRcounts, resConcFeats, resDiscFeats
    
    # Otherwise create it by processing
    # Create a database from the gff if it does not exist yet
    print('Checking for a database...')
    if not os.path.isfile(out_db):
        print('None had been created - creating a new database...')
        gffutils.create_db(
            gff,
            dbfn = out_db,
            force = True,
            merge_strategy = 'create_unique',
            id_spec = ['ID', 'Name'],
            verbose = True
            )
    print('The database is ready.')
    # Create the list of pieces to be processed in parallel
    allPieces = getPieces(base, fasta, cpus, cFR)
    print('Looking for SNRs across {} parallel processes...'.format(cpus))
    # Start measuring time
    timeStart = time.time()
    # Create a pool of workers with given the # of processes & max tasks
    #  for each one before it is restarted
    pool = Pool(processes = cpus, maxtasksperchild = maxtasks)
    # Queue up the processes in the order of the list. The results are sent
    #  to the collectResult() function and from there in turn to the global
    #  result variables.
    for piece in allPieces:
        pool.apply_async(
            func = findSNRs,
            args = (base, piece, out_db, temp, minFeatLen, minSNRlen),
            callback = collectResult
            )
    # Close the pool
    pool.close()
    # Join the processes
    pool.join()
    # Save the newly processed data
    savePKL(out_snrs, resSNRs)
    saveSNRcsv(out_csv, resSNRcounts)
    savePKL(out_concf, resConcFeats)
    savePKL(out_discf, resDiscFeats)
    
    return resSNRs, resSNRcounts, resConcFeats, resDiscFeats
