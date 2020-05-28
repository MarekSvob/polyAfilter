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
import numpy as np
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from random import randint


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


def measureGenome(fasta):
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
    # Initiate the genome length to be measured
    genomeLength = 0
    # Keep genome open only as long as necessary, while going over each record
    #  to measure the total genome length.
    with open(fasta, 'r') as genome:
        for ch in SeqIO.parse(genome, 'fasta'):
            genomeLength += len(ch.seq)
            
    return genomeLength


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
    global totalPieces, genomeLength
    
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
    # Get the genome length if not already done before
    if genomeLength == 0:
        genomeLength = measureGenome(fasta)
    
    # Initiate the range of piece sizes as a reusable tuple (which is a range
    #  of pre-set franctions of the genome divided by the # of cpus) & announce
    ran = (genomeLength//(cFR[0]*cpus), genomeLength//(cFR[1]*cpus))
    print('The total genome length is {:,} bp.'.format(genomeLength))
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
                                strand = '+' if b == base else '-'
                                ):
                            # Note that the region() query is a 1-based (a,b)
                            #  open interval, as documented here:
                            #  https://github.com/daler/gffutils/issues/129
                                # Save the feature type & add to the set
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
                                    # The gene is either added to the dict w/
                                    #  'False' or if it already exists, its
                                    #  bool is not over-written
                                # If this is an exon, make sure that the gene
                                #  it represents is in the dict with val 'True'
                                elif feat == 'exon':
                                    # If gene_id is N/A, get 'gene' (DM genome)
                                    [g] = ft.attributes['gene_id'] \
                                        if 'gene_id' in ft.attributes.keys() \
                                            else ft.attributes['gene']
                                    genes[g] = True
                            elems.append(feats)
                            elems.append(genes)
                        # Now add the new SNR to the dict (keeping itr 0-based)
                        lengthToSNRs[length].append(
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


def normalizeLabels(
        df,
        out_feats,
        out_db = None,
        fasta = None,
        exclusivePairs = pairs,
        other = 'Exon'
        ):
    """Function to calculate what proportion of the genome is covered by each
    of the respective labels and subsequenly normalize the df by these
    proportions. If the flattened feats have not been processed previously,
    they will be processed and saved (can take up to 5 hrs).
    
    Parameters
    ----------
    df : (dataframe)
        The df of measured SNR labels to be normalized.
    out_feats : (str)
        The location of the PKL file containing the dict of flattened feats.
    fasta : (str), optional
        The location of the FASTA reference sequence file. The default is None.
    out_db : (str), optional
        The location of the reference annotation database file. The default is
        None.
    exclusivePairs : (list), optional
        A list of tuples with features and the respective labels assigned to
        mark their absence, such that [ (feature, label) ]. The default is
        [('gene','Intergenic'),('transcript','Regulatory'),('exon','Intron')].

    Returns
    -------
    df_norm : (dataframe)
        Normalized dataframe
    """
    
    global genomeLength
    
    regions = collections.defaultdict(int)
    # First, if not done yet, flatten the relevant features, using the
    #  simplified peak-merging strategy from Biostats Final Project
    # Note that GFF features are 1-based, closed intervals
    if os.path.isfile(out_feats):
        flatFeats = loadPKL(out_feats)
        for k,feats in flatFeats.items():
            # 1-based, closed intervals [start, end]
            regions[k[:-1]] += sum([feat[2] - feat[1] + 1 for feat in feats])
    else:
        # Connect the db
        db_conn = gffutils.FeatureDB(out_db, keep_order = True)
        
        # Initiate the variables needed
        strands = ('+', '-')
        featsOfInterest = [p[0] for p in exclusivePairs]
        flatFeats = {}
        
        for featType in featsOfInterest:
            # Go over each strand separately
            for strd in strands:
                print('Going over {}{}s'.format(strd, featType))
                featList = []
                # Iterate through ALL features of this type, for each strand.
                for feat in db_conn.all_features(
                        featuretype = featType,
                        strand = strd
                        ):
                    # Save in a tuple: (ref, start, end)
                    nFeat = (feat.seqid, feat.start, feat.end)
                    # Go over the previous ones, check for overlap on the same
                    #  ref and save all the overlaps found in a list.
                    overlaps = []
                    for oFeat in featList:
                        # If the ref is the same and (
                        #  (the new Feat's beginning is within the old Feat)
                        #  or (the new Feat's end is within the oldFeat)
                        #  or (the newFeat spans the oldFeat)
                        #  ), add old Feat to the overlaps list
                        if nFeat[0] == oFeat[0] and (
                            (nFeat[1] >= oFeat[1] and nFeat[1] <= oFeat[2]) \
                            or (nFeat[2] >= oFeat[1] and nFeat[2] <= oFeat[2]) \
                            or (nFeat[1] < oFeat[1] and nFeat[2] > oFeat[2])
                            ):
                            overlaps.append(oFeat)
                    #  If overlaps have been found, merge with the new one
                    if overlaps != []:
                        # Initialize the start & end of the merged feature
                        #  using the new Feat (not in the overlaps list)
                        ref = nFeat[0]; start = nFeat[1]; end = nFeat[2]
                        # Go over all the overlaps & update the start & end of
                        #  the merged feature as necessary
                        for ft in overlaps:
                            if ft[1] < start:
                                start = ft[1]
                            if ft[2] > end:
                                end = ft[2]
                            # When done considering this feature, remove it
                            #  from the master list
                            featList.remove(ft)
                        # When done, add the new merged feature to the list
                        featList.append((ref, start, end))
                    # If no overlap has been found, add new as a new one
                    else:
                        featList.append(nFeat)
                # Once all features have been flattened, add up their lengths
                for feat in featList:
                    # 1-based, closed intervals [start, end]
                    regions[featType] += (feat[2] - feat[1] + 1)
                # Save the feats in the master list to be saved
                flatFeats['{}{}'.format(featType, strd)] = sorted(featList)
        # Save the flattened feats to speed up the process in the future
        savePKL(out_feats, flatFeats)
    
    # Measure the genome length, if not already done
    if genomeLength == 0:
        genomeLength = measureGenome(fasta)
    # Now calculate the proportion of genome covered by each label - e.g.,
    #  Intergenic = genome - genes
    #  Regulatory = genes - transcripts
    #  Intron = transcripts - exons
    #  Exon = exons

    # The first label in exclusivePairs depends on the length of the genome
    labelProps = np.array(
        [(genomeLength*2 - regions[exclusivePairs[0][0]]) / (genomeLength*2)]
        )
    # Add values for the subsequent labels in exclusivePairs
    for i in range(1,len(exclusivePairs)):
        labelProps = np.concatenate((labelProps, np.array([
            (regions[exclusivePairs[i-1][0]] - regions[exclusivePairs[i][0]]) \
                / (genomeLength*2)
            ])))
    # Add the last one, independent of other labels
    labelProps = np.concatenate(
        (labelProps, np.array([regions[exclusivePairs[-1][0]] / (genomeLength*2)]))
        )
    
    # Normalize the df using these proportions
    df_norm = pd.DataFrame(
        data = df.values / labelProps[np.newaxis, :],
        columns = df.columns,
        index = df.index
        )
    
    return df_norm


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
        sumKnownBases = bases['A'] + bases['a'] + bases['T'] + bases['t'] \
            + bases['C'] + bases['c'] + bases['G'] + bases['g']
        print(
            'Scanning finished, G/C content: {:.2%}'.format(
                round(gc / sumKnownBases, 2)
                )
            )
    
    return bases


def countsCheck(base, lengthToSNRcounts, out_bases, fasta = None):
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
    out_bases : (str)
        File path of the pre-calculated bases file.
    fasta : (str)
        File path of the fasta reference file. The default is None.

    Returns
    -------
    None.
    """
    
    # Sum the number of bases contained in the SNR counts
    SNRbases = 0
    for k,v in lengthToSNRcounts.items():
        SNRbases += int(k) * int(v)
    
    # Load the bases dictionary, if available
    if os.path.isfile(out_bases):
        bases = loadPKL(out_bases)
    # Otherwise, scan the genome to create it
    else:
        bases = getBaseComp(fasta)
        savePKL(out_bases, bases)
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


def SNRcountTable(base, lengthToSNRcounts, out_bases, fasta = None):
    """Obtain a df of SNR Length, Observed counts, and Observed/Expected
    ratios.
    Note: 'bases' has to be a (defaultdict) because if non-caps bases are not
    defined, a simple (dict) would return an error.

    Parameters
    ----------
    base : (str)
        The base in question.
    lengthToSNRcounts : (dict)
        { SNRlength : SNRcount }
    out_bases : (str)
        File path of the pre-calculated bases file.
    fasta : (str)
        File path of the fasta reference file. The default is None.

    Returns
    -------
    df : (pandas.dataframe)
        The table with calculated values.
    """
    
    # Load the bases dictionary, if available
    if os.path.isfile(out_bases):
        bases = loadPKL(out_bases)
    # Otherwise, scan the genome to create it
    else:
        bases = getBaseComp(fasta)
        savePKL(out_bases, bases)
    
    sumKnownBases = bases['A'] + bases['a'] + bases['T'] + bases['t'] \
        + bases['C'] + bases['c'] + bases['G'] + bases['g']
    
    # Calculate the frequency of the base in question and its complement
    pBase = (bases[capsDict[base][0]] \
             + bases[capsDict[base][1]]) / sumKnownBases
    pComp = (bases[capsDict[compDict[base]][0]] \
             + bases[capsDict[compDict[base]][1]]) / sumKnownBases
    
    # Obtain the table data as a list of 3 respective lists: SNRlength,
    #  Observed absolute # of SNRs, and O/E ratio
    polyLen = []
    obs = []
    oe = []
    for k,v in lengthToSNRcounts.items():
        polyLen.append(int(k))
        obs.append(int(v))
        oe.append(
            (int(v)/int(sumKnownBases)) \
                / ( (1-pBase) * pBase**int(k) * (1-pBase) \
                   + (1-pComp) * pComp**int(k) * (1-pComp) )
            )
    table_data = [polyLen,obs,oe]
    # Construct the df with the desired column names & types
    df = pd.DataFrame(data = pd.DataFrame(data = table_data).T)
    df.columns = ('SNR Length','Observed', 'O/E')
    df['Observed'] = df['Observed'].astype(int)
    df['SNR Length'] = df['SNR Length'].astype(int)
    df.sort_values(by = ['SNR Length'], inplace = True)
    
    return df