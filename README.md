# Introduction
**polyAfilter** is a python-based command line tool to process BAM files to remove sequencing read alignments that have likely resulted from internal priming in poly(dT) priming-based bulk and single cell RNA sequencing library preparation methods (QuantSeq, 10X, inDrop, Seq-Well, Drop-Seq, CEL-Seq2, etc.).

More details available in the manuscript: [link]

### Dependencies

At minimum, the following is required to run **polyAfilter** [tested version of each dependency is indicated in brackets]
* python [3.7.6]
* [samtools](http://www.htslib.org/) [1.10]
* [pysam](https://pysam.readthedocs.io/en/latest/) [0.16.0.1]
* [gffutils](http://daler.github.io/gffutils/) [0.10.1]
* [biopython](https://biopython.org/) [1.74]
* [numpy](https://numpy.org/doc/stable/) [1.17.0]
* [dill](https://dill.readthedocs.io/en/latest/) [0.3.1.1]

* pandas
* iPython
* sklearn
* scumi
* regex

### Installation

Clone this repository into your local directory using the following command:

```bash
$ git clone https://github.com/MarekSvob/polyAfilter.git
```

This will create a folder titled `polyAfilter` with the contents of this repository.
After this has been done once, you can run the following command from within this folder in order to download the latest updates:

```bash
$ git pull
```

# Running polyAfilter

While **polyAfilter** contains a wider set of tools useful for A-SNR analysis, the follwing commands (usually run in this order) are core to its functionality:

```
$ python polyAfilter.py createDB --help
```
```
usage: polyAfilter.py createDB [-h] [-v] GTF/GFF_FILE OUT_GTF/GFF_FILE

Wrapper around the create_db function from the gffutils package using predefined parameters to
ensure that the database is created properly for the puposes of the BAMfilter.

positional arguments:
  GTF/GFF_FILE      Location of the reference GTF/GFF file
  OUT_GTF/GFF_FILE  Location of the resulting annotation database

optional arguments:
  -h, --help        show this help message and exit
  -v, --verbose     Extra messages are logged.
```

```
$ python polyAfilter.py createTRANS --help
```
```
usage: polyAfilter.py createTRANS [-h] [-i] TRANS_FILE DB_FILE BAM_FILE

Creates a cache file with all transcripts that are expressed as per the associated BAM file.

positional arguments:
  TRANS_FILE     Location of a cache file that stores expressed transcripts
  DB_FILE        Location of the reference annotation database
  BAM_FILE       Location of the sorted and indexed BAM file

optional arguments:
  -h, --help     show this help message and exit
  -i, --introns  Instructs to include intronic coverage
```

```
$ python polyAfilter.py BAMfilter --help
```
```
usage:  polyAfilter.py BAMfilter [-h] [-d DB_FILE] [-b N] [-m MISM] [-o OUT_BAM_FILE] [-i]
        [-c CB_FILE] [-e OUT_CB_FILE] [-v] [-p P] COVLEN MINSNRLEN BAM_FILE FASTA_FILE TRANS_FILE

Filters a sorted and indexed BAM file to remove sparse alignments that likely resulted from internal
priming.

positional arguments:
  COVLEN                Maximal coverage distance
  MINSNRLEN             Minimal SNR length
  BAM_FILE              Location of the sorted and indexed BAM file
  FASTA_FILE            Location of the reference FASTA file
  TRANS_FILE            Location of a cache file that stores expressed transcripts; if it does not
                        already exist, it will be created

optional arguments:
  -h, --help            show this help message and exit
  -d DB_FILE, --out_db DB_FILE
                        Location of the reference annotation database; not necessary if TRANS_FILE
                        has already been created.
  -b N, --base N        DNA base constituting SNRs to be searched for
  -m MISM, --mismatches MISM
                        The maximum number of mismatches allowed in an SNR
  -o OUT_BAM_FILE, --out_bamfile OUT_BAM_FILE
                        Location of the resulting filtered bamfile; if not provided, the name of the
                        input bamfile name is used with ".filtered.bam" appended.
  -i, --introns         Instructs to include intronic coverage
  -c CB_FILE, --cbFile CB_FILE
                        If the scumi package is being used for alignment counting, provide the
                        location of the file with cell barcode count (output of scumi merge_fastq)
  -e OUT_CB_FILE, --out_cbFile OUT_CB_FILE
                        If the scumi package is being used for alignment counting, provide the
                        location for the modified cell barcode count file to be used as input for
                        scumi count_umi; if not provided in the presence of cbFile, cbFile name is
                        used with ".filtered.tsv" appended
  -v, --verbose         Extra messages are logged.
  -p P, --processes P   Set the maximum number of processes to use besides the main process. If not
                        provided, serial processing is used. Note that parallel processing will
                        produce the following warning for each temp file, which can be safely
                        ignored: "[E::idx_find_and_load] Could not retrieve index file for <file>"
```

Copyright (c) The Trustees of Dartmouth College
