# polyAfilter: Removal of internally primed alignments

## Introduction
**polyAfilter** is a python-based command line tool to filter BAM files to remove sequencing read alignments that have likely resulted from internal priming in poly(dT) priming-based bulk and single cell RNA sequencing library preparation methods (QuantSeq, 10X, inDrop, Seq-Well, Drop-Seq, CEL-Seq2, etc.). More details available in our bioRxiv pre-print:

Svoboda, M., Frost, H. R. & Bosco, G. **Internal oligo(dT) priming in bulk and single cell RNA sequencing.** _bioRxiv_ 2021.09.24.461289 (2021). [doi:10.1101/2021.09.24.461289](https://doi.org/10.1101/2021.09.24.461289)

## Dependencies

At minimum, the following is required to run **polyAfilter** [tested version of each dependency is indicated]:
* python [3.7.6]
* [samtools](http://www.htslib.org/) [1.10]
* [pysam](https://pysam.readthedocs.io/en/latest/) [0.16.0.1]
* [gffutils](http://daler.github.io/gffutils/) [0.10.1]
* [biopython](https://biopython.org/) [1.74]
* [numpy](https://numpy.org/doc/stable/) [1.17.0]
* [dill](https://dill.readthedocs.io/en/latest/) [0.3.1.1]

## Installation

Clone this repository into your local directory using the following command:

```bash
git clone https://github.com/MarekSvob/polyAfilter.git
```

This will create a folder titled `polyAfilter` with the contents of this repository.
After this has been done once, you can run the following command from within this folder in order to download the latest updates:

```bash
git pull
```

## Running polyAfilter

While **polyAfilter** contains a wider set of tools useful for A-SNR analysis, the follwing commands (usually run in this order) are core to its functionality and can be run from the command line:

### createDB

First, it is necessary to create a gff/gtf reference (`GTF/GFF_FILE`)-based database (`DB_FILE`), which provides fast access to annotation features. The `createDB` function is a wrapper around [gffutils.create_db](https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html), run with a specific set of parameters to ensure compatibility of the output with polyAfilter's other functions. Only one `DB_FILE` needs to be created for each `GTF/GFF_FILE`; this process may take several hours. _Make sure to use the same gtf/gff file that was used to create the BAM file alignment you would like to filter._
```
python polyAfilter.py createDB --help
```
```
usage: polyAfilter.py createDB [-h] [-v] GTF/GFF_FILE DB_FILE

Wrapper around the create_db function from the gffutils package using predefined parameters to
ensure that the database is created properly for the puposes of the BAMfilter.

positional arguments:
  GTF/GFF_FILE   Location of the reference GTF/GFF file
  DB_FILE        Location of the output annotation database

optional arguments:
  -h, --help     show this help message and exit
  -v, --verbose  Extra messages are logged.
```

### createTRANS

Next, a cache file (`TRANS_FILE`) is created, which contains information about all the transcripts expressed in the BAM file to be filtered. Only one `TRANS_FILE` needs to be created for each `BAM_FILE` and there are two options to do so:
- Run `createTRANS` command before running `BAMfilter` to create the associated `TRANS_FILE`. In this scenario, you can later omit the optional `--out_db DB_FILE` argument when running the `BAMfilter`. This process (`createTRANS`) will take a few hours to run but creation of the `TRANS_FILE` will be then skipped when later running the `BAMfilter`. This approach is especially useful if you are planning to run the `BAMfilter` on the the same `BAM_FILE` multiple times (perhaps to try several different values of the `COVLEN`, `MINSNRLEN`, and `MISM` parameters).
- Run `BAMfilter` directly. If an associated `TRANS_FILE` has not yet been created, the `--out_db DB_FILE` argument has to be supplied and the `TRANS_FILE` will be created first. This will extend the run time of `BAMfilter` the first time it is run but the subsequent runs, once a `TRANS_FILE` exists, will skip this step (even if the `--out_db DB_FILE` argument is supplied).

```
python polyAfilter.py createTRANS --help
```
```
usage: polyAfilter.py createTRANS [-h] [-i] DB_FILE BAM_FILE TRANS_FILE

Creates a cache file with all transcripts that are expressed as per the associated BAM file.

positional arguments:
  DB_FILE        Location of the reference annotation database
  BAM_FILE       Location of the sorted and indexed BAM file
  TRANS_FILE     Location of the output cache file

optional arguments:
  -h, --help     show this help message and exit
  -i, --introns  Instructs to include intronic coverage
```

### BAMfilter
Filter a `BAM_FILE` to remove alignments that likely resulted from internal poly(dT) priming up to `COVLEN` bp downstream of these alignments, onto genome-encoded poly(A) sequences (A-Single Nucleotide Repeats, A-SNRs) at least `MINSNRLEN` A's long, with up to `-m MISM` mismatches (the default is `0`). [More details in the associated [manuscript](https://doi.org/10.1101/2021.09.24.461289).] If a `TRANS_FILE` for this `BAM_FILE` does not yet exist, an `--out_db DB_FILE` needs to be provided to create one (as discussed above).

```
python polyAfilter.py BAMfilter --help
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

---
_Copyright (c) The Trustees of Dartmouth College_
