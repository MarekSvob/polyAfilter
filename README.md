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

While **polyAfilter** contains a wider set of tools useful for A-SNR analysis, the following commands (usually run in this order) are core to its functionality and can be run from the command line:

### createDB

First, it is necessary to create a gff/gtf reference (`GTF/GFF_FILE`)-based database (`DB_FILE`), which provides fast access to annotation features. The `createDB` function is a wrapper around [gffutils.create_db](https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html), run with a specific set of parameters to ensure compatibility of the output with polyAfilter's other functions. Only one `DB_FILE` needs to be created for each `GTF/GFF_FILE`; this process may take several hours.

_Make sure to use the same `GTF/GFF_FILE` that was originally used to create the `BAM_FILE` alignment you would like to filter._
```bash
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
  -v, --verbose  Extra messages are logged. (default: False)
```

### createTRANS

Next, a cache file (`TRANS_FILE`) is created, which contains information about all the transcripts expressed in the BAM file to be filtered. Only one `TRANS_FILE` needs to be created for each `BAM_FILE` and there are two options to do so:
- Run `createTRANS` command before running `BAMfilter` to create the associated `TRANS_FILE`. In this scenario, you can later omit the optional `--out_db DB_FILE` argument when running the `BAMfilter`. This process (`createTRANS`) will take a few hours to run but creation of the `TRANS_FILE` will be then skipped when later running the `BAMfilter`. This approach is especially useful if you are planning to run the `BAMfilter` on the the same `BAM_FILE` multiple times (perhaps to try several different values of the `COVLEN`, `MINSNRLEN`, and `MISM` parameters).
- Run `BAMfilter` directly. If an associated `TRANS_FILE` has not yet been created, the `--out_db DB_FILE` argument has to be supplied and the `TRANS_FILE` will be created first. This will extend the run time of `BAMfilter` the first time it is run but the subsequent runs, once a `TRANS_FILE` exists, will skip this step (even if the `--out_db DB_FILE` argument is supplied).

_Make sure to use the `DB_FILE` created from the same `GTF/GFF_FILE` that was originally used to create the `BAM_FILE` alignment you would like to filter._

```bash
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
  -i, --introns  Instructs to include intronic coverage (default: False)
```

### BAMfilter
Filter a `BAM_FILE` to remove alignments that likely resulted from internal poly(dT) priming up to `COVLEN` bp downstream of these alignments, onto genome-encoded poly(A) sequences (A-Single Nucleotide Repeats, A-SNRs) at least `MINSNRLEN` A's long, with up to `-m MISM` mismatches (the default is `0`). [More details in the associated [manuscript](https://doi.org/10.1101/2021.09.24.461289).] If a `TRANS_FILE` for this `BAM_FILE` does not yet exist, an `--out_db DB_FILE` needs to be provided to create one (as discussed [above](https://github.com/MarekSvob/polyAfilter#createtrans)).

_Make sure to use the `FASTA_FILE` that was originally used to create the `BAM_FILE` alignment you would like to filter. Similarly, only use the `TRANS_FILE` associated with this specific `BAM_FILE`._

```bash
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
                        has already been created. (default: None)
  -b N, --base N        DNA base constituting SNRs to be searched for (default: A)
  -m MISM, --mismatches MISM
                        The maximum number of mismatches allowed in an SNR (default: 0)
  -o OUT_BAM_FILE, --out_bamfile OUT_BAM_FILE
                        Location of the resulting filtered bamfile; if not provided, the name of the
                        input bamfile name is used with ".filtered.bam" appended. (default: None)
  -i, --introns         Instructs to include intronic coverage (default: False)
  -c CB_FILE, --cbFile CB_FILE
                        If the scumi package is being used for alignment counting, provide the
                        location of the file with cell barcode count (output of scumi merge_fastq)
                        (default: None)
  -e OUT_CB_FILE, --out_cbFile OUT_CB_FILE
                        If the scumi package is being used for alignment counting, provide the
                        location for the modified cell barcode count file to be used as input for
                        scumi count_umi; if not provided in the presence of cbFile, cbFile name is
                        used with ".filtered.tsv" appended (default: None)
  -v, --verbose         Extra messages are logged. (default: False)
  -p P, --processes P   Set the maximum number of processes to use besides the main process. If not
                        provided, serial processing is used. Note that parallel processing will
                        produce the following warning for each temp file, which can be safely
                        ignored: "[E::idx_find_and_load] Could not retrieve index file for <file>"
                        (default: None)
```

### Example
The following code shows an example of filtering the 10X [BAM file](https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/500_PBMC_3p_LT_Chromium_X_possorted_genome_bam.bam) (with an associated [BAM index](https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/500_PBMC_3p_LT_Chromium_X_possorted_genome_bam.bam.bai)) from the "500 Human PBMCs, 3' LT v3.1, Chromium X" dataset. According to the associated [web summary](https://cf.10xgenomics.com/samples/cell-exp/6.1.0/500_PBMC_3p_LT_Chromium_X/500_PBMC_3p_LT_Chromium_X_web_summary.html), this BAM file was created using the [refdata-gex-GRCh38-2020-A](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz) 10X Genomics reference.

The code below assumes that the BAM file, the BAM index file, and the reference folder are located inside the `TEST_DIR`. The BAM file is filtered to remove alignments up to `300` (`COVLEN`) bps upstream of all A-SNRs at least `10` (`MINSNRLEN`) A's long, with up to `1` (`-m MISM`) mismatch:

```bash
# DIRECTORIES [edit these variables with your paths]
WORKING_DIR=/path/to/polyAfilter
TEST_DIR=/path/to/test/folder/

# INPUT FILES
FASTA=$TEST_DIR"refdata-gex-GRCh38-2020-A/fasta/genome.fa"
GTF=$TEST_DIR"refdata-gex-GRCh38-2020-A/genes/genes.gtf"
BAM=$TEST_DIR"500_PBMC_3p_LT_Chromium_X_possorted_genome_bam.bam"

# OUTPUT FILES
DB=$TEST_DIR"gtfdb.db"
TRANS=$TEST_DIR"trans.pkl"

NTHREADS=4

# SCRIPT
cd $WORKING_DIR

# Create the GTF database
python polyAfilter.py createDB -v $GTF $DB

# Create the TRANS file
python polyAfilter.py createTRANS $DB $BAM $TRANS

# Filter the BAM file
python polyAfilter.py BAMfilter -m 1 -v -p $NTHREADS 300 10 $BAM $FASTA $TRANS
```
As a result of running this script, new `$DB`, `$TRANS` and `$BAM.filtered.bam` files should be created inside the `TEST_DIR`. The filtered BAM file (`$BAM.filtered.bam`) may be used for further downstream processing.

In case of 10X Genomics BAM files that are being processed using the [`cellranger`](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline (such as the one above), [`bamtofastq`](https://support.10xgenomics.com/docs/bamtofastq) may need to be used to convert the BAM file back to FASTQ files, which can be used as input to `cellranger count`. The following example specifically uses [`bamtofastq_linux`](https://github.com/10XGenomics/bamtofastq/releases/tag/v1.3.0), using the same variables defined in the script above:

```bash
# LOCATION OF BAMTOFASTQ [edit this variable with your path]
BAMTOFASTQ=/path/to/bamtofastq_linux

OUT=$TEST_DIR"bamToFastq"

# Convert filtered BAM file to FASTQ files
$BAMTOFASTQ --nthreads $NTHREADS $BAM".filtered.bam" $OUT

# Re-do the count analysis
module load cellranger

cellranger count \
  --localcores $NTHREADS \
  --fastqs $OUT"/500_PBMC_3p_LT_Chromium_X_0_1_HFFLJDSX2" \
  --id count \
  --transcriptome $TEST_DIR"refdata-gex-GRCh38-2020-A"
```

---
_Copyright (c) The Trustees of Dartmouth College_
