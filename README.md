# CAP-C and BL-CAP-C

The following instructions describe the pipeline used to process
all CAP-C and BL-CAP-C libraries in the manuscript. For BL-CAP-C
libraries, additional steps have been taken to merge overlapping
read pairs longer than fragments before removal of the bridge
linker adapter.

# Software Pre-requisites

## Dependencies and Prerequisites
The following software prerequisites are required to run the Snakemake pipeline
+ Python 2.7
+ numpy (>=1.13.3)
+ pysam (>=0.8.4)
+ Snakemake
+ Juicer Tools
+ Samtools
+ Bedtools
+ GNU sort (with --parallel option)
+ BWA MEM (>=0.7.12)
+ PICARD
+ SeqPrep (BL-CAP-C libraries)
+ cutadapt (BL-CAP-C libraries)

# General Instructions to setup pipeline

## 1. Point to the correct executables
Edit the Snakefile by pointing to the correct path of each program's
executable or JAR file listed above.

## 2. Naming convention of FASTQ
Place FASTQ files under the `data/fastq` folder and name them as
`<example1>_P1.fastq.gz` and `<example1>_P2.fastq.gz`

If the FASTQ files are from BL-CAP-C libraries, place them
under the `data/BLCAPC` folder instead. Naming conventions remain
the same as above.

## 3. Edit the targets
Edit the Snakefile by changing the sample names in targets_pre and
targets_hic to the FASTQ names used, i.e. `<example1>_P1.fastq.gz`.
We will change the value of size (filtering value of orientation analysis)
in targets_hic or `GENOMIC_DISTANCE_TO_FILTER` variable after running
the first stage of the pipeline.

## 4. Once the above is specified, run the pipeline
If you are processing BL-CAP-C libraries, run the BL-CAP-C pipeline first
before running the CAP-C pipeline described in the next subsection. If you are
processing CAP-C libraries, go straight to the CAP-C pipeline.

### The BL-CAP-C pipeline consists of three additional stages.
#### Stage 1. Split files into chunks of 5,000,000 reads per file

`snakemake -j <no_of_jobs> -s <Snakefile> --latency-wait 120 targets_split`


#### Stage 2. Gzip chunks

`snakemake -j <no_of_jobs> -s <Snakefile> --latency-wait 120 targets_zip`


#### Stage 3. Merge overlapping read pairs and remove bridge linkers

`snakemake -j <no_of_jobs> -s <Snakefile> --latency-wait 120 targets_prep`

Once done, continue onto CAP-C pipeline.

### The CAP-C pipeline consists of two stages:

#### Stage 1. Preprocessing to get the PRE file and perform orientation analysis to get the minimum distance in which the proportions of RF, FR, RR and FF reads are equal (+/-1%)

`snakemake -j <no_of_jobs> -s <Snakefile> --latency-wait 120 targets_pre`

#### Stage 2. Filter PRE file based on the minimum distance supplied and convert PRE to HiC file.

Edit the Snakefile by changing the size parameter or `GENOMIC_DISTANCE_TO_FILTER` variable
in rule targets_hic to the value stated in the first line of the orietnation analysis results
file under `analysis/orientation/\<sample>/\<sample>.txt`
If you do not want to do this, keep the default parameter of 1000.

Next run the following bash command to convert the PRE file to a HiC file

`snakemake -j <no_of_jobs> -s <Snakefile> --latency-wait 120 targets_hic`

