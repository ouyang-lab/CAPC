========================
Software Pre-requisites
========================
1) Python 2.7
  - numpy (>=1.13.3)
  - pysam (>=0.8.4)
2) Snakemake
3) Juicer Tools
4) Samtools
5) Bedtools
6) GNU sort (with --parallel option)
7) BWA MEM (>=0.7.12)
8) PICARD

Instructions to setup pipeline
1. Point to the correct executables
Edit the Snakefile by pointing to the correct path of each program's
executable or JAR file listed above.

2. Naming convention of FASTQ
Place FASTQ files under the data/fastq folder and name them as
<example1>_P1.fastq.gz and <example1>_P2.fastq.gz

3. Edit the targets
Edit the Snakefile by changing the sample names in targets_pre and
targets_hic to the FASTQ names used, i.e. (<example1>_P1.fastq.gz).
We will change the value of size (filtering value of orientation analysis)
in targets_hic after running the first stage of the pipeline.

4. Once the above is specified, run the pipeline
The pipeline consists of two stages
    1) Preprocessing to get the PRE file and perform orientation
       analysis to get the minimum distance in which the proportions
       of RF, FR, RR and FF reads are equal (+/-1%)
    2) Filter PRE file based on the minimum distance supplied and
        convert PRE to HiC file

Stage 1:
First run the following command to preprocess the CAP-C dataset

snakemake \
    -j <no_of_jobs> \
    -s <Snakefile> \
    --latency-wait 120 \
    targets_pre


Stage 2:
Edit the Snakefile by changing the size parameter in rule targets_hic
to the value stated in the first line of the orietnation analysis results
file under analysis/orientation/<sample>/<sample>.txt
If you do not want to do this, keep the default parameter of 1000.

Next run the following bash command to convert the PRE file to a HiC file

snakemake \
    -j <no_of_jobs> \
    -s <Snakefile> \
    --latency-wait 120 \
    targets_hic

