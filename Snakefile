# ===================
# Executables
# ===================
JAVA = "java"
SAMTOOLS = "samtools"
BEDTOOLS = "bedtools"
GNU_SORT = "~/software/coreutils-8.27/src/sort"
BWA = "bwa"
PYTHON = "python"
SEQPREP = "~/software/SeqPrep/SeqPrep"
CUTADAPT = "cutadapt"

# ===================
# Java JAR files
# ===================
PICARD = "~/software/picard-tools-2.2.2/picard.jar"
JUICER_TOOLS = "~/software/juicer/juicer_tools_0.7.0.jar"

# ===================
# Additional files
# ===================
BWA_INDEX = "/path/to/index"
FRAGS_BED = "data/metadata/mboi_mm10_ordered_example.bed"

# ===================
# General Parameters
# ===================
SAMPLE_ID = ["example1", "example2"]
GENOMIC_DISTANCE_TO_FILTER = 1000
QUALITY = ["1", "30"]

# ===================
# BL-CAP-C Parameters
# ===================
FORWARD_ADAPTER = "ACGCGATATCTTATCTGACT"
REVERSE_ADAPTER = "AGTCAGATAAGATATCGCGT"

READ_LENGTH = 150

# ===================
# Targets (BL-CAP-C)
# ===================
rule targets_split:
    input:
        expand("data/splits/{sample}/split.out",
               sample=SAMPLE_ID)

rule targets_zip:
    input:
        expand("data/splits/{sample}.fastq.gz",
               sample=SAMPLE_ID)

rule targets_prep:
    input:
        expand("data/fastq/{sample}_R1.fastq.gz",
               sample=SAMPLE_ID)


# ===================
# Targets (CAP-C)
# ===================
rule targets_pre:
    input:
        expand("data/contacts/mm10/pre/{sample}.txt.gz",
               sample=SAMPLE_ID)

rule targets_hic:
    input:
        expand("data/hic/filter/{size}/{sample}_MAPQ{quality}.hic",
               sample=SAMPLE_ID,
               size=GENOMIC_DISTANCE_TO_FILTER,
               quality=QUALITY)


# ===================
# Rules
# ===================
# ===================
# BL-CAP-C Specific
# ===================
rule split_files:
    input:
        P1 = "data/BLCAPC/{sample}_R1.fastq.gz",
        P2 = "data/BLCAPC/{sample}_R2.fastq.gz"
    output: "data/splits/{sample}/split.out"
    params:
        prefix = "data/splits/{sample}/{sample}"
    shell:
        """
        zcat {input.P1} | split -l 20000000 - {params.prefix}_R1_
        zcat {input.P2} | split -l 20000000 - {params.prefix}_R2_
        for i in {params.prefix}*; do
            mv ${{i}} ${{i}}.fastq
        done
        touch {output}
        """

rule gzip_splits:
    input:
        R1 = "data/splits/{sample}/{sample}_R1_{split}.fastq",
        R2 = "data/splits/{sample}/{sample}_R2_{split}.fastq"
    output:
        R1 = "data/splits/{sample}/{sample}_R1_{split}.fastq.gz",
        R2 = "data/splits/{sample}/{sample}_R2_{split}.fastq.gz"
    shell:
        """
        gzip {input.R1}
        gzip {input.R2}
        """

def get_all_splits(wildcards):
    splits, = glob_wildcards("data/splits/%s/{fastq}" % wildcards.sample)
    splits = [ "data/merge/%s/%s" % (wildcards.sample, split)
               for split in splits if split != "split.out"]
    return splits


rule combine_all:
    input: get_all_splits
    output:
        P1 = "data/fastq/{sample}_R1.fastq.gz",
        P2 = "data/fastq/{sample}_R2.fastq.gz",
        F1 = "failed/{sample}_R1.fastq.gz",
        F2 = "failed/{sample}_R2.fastq.gz",
        S1 = "data/stats/{sample}.txt"
    params:
        merge = "data/merge/{sample}",
        stats = "data/stats/{sample}"
    shell:
        """
        cat {params.merge}/*{{OK,UN,US}}_R1*.fastq.gz > {output.P1}
        cat {params.merge}/*{{OK,UN,US}}_R2*.fastq.gz > {output.P2}
        cat {params.merge}/*MN_R1*.fastq.gz > {output.F1}
        cat {params.merge}/*MN_R2*.fastq.gz > {output.F2}
        cat {params.stats}/*.txt > {output.S1}
        """

rule remove_bridge_linker_cutadapt:
    input:
        P1 = "data/splits/{sample}/{sample}_R1_{part}.fastq.gz",
        P2 = "data/splits/{sample}/{sample}_R2_{part}.fastq.gz"
    output:
        merge_T1 = "data/cutadapt/{sample}/merge_T1_{part}.fastq.gz",
        merge_T2 = "data/cutadapt/{sample}/merge_T2_{part}.fastq.gz",
        unmerge_T1 = "data/cutadapt/{sample}/unmerge_T1_{part}.fastq.gz",
        unmerge_T2 = "data/cutadapt/{sample}/unmerge_T2_{part}.fastq.gz",
        unmerge_A1 = "data/cutadapt/{sample}/unmerge_A1_{part}.fastq.gz",
        unmerge_G1 = "data/cutadapt/{sample}/unmerge_G1_{part}.fastq.gz",
        unmerge_A2 = "data/cutadapt/{sample}/unmerge_A2_{part}.fastq.gz",
        unmerge_G2 = "data/cutadapt/{sample}/unmerge_G2_{part}.fastq.gz",
    params:
        merge = "data/SeqPrep/{sample}/merge_P1_{part}.fastq.gz",
        unmerge_P1 = "data/SeqPrep/{sample}/unmerge_P1_{part}.fastq.gz",
        unmerge_P2 = "data/SeqPrep/{sample}/unmerge_P2_{part}.fastq.gz",
        forward = FORWARD_ADAPTER,
        reverse = REVERSE_ADAPTER,
        seqprep_dir = "data/SeqPrep/{sample}",
        cutadapt_dir = "data/cutadapt/{sample}",
        seqprep = SEQPREP,
        cutadapt = CUTADAPT
    shell:
        """
        mkdir -p {params.seqprep_dir}
        mkdir -p {params.cutadapt_dir}
        {params.seqprep} \
            -f {input.P1} \
            -r {input.P2} \
            -s {params.merge} \
            -1 {params.unmerge_P1} \
            -2 {params.unmerge_P2}
        {params.cutadapt} -n 1 --overlap 10 \
            -a forward={params.forward} \
            -a reverse={params.reverse} \
            -o {output.merge_T1} \
            {params.merge}
        {params.cutadapt} -n 1 --overlap 10 \
            -g forward={params.forward} \
            -g reverse={params.reverse} \
            -o {output.merge_T2} \
            {params.merge}
        {params.cutadapt} -n 1 --overlap 10 \
            -a forward={params.forward} \
            -a reverse={params.reverse} \
            -o {output.unmerge_T1} \
            {params.unmerge_P1}
        {params.cutadapt} -n 1 --overlap 10 \
            -g forward={params.forward} \
            -g reverse={params.reverse} \
            -o {output.unmerge_T2} \
            {params.unmerge_P1}
        {params.cutadapt} -n 1 --overlap 10 \
            -a forward={params.forward} \
            -a reverse={params.reverse} \
            -o {output.unmerge_A1} \
            --mask-adapter \
            {params.unmerge_P1}
        {params.cutadapt} -n 1 --overlap 10 \
            -g forward={params.forward} \
            -g reverse={params.reverse} \
            -o {output.unmerge_G1} \
            --mask-adapter \
            {params.unmerge_P1}
        {params.cutadapt} -n 1 --overlap 10 \
            -a forward={params.forward} \
            -a reverse={params.reverse} \
            -o {output.unmerge_A2} \
            --mask-adapter \
            {params.unmerge_P2}
        {params.cutadapt} -n 1 --overlap 10 \
            -g forward={params.forward} \
            -g reverse={params.reverse} \
            -o {output.unmerge_G2} \
            --mask-adapter \
            {params.unmerge_P2}
        """

rule remove_bridge_linker_modified:
    input:
        merge_T1 = "data/cutadapt/{sample}/merge_T1_{part}.fastq.gz",
        merge_T2 = "data/cutadapt/{sample}/merge_T2_{part}.fastq.gz",
        unmerge_G1 = "data/cutadapt/{sample}/unmerge_G1_{part}.fastq.gz",
        unmerge_A1 = "data/cutadapt/{sample}/unmerge_A1_{part}.fastq.gz",
        unmerge_G2 = "data/cutadapt/{sample}/unmerge_G2_{part}.fastq.gz",
        unmerge_A2 = "data/cutadapt/{sample}/unmerge_A2_{part}.fastq.gz"
    output:
        P1 = "data/merge/{sample}/{sample}_R1_{part}.fastq.gz",
        P2 = "data/merge/{sample}/{sample}_R2_{part}.fastq.gz",
        e = "data/stats/{sample}/{sample}_{part}.txt"
    params:
        prefix = "data/merge/{sample}/{sample}_{part}",
        rlen = READ_LENGTH,
        python = PYTHON
    shell:
        """
        paste <(zcat {input.unmerge_G1}) <(zcat {input.unmerge_A1}) \
              <(zcat {input.unmerge_G2}) <(zcat {input.unmerge_A2}) | \
            {params.python} src/filter_N.py \
                -s {params.rlen} -m1 {input.merge_T1} -m2 {input.merge_T2} \
                -prefix {params.prefix} -e {output.e}
        touch {output.P1}
        touch {output.P2}
        """


# ===================
# CAP-C Specific
# ===================
# ===============================================
# TXT to Pre
# ===============================================
rule pre_to_hic:
    input: "data/contacts/mm10/pre/{sample}.txt.gz"
    output: "data/hic/raw/{sample}_MAPQ{quality}.hic"
    params:
        quality = "{quality}",
        userpref = "/tmp/{sample}_{rep}_MAPQ{quality}",
        java = JAVA,
        juicer = JUICER_TOOLS
    shell:
        """
        mkdir -p {parans.userpref}
        {params.java} -Duser.homer={params.userpref} -jar {params.juicer} \
            pre \
            -r 500,1000,2000,5000,10000,25000,40000,50000,100000,250000,500000,1000000,2500000 \
            -q {params.quality} \
            {input} \
            {output} \
            mm10
        """

rule filter_pre_to_hic:
    input: "data/contacts/mm10/pre_filter/{size}/{sample}_{rep}.txt.gz"
    output: "data/hic/filter/{size}/{sample}_{rep}_MAPQ{quality}.hic"
    params:
        quality = "{quality}",
        userpref = "/tmp/{sample}_{rep}_MAPQ{quality}",
        java = JAVA,
        juicer = JUICER_TOOLS
    shell:
        """
        mkdir -p {params.userpref}
        {params.java} -Duser.home={params.userpref} -jar {params.juicer} \
            pre \
            -r 500,1000,2000,5000,10000,25000,40000,50000,100000,250000,500000,1000000,2500000 \
            -q {params.quality} \
            {input} \
            {output} \
            mm10
        """

rule filter_pre:
    input: "data/contacts/mm10/pre/{sample}.txt.gz"
    output: "data/contacts/mm10/pre_filter/{size}/{sample}.txt.gz"
    params:
        size = "{size}",
        prefix = "data/contacts/mm10/pre_filter/{size}"
    shell:
        """
        mkdir -p {params.prefix}
        zcat {input} \
            | awk '$5!=$9' \
            | awk '{{if($3==$7){{if($4-$8>{params.size}){{print $0}}else if($4-$8<-{params.size}){{print $0}}}}else{{print $0}}}}' | gzip -c > {output}
        """

# ===============================================
# BAM to TXT
# ===============================================
rule extract_reads:
    input:
        bam = "data/bam/{sample}_byname.bam"
    output:
        p1 = "data/contacts/{genome}/reads/{sample}_P1.txt.gz",
        p2 = "data/contacts/{genome}/reads/{sample}_P2.txt.gz"
    params:
        python = PYTHON
    shell:
        """
        {params.python} src/get_reads.py {input.bam} {output.p1}.tmp {output.p2}.tmp
        cat {output.p1}.tmp | gzip -c > {output.p1}
        cat {output.p2}.tmp | gzip -c > {output.p2}
        rm {output.p1}.tmp
        rm {output.p2}.tmp
        """

rule reads_to_frags:
    input:
        txt = "data/contacts/{genome}/reads/{sample}.txt.gz",
        frag = "data/metadata/mboi_mm10_ordered.bed"
    output:
        bed = "data/contacts/{genome}/frags/{sample}.txt.gz"
    threads: 4
    params:
        bedtools = BEDTOOLS,
        sort = GNU_SORT
    shell:
        """
        {params.bedtools} intersect \
            -a <(zcat {input.txt}) \
            -b {input.frag} -wb \
            | cut -f1,2,3,4,5,6,7,11 \
            | {params.sort} -k4,4 --parallel={threads} \
            | gzip -c > {output.bed}
        """

rule frags_to_pre:
    input:
        p1 = "data/contacts/{genome}/frags/{sample}_P1.txt.gz",
        p2 = "data/contacts/{genome}/frags/{sample}_P2.txt.gz"
    output: "data/contacts/{genome}/pre/{sample}.txt.gz"
    threads: 4
    params:
        orientation = "analysis/orientation/{sample}",
        txt = "analysis/orientation/{sample}/{sample}.txt",
        python = PYTHON,
        sort = GNU_SORT
    shell:
        """
        paste <(zcat {input.p1}) <(zcat {input.p2}) \
            | {params.python} src/frags_to_pre.py medium \
            | {params.sort} -k3,3d -k7,7d --parallel={threads} \
            | gzip -c > {output}
        mkdir -p {params.orientation}
        {params.python} src/measure_orientation.py {output} {params.orientation}
        {params.python} src/calc_orientation.py {params.orientation}/*.txt.gz \
            > {params.txt}
        """


# ===============================================
# Map, Combine, Remove Duplicates, Sort By Name
# ===============================================
rule align_with_bwa:
    input: "data/fastq/{sample}.fastq.gz"
    output: "data/bam/individual/{sample}.bam"
    params:
        prefix = BWA_INDEX,
        bwa = BWA,
        samtools = SAMTOOLS
    threads: 8
    shell:
        """
        {params.bwa} mem -M -t {threads} {params.prefix} {input} \
            | {params.samtools} view -bS - > {output}
        """

rule combine_bam:
    input:
        P1 = "data/bam/individual/{sample}_P1.bam",
        P2 = "data/bam/individual/{sample}_P2.bam"
    output:
        bam = "data/bam/{sample}.sorted.bam"
    threads: 4
    params:
        samtools = SAMTOOLS,
        java = JAVA,
        python = PYTHON
    shell:
        """
        cat <({params.samtools} view -H {input.P1}) \
            <(paste -d "\\n" \
                <({params.samtools} view -F0x100 {input.P1}) \
                <({params.samtools} view -F0x100 {input.P2}) \
                | {params.python} src/combine_bam.py) \
            | {params.samtools} view -bS - \
            > {output.bam}.tmp
        {params.samtools} sort -@ {threads} {output.bam}.tmp -o {output.bam}
        rm {output.bam}.tmp
        """

rule remove_duplicates:
    input:
        bam = "data/bam/{sample}.sorted.bam"
    output:
        bam = "data/bam/{sample}.dedup.bam",
        txt = "data/statistics/dedup/{sample}.dedup.txt"
    params:
        samtools = SAMTOOLS,
        picard = PICARD,
        java = JAVA
    shell:
        """
        {params.java} -jar {params.picard} MarkDuplicates \
            REMOVE_DUPLICATES=true \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.txt}
        {params.samtools} index {output.bam}
        """

rule sort_bam_byname:
    input:
        bam = "data/bam/{sample}.dedup.bam"
    output:
        bam = "data/bam/{sample}_byname.bam"
    threads: 4
    params:
        samtools = SAMTOOLS
    shell:
        """
        {params.samtools} sort -n -@ {threads} {input.bam} -o {output.bam}
        """

rule flagstat:
    input: "data/bam/{sample}_{rep}_byname.bam"
    output: "data/bam/stats/{sample}_{rep}.log"
    params:
        samtools = SAMTOOLS
    shell:
        """
        {params.samtools} flagstat {input} > {output} 2>&1
        """
