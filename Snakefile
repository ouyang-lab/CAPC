# ==================
# Executables
# ==================
JAVA = "java"
SAMTOOLS = "samtools"
BEDTOOLS = "bedtools"
GNU_SORT = "~/software/coreutils-8.27/src/sort"
BWA = "bwa"
PYTHON = "python"

# ==================
# Java JAR files
# ===================
PICARD = "~/software/picard-tools-2.2.2/picard.jar"
JUICER_TOOLS = "~/software/juicer/juicer_tools_0.7.0.jar"


# ===================
# Targets
# ===================
rule targets_pre:
    input:
        expand("data/contacts/mm10/pre/{sample}.txt.gz",
               sample=["example1", "example2"])

rule targets_hic:
    input:
        expand("data/hic/filter/{size}/{sample}_MAPQ{quality}.hic",
               sample=["example1", "example2"],
               size=["1000"],
               quality=["1", "30"])


# ===================
# Rules
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
        prefix = "/data/acheng/data/reference/mm10/bwa/mm10",
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
    shell:
        """
        {params.samtools} flagstat {input} > {output} 2>&1
        """
