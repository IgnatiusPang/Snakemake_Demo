rule bwa:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped/A.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.9
        bwa mem {input} | samtools view -Sb - > {output}
        """