rule bwa:
  input:
    "data/genome.fa",
    "data/samples/{sample}.fastq"
  output:
    "mapped/{sample}.bam"
  conda: 
    "envs/mapping.yaml"
  shell:  
    """
    module load bwa/0.7.17
    module load samtools/1.9
    bwa mem -t {threads} {input} | samtools view -Sb - > {output}
    """


rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.9
        samtools sort -o {output} {input}
        """
        
        
        
