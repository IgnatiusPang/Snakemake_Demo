samples = ["A", "B", "C"] ## Put this at the top of the Snake file 


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

rule call:
  input:
      fa="data/genome.fa",
      bam=expand("mapped/{sample}.sorted.bam", sample=samples)
  output:
      "calls/all.vcf"
  conda:
      "envs/calling.yaml"
  shell:
      """    
      ## Note that {input.bam} includes mapped/A.sorted.bam, mapped/B.sorted.bam and mapped/C.sorted.bam
      samtools mpileup -g -f {input.fa} {input.bam} | \
      bcftools call -mv - > {output}
      """      
      

rule stats:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    conda:
        "envs/stats.yaml"
    script:
        "scripts/plot-quals.py"
        
        
      

      
      
