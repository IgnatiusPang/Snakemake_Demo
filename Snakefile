# Tutorial source
# https://snakemake.readthedocs.io/en/stable/tutorial/short.html

## Make sure I start python 
# module load python/3.7.7

# qsub -I -N SmkDemo -l nodes=1:ppn=12,walltime=12:00:00,mem=60gb


### Transfer files with these commands 
# rsync -artuvz /home/ignatius/PostDoc/2020/StaphCLASH2020/Source/Test_Snakemake/*  clive:/home/z3371724/snakemake-demo
# rsync -artuvz /home/ignatius/PostDoc/2020/StaphCLASH2020/Source/Test_Snakemake/*  kdm:/home/z3371724/snakemake-demo

# rsync -artuvz clive:/home/z3371724/snakemake-demo/dag3.pdf  /home/ignatius/PostDoc/2020/StaphCLASH2020/Source/Test_Snakemake/
# rsync -artuvz clive:/home/z3371724/snakemake-demo/rulegraph.pdf  /home/ignatius/PostDoc/2020/StaphCLASH2020/Source/Test_Snakemake/

## Dry run 
# snakemake --use-conda -n mapped/A.bam

## Running with specified number of cores 
# snakemake --use-conda --cores 4  mapped/A.bam

## -j controls the number of jobs that can be submitted at one time 
#  snakemake -j 2 --cluster-config cluster.json --cluster "qsub -N {cluster.N} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime},mem={cluster.mem} -M {cluster.email} -m ae -j oe"
#  snakemake -j 3 --cluster-config cluster.json --cluster "qsub -N {cluster.N} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime},mem={cluster.mem} -M {cluster.email} -m ae -j oe"
# snakemake --report report.html

samples = ["A", "B", "C"]

rule all:
    input:
      "calls/all.vcf",
      "plots/quals.svg"

# localrules: stats

rule bwa:
  input:
    "data/genome.fa",
    "data/samples/{sample}.fastq"
  output:
    temp("mapped/{sample}.bam"),  ## Make it a temporary file, will be deleted automatically when it is no longer needed
  conda: 
    "envs/mapping.yaml"
  threads: 8  ## Give it eight threads 
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
  shell:
    """
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
    "envs/call.yaml"  
  shell:
    """
    module load samtools/1.9
    samtools mpileup -g -f {input.fa} {input.bam} | bcftools call -mv - > {output}
    """


rule stats:
  input:
    "calls/all.vcf"
  output:
    report("plots/quals.svg", caption="report/calling.rst")  ## Create report 
  conda:
    "envs/stats.yaml"
  script:
    "scripts/plot-quals.py"
  
    


    
