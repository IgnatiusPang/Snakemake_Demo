# Snakemake_Demo
Test the use of Snakemake on a PBSPro cluster (e.g. Katana).
These codes are based on the Snakemake short tutorial available on [here](https://snakemake.readthedocs.io/en/stable/tutorial/short.html).


# Setup 
Need to download the scripts and then the test data.
```bash
git clone git@github.com:IgnatiusPang/Snakemake_Demo.git
cd Snakemake_Demo
wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz
tar --wildcards -xf v5.4.5.tar.gz --strip 1 "*/data"
```

# Run 
* -j controls the number of jobs that can be submitted at one time 
* The `Snakefile` contains the entire demo pipeline
* The Step*.smk includes the steps outlined in the tutorial. 
* The following command will run the entire pipeline, including installation of software, on Katana. Jobs will be automatically submitted by snakemake
```bash
snakemake -j 3 --cluster-config cluster.json --cluster "qsub -N {cluster.N} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime},mem={cluster.mem} -M {cluster.email} -m ae -j oe"
```

# Generate Report
```bash
snakemake --report report.html
```

# Generate Directed Acyclic Graph of the Pipeline
```bash
bash make_dag.sh
```

# Generate the Rule Graph of the Pipeline
```bash
bash make_rulegraph.sh
```


