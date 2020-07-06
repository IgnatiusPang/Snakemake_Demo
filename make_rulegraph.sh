

snakemake \
    --rulegraph \
    | dot -Tpdf \
    > rulegraph.pdf
    
    