# Extending O. sinensis genes with FLAM-seq data  

This folder contains a Snakemake pipeline to extend genes with FLAM-seq data.  
It starts with FLAM-seq data (`data/`) and does the follwing:  

1. Extracting cleavage sites from FLAM-seq data. 
2. Assigning cleavage sites to the genes. 
3. Creating artificial last exons and associated mRNA molecules in the genome annotation.  


To run a pipeline, you have to have `snakemake` installed.   
Clone repository: 
```
git clone https://github.com/rajewsky-lab/octopus_microRNAs.git
cd octopus_microRNAs/gene_extension 
```
First, try creating conda enviroment:  
```
conda env create --name geneext --file=workflow/envs/env.yaml
conda activate geneext
```
If conda environment has been generated successfully, then try executing the pipeline:  
```
snakemkae -p --use-conda --threads 1 
```
