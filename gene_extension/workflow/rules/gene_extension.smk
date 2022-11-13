# These rules take as input gene ends file and produce genome annotation with extended gene regions

# CAVE: the extensions can not overlap other genes ( in principle :)

###### Extend the genes ######
        #input: "results/all_ends.bed", "resources/genes.gtf"

# 10.01.21 - changed -id -s -D ref to  -id -s -D a - all downstream features are ignored, all the distances of upstream elements should be 0 or negative.
rule get_closest_upstream_genes_for_ends:
      	input: "results/gene_ends/all_ends.bed",config['GENOMEANNO']
	output: temp(".tmp/closest_genes_plus"),temp(".tmp/closest_genes_minus"), ".tmp/closest_genes"
    conda:  "../envs/env.yml"
        shell: """
awk '$6=="+"' {input[0]} > .tmp/ends_plus.bed
awk '$6=="-"' {input[0]} > .tmp/ends_minus.bed

# get gene regions from the file:

awk '$3=="gene"' {input[1]} > .tmp/gregs.gtf
bedtools sort -i .tmp/gregs.gtf > tmp; mv tmp .tmp/gregs.gtf

bedtools closest -id -s -D a -a .tmp/ends_plus.bed -b .tmp/gregs.gtf > {output[0]}
bedtools closest -id -s -D a -a .tmp/ends_minus.bed -b .tmp/gregs.gtf > {output[1]}

cat {output[0]} {output[1]} > {output[2]}
bedtools sort -i {output[2]} > tmp; mv tmp {output[2]}
rm .tmp/gregs.gtf
"""

rule filter_flam_ends:
	input: config['GENOMEANNO'] , ".tmp/closest_genes"
	output: '.tmp/genes_ends.tsv' 
	params:
		count_t = config["COUNT_T"],
		abs_len_t = config["ABS_LEN_T"]
    conda:  "../envs/env.yml"
	shell: "mkdir -p results/gene_extension/; Rscript workflow/scripts/filter_flam_ends.R {params.count_t} {input[0]} {params.abs_len_t} {input[1]} > results/gene_extension/extension.log"
# add fake exons to each gene so their ends will correspond to the furtherst flamseq end for the gene
rule update_genome_artficial_exons:
	input: 
		annotation = config['GENOMEANNO'], 
		ends = ".tmp/genes_ends.tsv"
	output: "results/gene_extension/genome_extended_wexons.gtf","results/gene_extension/genes_ends.tsv"
    conda:  "../envs/env.yml"
	shell: """ 

mkdir -p results/gene_extension;
python workflow/scripts/three_extend_newmrna.py -g {input.annotation} -e  {input.ends} -o {output[0]};
sed -i  's/;  "";//g' {output[0]}
mv {input.ends} {output[1]}
"""

