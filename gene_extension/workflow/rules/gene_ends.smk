

# those rules use flamseq alignment files to capture the end of the genes 


# get ends from the alignment:
rule get_ends:
	input: 'data/{sample}.sam'
	output: temp(".tmp/{sample}.ends.bed")
	shell: "python workflow/scripts/get_ends_fromsam_wrange2.py -i {input} -o {output}"	

# collapse ends:
rule collapse_ends:
	input: ".tmp/{sample}.ends.bed"
	output: temp(".tmp/{sample}.ends.bed.sorted"),"results/gene_ends/ends_by_library/{sample}.ends.counted.bed"
	params:
		distance_threshold = config['MERGE_DIST']
	shell: """
# input files are reads to ends mapping bed files that contain readIDs, read ends and the length of the longest 
# unmapped block at the end

# sort the file - bedtools doesn't allow stdin input 
bedtools sort -i {input} > {output[0]}
# -s -c 6 -o distinct
# -d -1 merge ends overlapping by 1 base at least

# change from reporting unique read ids to all - multimappers
#bedtools merge -d -1 -s -c 1,6,5,4 -o count,distinct,max,collapse -i {output[0]} > {input}.counted
bedtools merge -d -1 -s -c 1,6,5,4,7 -o count,distinct,max,distinct,max -i {output[0]} > {input}.counted



# the width of the ends is not different from 1 - depends on the merging distance
# now append the unique name

# the last column of a file stores the longest ungapped block per end:
awk -F "\t" -v a={wildcards.sample} '{{print $1"\t"$2"\t"$3"\t"a":FL"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' {input}.counted > {output[1]}

"""
# merge ends from all libraries:

rule merge_ends:
	input: expand('results/gene_ends/ends_by_library/{sample}.ends.counted.bed',sample = config['libraries'])
	output: 'results/gene_ends/all_ends.bed','results/gene_ends/end_support.tab'
	params: 
		merge_dist = -1
	shell: """ 
mkdir -p .tmp
cat {input} | bedtools sort > .tmp/all_ends.bed
bedtools merge -s -d {params.merge_dist} -c 4,5,6,7,9 -o collapse,sum,distinct,max,max -i .tmp/all_ends.bed > .tmp/all_ends_longids.bed
# get short flamend ids --> store the mapping
awk '{{print $1"\t"$2"\t"$3"\tfl"NR"\t"$5"\t"$6"\t"$7"\t"$8}}' .tmp/all_ends_longids.bed > {output[0]}
# this table contains final flamend IDs and the names obtained by merging identical ends from multiple libraries - it can be used to parse end read support afterwards
awk '{{print $4"\tfl"NR}}' .tmp/all_ends_longids.bed > {output[1]}

"""
