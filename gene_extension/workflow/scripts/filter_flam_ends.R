options(max.print = 100)
message("Usage: script.R < minimum of FLAMseq reads supporting the end > < genes.gtf file > < output file > < absolute length threshold, bp > <closest genes file>")
# for interactive run
#setwd('~/octopus/ovu/flamseq_clean')
#input: config['GENOMEANNO'] , ".tmp/closest_genes"
#output: temp(".tmp/ext_genes.gtf"),'.tmp/extension_bp.tsv','.tmp/genes_ends.tsv'
#Rscript workflow/scripts/add_flam_ends.R {params.count_t} {input[0]} {output[0]} {params.abs_len_t} {input[1]} > extension.log
#options = c('1','.tmp/genes.gtf','.tmp/ext_genes.gtf','60000','.tmp/closest_genes')
#setwd('~/octopus/ovu/3ends_flamseq/.tmp/')
#setwd('~/octopus/ovu/flamseq_clean')
#options = c('1','.tmp/genes.gtf','60000','.tmp/closest_genes')


# the next script actually reads .tmp/genes_ends.tsv file ...
# whaaaaat the fuck ????

# Load .tmp/closest_genes file that specifies the closest gene per each end.
# Filter the file and write the extensions into a new data frame.
suppressPackageStartupMessages(library(dplyr))
options <- commandArgs(trailingOnly = TRUE)
if(length(options)<4){
	stop("Not enough arguments!")
}
options(stringsAsFactors = F)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(stringr))
print(options)
coverage_t = as.integer(options[1])
genes_file = options[2]
# sets the length to an absolute deviation from a gene end to exlcude novel genes (if there is a missed gene, the distance will be very high)
abs_length_t = as.integer(options[3])
closest_g_file = options[4]

cl = read.table(closest_g_file,sep ='\t')
message(paste0(nrow(cl),' ends loaded.'))
colnames(cl)[ncol(cl)] = 'distance'
print(dim(cl))
message('Distance summary:')
summary(cl$distance)
# remove ends without hits (can represent novel genes!):
message(sprintf('Ends without upstrean gene: %s',sum(cl[,ncol(cl)]==-1)))
#cl$gene_id = sub('(gene_id )|ID ','',sapply(str_split(cl$V15,';'),FUN = function(x) x[1]))
cl = cl[cl$V17!='.',]
cl$gene_id = gsub('gene_id ','',str_extract(cl$V17,'gene_id [A-Z,a-z,_,0-9]+'))
message(sprintf('%s ends are mapping to %s genes.',nrow(cl),length(unique(cl$gene_id ))))
# filter by the coverage:
clf = cl[cl$V5>=coverage_t,]
print(clf$gene_id)
message(sprintf("%s ends passed coverage threshold of %s reads",nrow(clf),coverage_t))
# V16 stores the distance from flamseq end to the gene, append this distance as a separate feature 
clf = clf[abs(as.integer(clf$distance))<=abs_length_t,]
# how will you tell upstream from downstream ends? 

# CAVE: the distance of upstream features will be negative - reverse
clf$distance = abs(clf$distance)


print('after filtering by distance')
print(dim(clf))
print('distance == 0')
print(table(clf$distance==0))
# only 10k genes - what the fuck, grisha
message(sprintf("%s ends passed distance threshold <= %skb",nrow(clf),abs_length_t/1e3))
# now, select maximum one end per gene by distance 
# select the most upstream gene per end:
message('Selecting the most downstream ends per gene...')
clf = clf%>%dplyr::group_by(gene_id)%>%dplyr::slice(which.max(abs(distance)))%>%as.data.frame()
clf = clf[clf$distance>0,]
# what are the negative distances???
message(sprintf('Most-dowstream ends for %s genes retained.',length(unique(clf$gene_id))))
clf=clf[!is.na(clf$gene_id),]
rownames(clf) = clf$gene_id
# remove mitochondrial genes
clf = clf[clf$V1!="NC_006353.1",]
# now, keep in mind that GTF file is 1-based as opposed to bed and bam files.???
clf$V2 = clf$V2+1
clf$V3 = clf$V3+1
# remove flamseq ends that overlap genomic annotation - ends are contained withing genomic range
clf = clf[clf$distance>0,]
message(paste0(nrow(clf),' genes with downstream ends.'))
# ---- Write output table -------
# originally, the output was genes_ends.tsv - all genes vs all ends - makes sense for the future analysis 
# Output table
#gene end extension 
#COX3 fl151 355
#ND3 fl173 70
ext_table = data.frame(gene = clf$gene_id,end = clf$V4,extension = abs(clf$distance))
write.table(ext_table,'.tmp/genes_ends.tsv',row.names =F, quote = F,sep = '\t')
