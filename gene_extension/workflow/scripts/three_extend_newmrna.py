# Extend genome annotation - create novel exons
verbose = True

# Given a genome annotation with extended gene ranges, extend the last exon of one transcript belonging to the gene

# TODO: simplify the script - it will take only the table with the extensions gene_id \t ext 

# prefix appended to the fake transcript IDs.
tag = 'FLAMseq'

import argparse
import os
import pandas as pd
import gffutils
import os.path
from os import path

# parser part belongs here 
parser = argparse.ArgumentParser(description='Add fake 3-prime exon transcripts to the genes based on gene ends file.')
parser.add_argument('-g', default= None,help = 'Genome .gtf/.gff file' )
parser.add_argument('-e', default= None,help = 'Gene ends table' )
parser.add_argument('-o', default = None, help = 'Output annotation')

args = parser.parse_args()


genome_file = args.g
gends_file = args.e
output_file = args.o

def transform(d):
    if d['gene_id']:
        if 'Trn' in d['gene_id']:
            print('Encountered tRNA gene, passing...')
            pass
        else:
            #print(d['gene_id'])
            d['gene_id'] = str(d['gene_id'][0]).replace(';','').replace('"','').replace('";','')
            d['gene_name'] = d['gene_id']
        try: 
            if d['transcript_id']:
                d['transcript_id']  = str(d['transcript_id'][0]).replace(';','').replace('"','').replace('";','')
        except KeyError:
            pass       
    return(d)
print('Loading %s into a database ... ' % genome_file)

# Add loading from the database for debugging purposes:
dbname = 'results/gene_extension/genes.db'
if not path.exists(dbname):
    print('db not found, creating ...')
    db = gffutils.create_db(genome_file, dbname,transform = transform ,disable_infer_genes=True,disable_infer_transcripts=True, merge_strategy = 'create_unique',id_spec = {'gene': 'gene_id', 'transcript': 'transcript_id'})
    print('done.')
else:
    db = gffutils.FeatureDB(dbname, keep_order=True)


print('DEBUG')
print(dir(db))
print('Feature types loaded:\n' + "\n".join([x for x in db.featuretypes()]))
print('Number of genes: %s' % str(db.count_features_of_type("gene")))
#for feature in db.features_of_type('transcript'):
#    print(feature.id)
print('done')
# later replace by this 
#db = gffutils.create_db(genome_file, ":memory:",transform = transform ,disable_infer_genes=True,disable_infer_transcripts=True, merge_strategy = 'create_unique')



#def transform(d):
#    print(d)
#    if d['gene_id']:
#        d['gene_id'] = str(d['gene_id'][0]).replace(';','').replace('"','').replace('";','')
#        #print(d['gene_id'])[0]
#    try: 
#        if d['transcript_id']:
#    #        print('is trid')
#    #        print(d['transcript_id'])
#            d['transcript_id']  = str(d['transcript_id'][0]).replace(';','').replace('"','').replace('";','')
#    except KeyError:
#        pass         
#    return(d)
#print('Loading the annotation')
#dbname = 'results/gene_extension/genes.db' 
#if not path.exists(dbname):
#    db = gffutils.create_db(genome_file, dbfn = dbname,transform = transform,disable_infer_genes=True,disable_infer_transcripts=True,merge_strategy = "create_unique" )
#    print('done.')
#else:
#    db = gffutils.FeatureDB(dbname, keep_order=True)
#db = gffutils.create_db(genome_file, dbfn = ":memory:",transform = transform,disable_infer_genes=True,disable_infer_transcripts=True,merge_strategy = "create_unique" )
print('done.')




# This dictionary contains extensions - just the number of downstream bases to extend the gene to
# TODO: the script should also be able to accept gene_end.bed will all the ends per gene -- filtering, extension
gene_ends = pd.read_csv(gends_file,sep = '\t')
print("Gene ends loaded")
# select the maximum extension.  
gene_ends = gene_ends.groupby(['gene'], sort=False)['extension'].max()
#gene_ends = dict(zip(gene_ends.gene,gene_ends.))
gene_ends = gene_ends[gene_ends!=0]
gene_ends = gene_ends.to_dict()
print('Extending %s genes.' % len(gene_ends))

print("Writing output...")
# you have to also define an mRNA corresponding
# for some reason, EVM genes ids are wrongly parsed


# New loop - iterate over genes
# so exons are not assigned 
#cnt=0
#for gene in db.features_of_type('gene'):
#    print(cnt)
#    if cnt > 3:
#        break
#    if gene.chrom !="NC_006353.1" :
#        print(gene.id)
#        # chilren are not inferred properly.
#        for feature in db.children(gene.id):
#            print("%s:%s"  % (feature.featuretype,feature.id))
#        cnt+=1
#    else:
#        pass

# children seeem to be inferred properly
print('Looping through the genome:')
cnt = 0
verbose = True

import numpy as np
print(np.unique([db[x].strand for x in gene_ends.keys()],return_counts = True))
genes_left = list(gene_ends.keys())
genes_noex = [] # container for genes without detected exons
with open(output_file, 'w') as fout:
    for d in db.directives:
        fout.write('##{0}\n'.format(d))
    for i,feature in enumerate(db.features_of_type("gene")):
        #if i> 100:
        #    exit()
        if feature['gene_id']:
            feature['gene_id'] = feature['gene_id'][0]
        else:
            print(feature)
            exit("no gene_id attribute found for gene %s! ( see the feature line above)" % feature.id)
        # remove annoying "note" attribute:
        if 'note' in [x for x in feature.attributes]: 
            last_exon['note'].pop()
        # Now, this is a gene, check whether there are any exons:
        n_exons = len([x for x in db.children(db[feature.id],featuretype='exon')])
        gene = feature.id
        if n_exons and feature.id in gene_ends.keys(): # if no exons found for a gene - it's not written CAVE MT genes!
            print("%s exons found" % n_exons)
            # Select the most downstream exon:
            # 3'most exon
            if verbose:
                print("    %s strand." % db[gene].strand)
            if db[gene].strand == '+':
                # update gene bound:
                # gene end is probably not updated properly
                feature.end = db[gene].end + gene_ends[feature.id]
                max_end = max([f.end for f in db.children(db[gene],featuretype='exon')])
                last_exon = [x for x in db.children(db[gene],featuretype='exon') if x.end == max_end ][0]
                last_exon.end = last_exon.end + gene_ends[feature.id]
                if verbose:
                    print("Gene end chage: %s --> %s; %s" % (str(last_exon.end),str(last_exon.end+gene_ends[feature.id]),str(gene_ends[feature.id])))
            elif db[gene].strand == '-':
                feature.start = db[gene].start - gene_ends[feature.id]
                max_end = min([f.start for f in db.children(db[gene],featuretype='exon')])
                last_exon = [x for x in db.children(db[gene],featuretype='exon') if x.start == max_end ][0]
                last_exon.start = last_exon.start - gene_ends[feature.id]
                if verbose:
                    print("Gene end chage: %s --> %s; %s" % (str(last_exon.start),str(str(last_exon.start-gene_ends[feature.id])),str(gene_ends[feature.id])))
            # Select, the most downstream exon and update:

            mrna_id = last_exon['transcript_id'][0]
            print("mRNA with the downstreamest exon: %s" % mrna_id)
            # change last_exon id 
            last_exon.id = '%s_lastexon~' % tag +last_exon.id
            # now, change the parent
            last_exon['transcript_id'][0] = 'FLAMseq_ext~'+last_exon['transcript_id'][0]
            # duplicate last exon, change the metadata:
            last_exon.source = last_exon.source + '~' + tag
            if 'model_evidence' in [x for x in last_exon.attributes]:
                last_exon['model_evidence'].pop()
            if 'model_evidence' in [x for x in last_exon.attributes]:
                last_exon['product'].pop()
            #last_exon['transcript_id'] = last_exon['transcript_id'][0]+'_FLAMseq_ext'
            last_exon['three_prime_ext'] = str(gene_ends[gene])
            ###### Create fictional mRNA ########
            # CAVE: mrna ranges should be also updated
            mrna = db[mrna_id]
            mrna['three_prime_ext'] = str(gene_ends[gene])
            mrna_old = mrna

            if mrna.strand == "+":
                mrna.end = last_exon.end
            elif mrna.strand == "-":
                mrna.start = last_exon.start
            mrna.source = mrna.source+'~'+tag
            # create fictional mrnaid
            mrna.id = 'FLAMseq_ext~' + mrna.id
            mrna['transcript_id'][0] = mrna.id
            ####### Write down the gene and extended mrna ######
            fout.write(str(feature)+'\n')
            fout.write(str(mrna)+'\n')

            print('adding %s: [%s/%s]' % (mrna.id,cnt,len(gene_ends)))
            cnt += 1

            for child_exon in db.children(mrna_id,featuretype = 'exon'):
                if not child_exon.id == last_exon.id:
                    #print(child_exon.id)
                    child_exon.id = '%s_exon~' % tag +child_exon.id
                    child_exon['transcript_id'][0] = 'FLAMseq_ext~'+mrna_id
                    child_exon.source = child_exon.source + '~'+tag
                    #print(child_exon.id)
                    #print(child_exon['transcript_id'][0])
                    fout.write(str(child_exon)+'\n')
            fout.write(str(last_exon) + '\n')
            # now, write remainining transcripts:
            for transcript in db.children(gene,featuretype = 'transcript'): 
                fout.write(str(transcript) + '\n')
                for exon in db.children(transcript.id):
                    if exon.featuretype in ['exon','CDS','start_codon','stop_codon']:
                    # write down the exon
                        fout.write(str(exon) + '\n')
            # remove the gene from the list:
            genes_left.remove(gene)
        elif not n_exons:
            print("no exons found for %s!" % gene)
            genes_noex = genes_noex + [gene]
        ######### Write to the file as is ############
        elif not gene in gene_ends.keys(): # write the gene and all children as they are in the file:
            #print("gene shoudn't be extended - omitting...")
            fout.write(str(feature) + '\n') # genes and transcripts are children of themselves
            for transcript in db.children(gene):
                if transcript.featuretype == 'transcript':
                    fout.write(str(transcript) + '\n')
                    for exon in db.children(transcript.id):
                        if exon.featuretype in ['exon','CDS','start_codon','stop_codon']:
                        # write down the exon
                            fout.write(str(exon) + '\n')  
        else:
            continue  
if len(genes_left)+len(genes_noex)>0:
    print("Genes have not been written: ")
    print("Genes that should have been extended but not seen?: %s" % len(genes_left))
    print("Genes without exons: %s" % len(genes_noex))
