# For the alignment file, get read ends in a bed file 
import pysam
import argparse

parser = argparse.ArgumentParser(description='For the .sam file, get read ends in a .bed file')
parser.add_argument('-i', default= None,help = 'Input sam' )
parser.add_argument('-o', default = None, help = 'Output bed')

args = parser.parse_args()


infile = args.i
outfile = args.o

#infile = '/home/gzolota/octopus/ovu/3ends_flamseq/data/flamseq_wpolya.sam'
#outfile = '/home/gzolota/octopus/utrs/flamends/toy_ends.bed'

samfile = pysam.AlignmentFile(infile, "rb")
out = open(outfile, "w")
# For each read, get the alignment end
# CAVE: reference end points to the one position past the last aligned base 
# maybe, polya lengths are in play
for i,read in enumerate(samfile):
    if i%100000 == 0:
        print('%s reads processed ...' % i)
    end = {'chrom':None,'start':None,'end':None,'strand':None,'name':None}
    if not read.is_unmapped:
        #print(read.query_name)
        if read.is_reverse:
            #print('plus strand')
            # get the reference end 
            end['chrom'] = read.reference_name
            end['start'] = read.reference_end-1
            end['end'] = read.reference_end
            end['strand'] = '+'
            end['name'] = read.query_name
            #print(end)
        elif not read.is_reverse:
            #print('minus strand')
            end['chrom'] = read.reference_name
            end['start'] = read.reference_start
            end['end'] = read.reference_start+1
            end['strand'] = '-'
            end['name'] = read.query_name
            #print(end)
        # write ends to the file 
        out.write(str(end['chrom'])+"\t"+str(end['start'])+"\t"+str(end['end'])+"\t"+str(end['name'])+"\t"+str(40)+"\t"+str(end['strand'])+'\n')
out.close()        
print('done')

