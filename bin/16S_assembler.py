#!/usr/bin/env python

### - - - Script for extracting clustered 16S rRNA reads from datasets - - - ###

###- Importing modules -#
import sys
import code
try:
    from Bio import SeqIO
except ImportError:
    print "Please install Biopython first"
    exit(1)

###- Determine inputs -#

guppy_file = sys.argv[1]
aligned_reads_file = sys.argv[2]

###- Store the aligned sequence file into a set. Also store the length of the gene of interest -#

aligned_reads_set = set(SeqIO.parse(open(aligned_reads_file, 'r'), 'fasta'))

gene_length = len(list(SeqIO.parse(open(aligned_reads_file, 'r'), 'fasta'))[1].seq)
sequence_list = []
###- Creating dictionaries where each entry is a read, and it's item it the taxonomy assigned by pplacer -#

read_tax_dict = {}

###- Creating lists with the range of genera etc within each sample -#
taxonomy_div = {}


###- Creating dictionary with all reads, and taxonomy assigned to each -#

for line in open(guppy_file,'r'):
    lst = list(line.rstrip().split())
    if lst[0] != 'name' and lst[1] == lst[2] and float(lst[len(lst)-2]) > float(0.5):
        if lst[3] not in taxonomy_div:
            taxonomy_div[lst[3]] = []
        if lst[0] not in read_tax_dict:
            read_tax_dict[lst[0]] = []
            read_tax_dict[lst[0]].append(lst[3])
        else:
            read_tax_dict[lst[0]].append(lst[3])
print taxonomy_div
###- Stratifying this main dictionary into dictionaries by different taxonomic rank. At the same time, building -#


for taxonomy_index in reversed(range(1,8)):
    for type in taxonomy_div:
        print '>'+type
        for read_name, taxonomy in read_tax_dict.iteritems():
            if len(taxonomy) == taxonomy_index:
                if type == taxonomy[taxonomy_index-1]:
                    for x in aligned_reads_set:
                        if x.id == read_name:
                            print x.seq    
            

exit(1)