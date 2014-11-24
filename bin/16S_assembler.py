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

###- Creating dictionaries where each entry is a read, and it's item it the taxonomy assigned by pplacer -#

read_tax_dict = {}
read_tax_dict_G = {}
read_tax_dict_F = {}
read_tax_dict_O = {}
read_tax_dict_C = {}
read_tax_dict_P = {}
read_tax_dict_K = {}

###- Creating lists with the range of genera etc within each sample -#
genera = []
families = []
orders = []
classes = []
phyla = []
kingdoms = []


###- Creating dictionary with all reads, and taxonomy assigned to each -#
for line in open(guppy_file,'r'):
    lst = list(line.rstrip().split())
    test = 1
    if lst[0] != 'name' and lst[1] == lst[2] and float(lst[len(lst)-2]) > float(0.5):
        if lst[0] not in read_tax_dict:
            read_tax_dict[lst[0]] = []
            read_tax_dict[lst[0]].append(lst[3])
        else:
            read_tax_dict[lst[0]].append(lst[3])

###- Stratifying this main dictionary into dictionaries by different taxonomic rank. At the same time, building -#

for read_name, taxonomy in read_tax_dict.iteritems():
    if len(taxonomy) == 7:
        read_tax_dict_G[read_name] = taxonomy[6]
        if taxonomy[6] not in genera:
                genera.append(taxonomy[6])

    if len(taxonomy) == 6:
        read_tax_dict_F[read_name] = taxonomy[5]
        if taxonomy[5] not in families:
                families.append(taxonomy[5])

    if len(taxonomy) == 5:
        read_tax_dict_O[read_name] = taxonomy[4]
        if taxonomy[4] not in orders:
                orders.append(taxonomy[4])

    if len(taxonomy) == 4:
        read_tax_dict_C[read_name] = taxonomy[3]
        if taxonomy[3] not in classes:
                classes.append(taxonomy[3])

    if len(taxonomy) == 3:
        read_tax_dict_P[read_name] = taxonomy[2]
        if taxonomy[2] not in phyla:
                phyla.append(taxonomy[2])

    if len(taxonomy) == 2:
        read_tax_dict_K[read_name] = taxonomy[1]
        if taxonomy[1] not in kingdoms:
                kingdoms.append(taxonomy[1])

###- For each member of the genus list, match up the reads that have this classification, and assemble the reads. -#
code.interact(local=locals())

for genus in genera:
    for read_id, classification in read_tax_dict_G.iteritems():
        if classification == genus:
            for x in aligned_reads_set:
                if x.id == read_id:
                    print x.seq
        else:
            continue
    exit(1)
exit(1)
