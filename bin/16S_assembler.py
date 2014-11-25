#!/usr/bin/env python


# ============================================================== #
#  Script for extracting clustered 16S rRNA reads from datasets  #
#                              v.0.0.1                           #
# ============================================================== #

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2014"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
__version__ = "0.0.1"

# --- Importing modules
import sys
import code
import subprocess
from datetime import datetime
try:
    from Bio import SeqIO
except ImportError:
    print "Please install Biopython first"
    exit(1)



# --- Determine inputs

guppy_file = sys.argv[1]
reads_file = sys.argv[2]




# --- Store the aligned sequence file into a set. Also store the length of the gene of interest 

aligned_reads_set = set(SeqIO.parse(open(reads_file, 'r'), 'fasta'))
gene_length = len(list(SeqIO.parse(open(reads_file, 'r'), 'fasta'))[1].seq)



# --- Defining lists

hits_dictionary = {} #- Dictionary of taxonomic ranks with > 1 sequence assigned to it
sequence_list = []
read_tax_dict = {}

tax_list = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]



# --- Function that sends message with date and time

def message(message):
    time = datetime.now().strftime('%H:%M:%S')
    print '[' + time + ']: ' + str(message)



# --- Creating dictionary with all reads, and taxonomy assigned to each

for line in open(guppy_file,'r'):
    
    lst = list(line.rstrip().split())
    
    if lst[0] != 'name' and lst[1] == lst[2] and float(lst[len(lst)-2]) > float(0.5):
        
        if lst[0] not in read_tax_dict:
            read_tax_dict[lst[0]] = []
            read_tax_dict[lst[0]].append(lst[3])
        
        else:
            read_tax_dict[lst[0]].append(lst[3])



# --- Finding list of diversity

for read_name, taxonomy in read_tax_dict.iteritems():
    taxonomy.append(tax_list[len(taxonomy)-2])



# --- Stratifying this main dictionary into dictionaries by different taxonomic rank. Write them to a file.

for taxonomy_index in reversed(range(1,8)):
    
    if taxonomy_index > 1:
        message('Searching at the %s taxonomic rank' % tax_list[taxonomy_index-2])   
        
        for read_name, taxonomy in read_tax_dict.iteritems():
            
            if taxonomy[len(taxonomy)-1] == tax_list[taxonomy_index-2]:                  

                if taxonomy[len(taxonomy)-2] not in hits_dictionary:
                    hits_dictionary[taxonomy[len(taxonomy) - 2]] = []                                  

                for sequence in aligned_reads_set:
                    
                    if sequence.id == read_name:
                        hits_dictionary[taxonomy[len(taxonomy) - 2]].append(sequence)             
                         
                         
                                                
# --- Assembly using CAP3 (for now)

message('Writing parsed sequences to files')

for tax_name, sequences in hits_dictionary.iteritems():
    output_sequence_file_path = open('%s_sequences.fa' % tax_name, 'w')
    
    for sequence in sequences:
        output_sequence_file_path.write('>'+sequence.id+'\n')
        output_sequence_file_path.write(str(sequence.seq+'\n'))
        
    
    