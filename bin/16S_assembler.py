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


import argparse
import subprocess
import code
import tempfile
from datetime import datetime
try:
    from Bio import SeqIO
except ImportError:
    print "Please install Biopython first"
    exit(1)
    
intro = '''
       #######################################
       ##         16S assembler %s       ##
       ##               beta                ##
       ##             Joel Boyd             ##
       #######################################  
       
''' % __version__

class Messenger:
    
    def __init__(self, message):
        self.message = message
        self.time = datetime.now().strftime('%H:%M:%S')
    
    
    
    def date_n_message(self):
        return '  [%s]: %s' % (self.time, self.message)
    
    def error_message(self):
        return '  [%s] == ERROR == : %s' % (self.time, self.message)



class TaxoGroup:       
        
    def guppy_splitter(self, guppy_file):  
        
        print Messenger("Parsing guppy file").date_n_message()
        read_tax_dict = {}
        
        for line in open(guppy_file,'r'):
        
            lst = list(line.rstrip().split())
            
            if lst[0] != 'name' and lst[1] == lst[2] and float(lst[len(lst)-2]) > float(0.5):
                
                if lst[0] not in read_tax_dict:
                    read_tax_dict[lst[0]] = []
                    read_tax_dict[lst[0]].append(lst[3])
                
                else:
                    read_tax_dict[lst[0]].append(lst[3])
        
        return read_tax_dict
          
     
    
    def rank_grouper(self, dictionary, sequences_file):
        
        print Messenger("Sorting sequences").date_n_message()
        sequences = set(SeqIO.parse(open(sequences_file, 'r'), 'fasta'))
        tax_list = [{}, {}, {}, {}, {}, {}]

        for read_name, classification in dictionary.iteritems():
            
            tax_group = classification[len(classification) - 1]
            tax_dict = tax_list[len(classification) - 2]
            
            if tax_group not in tax_dict:
                    tax_dict[tax_group] = []
                    
            for sequence in sequences:
                
                if sequence.id == read_name:
                    tax_dict[tax_group].append(sequence)
                    
        return tax_list
        
    
    

    def assembler(self, tax_rank_list, output_directory, assembly_type):
        
        rank = ['K','P','C','O','F','G']
        cmd = 'mkdir %s ' % output_directory
        
        try:
            subprocess.check_call(cmd, shell = True)
        
        except:
            print Messenger("A folder with the name '%s' already exists... This is awkward... \n" % output_directory).error_message()    
            exit(1)
        
        print Messenger("Assembling").date_n_message()
        
        for idx, taxonomic_rank in enumerate(tax_rank_list):
            
            for tax_name, sequences in taxonomic_rank.iteritems():
                number_of_sequences = len(list(sequences))
                
                if number_of_sequences > 10:
                    
                    output_sequence_file_path = '%s_%s_%s_assembly.fa' % (rank[idx], tax_name, str(number_of_sequences))
                    
                    if assembly_type == "cap3":
                        
                        with tempfile.NamedTemporaryFile() as tmp:
                            
                            SeqIO.write(sequences, tmp.name, "fasta")
                            cmd = 'cap3 %s > %s.junk' % (tmp.name, tmp.name)
                            subprocess.check_call(cmd, shell = True)

                            cmd = 'mv %s.cap.contigs ./%s/%s' % (tmp.name, output_directory, output_sequence_file_path)
                            subprocess.check_call(cmd, shell = True)   
                             
                    elif assembly_type == "phrap":
                        
                        with tempfile.NamedTemporaryFile() as tmp:
                                    
                            SeqIO.write(sequences, tmp.name, "fasta")
                            cmd = 'phrap  %s >/dev/null 2>&1' % (tmp.name)
                            subprocess.check_call(cmd, shell = True)
                                  
                            cmd = 'mv %s.contigs ./%s/%s' % (tmp.name, output_directory, output_sequence_file_path)
                            subprocess.check_call(cmd, shell = True)
                    
                    elif assembly_type == "hmm":
                        
                        gene_size = len(list(sequences[0].seq))
                        nucleotide_composition = [{"A": [0]*gene_size}, {"T": [0]*gene_size}, {"G": [0]*gene_size}, {"C": [0]*gene_size}]
                        total = [0]*gene_size
                        final_sequence = []

                        for sequence in sequences:
                            sequence_list = list(sequence)
                            
                            for idx, nucleotide in enumerate(sequence_list):
                                if nucleotide == "A":
                                    nucleotide_composition[0]["A"][idx] += 1
                                
                                elif nucleotide == "T":
                                    nucleotide_composition[1]["T"][idx] += 1
                                    
                                elif nucleotide == "G":
                                    nucleotide_composition[2]["G"][idx] += 1
                                    
                                elif nucleotide == "C":
                                    nucleotide_composition[3]["C"][idx] += 1
                                    
                                else:
                                    Messenger("Programming error... \n").error_message()
                        
                        
                        for i in range(0, gene_size):
                            idx = i - 1
                            
                            for position in total:
                                
                                for nucleotide_dictionary in nucleotide_composition:
                                    
                                    for key, item in nucleotide_dictionary.iteritems():
                                        
                                       total[idx] += item[idx]
                                       
                        for idx, position in enumerate(total):
                                                              
                            for key, item in nucleotide_dictionary.iteritems():
                                
                                if item[idx] > 0:
                                    print item[idx]
                                    print total[idx]
                                    print item[idx]/total[idx]
                                    item[idx] = item[idx]/total[idx] 
                                    
                                else:
                                    continue
                        
                        print nucleotide_composition
                            
                                    
                    
                    else:
                        print Messenger("Programming error").error_message()
                        exit(1)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='''--- 16S assembler %s --- Groups and assembles 16S reads from a guppy file by taxonomic rank and classification.''' % __version__
                                , epilog='==============================================================================')
    parser.add_argument('--guppy', metavar='guppy file produced by graftM',  nargs=1 , help='guppy file produced by graftM', required=True)
    parser.add_argument('--reads_file', metavar='reads file',  nargs=1, help='unaligned reads file produced by graftM', default=argparse.SUPPRESS)
    parser.add_argument('--output', metavar='output directory', nargs=1, help='Directory within which to place assembled sequences', required=True)
    parser.add_argument('--assembly_type', metavar='type of assembly', nargs=1, help='cap3 or native hmmalignment assembly',choices = ['cap3','hmm', 'phrap'], required=True)
    args = parser.parse_args()
    
    print intro
    
    split_guppy = TaxoGroup().guppy_splitter(args.guppy[0])
    grouped_sequences = TaxoGroup().rank_grouper(split_guppy, args.reads_file[0])
    TaxoGroup().assembler(grouped_sequences, args.output[0], args.assembly_type[0])
            
    exit(1)    





    

    