#!/usr/bin/env python
import cmd
from cmd import Cmd


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
        return '  [%s]: == ERROR == %s' % (self.time, self.message)



class TaxoGroup:       
        
    def guppy_splitter(self, guppy_file):  
        
        print Messenger("Parsing guppy file").date_n_message()
        read_tax_dict = {}
        
        for line in open(guppy_file,'r'):
            
            lst = list(line.rstrip().replace('Candidatus ', 'Candidatus_').split())
            
            if lst[0] != 'name' and lst[1] == lst[2] and float(lst[len(lst)-2]) > float(0.5):
                
                if lst[0] not in read_tax_dict:
                    read_tax_dict[lst[0]] = []
                    read_tax_dict[lst[0]].append(lst[3])
                
                else:
                    read_tax_dict[lst[0]].append(lst[3])
        
        return read_tax_dict
          

    def rank_grouper(self, dictionary, sequences_file):
        
        print Messenger("Sorting 16S sequences").date_n_message()
        sequences = set(SeqIO.parse(open(sequences_file, 'r'), 'fasta'))
        tax_list = [{}, {}, {}, {}, {}, {}]
        reverse_reads = []
        
        for read_name, classification in dictionary.iteritems():
            
            tax_group = classification[len(classification) - 1]
            tax_dict = tax_list[len(classification) - 2]
            
            if tax_group not in tax_dict:
                    tax_dict[tax_group] = []
            
            for sequence in sequences:
                
                if sequence.id == read_name.replace('_',':'):
                    reverse_reads.append(sequence.id)
                    tax_dict[tax_group].append(sequence)
        
        return tax_list, reverse_reads
        
    
    def all_read_extractor(self, sequences, reads):
        
        print Messenger("Extracting reverse reads to join assembly").date_n_message()
        
        with open('tmp_allreads.txt', 'w') as ids:
            
            for read in reads:
                ids.write(read+'\n')
        
        cmd = "fxtract -X -H -f tmp_allreads.txt %s | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > tmp_allreads.fa" % (sequences)
                
        subprocess.check_call(cmd, shell = True)
        
    def read_extractor(self, sequences):
        
        with open('tmp_grouped.txt', 'w') as ids:
            
            for sequence in sequences:
                ids.write(sequence.id+'\n')
        
        cmd = "fxtract -X -H -f tmp_grouped.txt tmp_allreads.fa > tmp_grouped.fa" 
        subprocess.check_call(cmd, shell = True)

        return
    
    def assembler(self, tax_rank_list, output_directory, assembly_type):
        
        rank = ['K','P','C','O','F','G']
        cmd = 'mkdir %s 2>/dev/null' % output_directory
        
        try:
            subprocess.check_call(cmd, shell = True)
        
        except:
            print Messenger("A folder with the name '%s' already exists... This is awkward... \n" % output_directory).error_message()    
            exit(1)
        
        print Messenger("Assembling using %s" % assembly_type).date_n_message()
        
        for idx, taxonomic_rank in enumerate(tax_rank_list):
            
            for tax_name, sequences in taxonomic_rank.iteritems():
                number_of_sequences = len(list(sequences))
                
                if number_of_sequences > 10:
                    output_sequence_file_path = '%s_%s_%s_assembly.fa' % (rank[idx], tax_name, str(number_of_sequences))
                    
                    if assembly_type == "cap3":
                                                
                        with tempfile.NamedTemporaryFile() as tmp:
                            SeqIO.write(sequences, tmp.name, "fasta")
                            cmd = 'cap3 %s -o 16 -p 98 > %s.junk' % (tmp.name, tmp.name)
                            subprocess.check_call(cmd, shell = True)

                            cmd = 'mv %s.cap.contigs ./%s/%s' % (tmp.name, output_directory, output_sequence_file_path)
                            subprocess.check_call(cmd, shell = True)  
                             
                    elif assembly_type == "omega":
                        
                        with tempfile.NamedTemporaryFile() as tmp:
                            SeqIO.write(sequences, tmp.name, "fasta")
                            cmd = 'omega -l 25 %s' % (tmp.name)
                            subprocess.check_call(cmd, shell = True)

                            cmd = 'mv %s.cap.contigs ./%s/%s' % (tmp.name, output_directory, output_sequence_file_path)
                            subprocess.check_call(cmd, shell = True) 
                        
                    elif assembly_type == "phrap":

                        with tempfile.NamedTemporaryFile() as tmp:
                            
                            SeqIO.write(sequences, tmp.name, "fasta")
                            cmd = 'phrap -minscore 40 -minmatch 45 -revise_greedy -forcelevel 0 -repeat_stringency 0.99 %s >/dev/null' % (tmp.name)

                            subprocess.check_call(cmd, shell = True)
                            
                            cmd = 'mv %s.contigs ./%s/%s' % (tmp.name, output_directory, output_sequence_file_path)
                            subprocess.check_call(cmd, shell = True)
                    
                    elif assembly_type == "felvet":
                                                   
                        with tempfile.NamedTemporaryFile() as tmp:
                            
                            if args.reverse_reads:
                                
                                TaxoGroup().read_extractor(sequences)
                                SeqIO.write(sequences, tmp.name, "fasta")
                                                                
                                cmd = 'cat tmp_grouped.fa %s > %s' % (tmp.name, output_sequence_file_path.replace('_assembly.fa', ''))
                                subprocess.check_call(cmd, shell = True)
                                                       
                                cmd = 'finishm visualise --quiet --assembly-coverage-cutoff 1.5 --assembly-kmer %s --velvet-directory %s/%s --fasta %s --assembly-svg %s.svg' % (args.kmer, output_directory, output_sequence_file_path, output_sequence_file_path.replace('_assembly.fa', ''), output_sequence_file_path)
                                subprocess.check_call(cmd, shell = True)
                                
                                cmd = 'mv *.svg %s' % (output_directory)
                                subprocess.check_call(cmd, shell = True)
                                
                                cmd = 'mv %s %s' % (output_sequence_file_path.replace('_assembly.fa', ''), output_directory)
                                subprocess.check_call(cmd, shell = True)
                                
                            else:
                                SeqIO.write(sequences, tmp.name, "fasta")
                            
                                subprocess.check_call(cmd, shell = True)
                                
                                cmd = 'mv *.svg %s' % (output_directory)
                                subprocess.check_call(cmd, shell = True)
                            
                    elif assembly_type == "hmm":
                        
                        gene_size = len(list(sequences[0].seq))
                        nucleotide_composition = [{"A": [0]*gene_size}, {"T": [0]*gene_size}, {"G": [0]*gene_size}, {"C": [0]*gene_size}]
                        final_sequence = []
                        liklihood = []
                        total = [0]*gene_size
                        n_sequences = len(sequences)
                        output_file = '%s/%s' %(output_directory, output_sequence_file_path)
                        
                        for sequence in sequences:
                            sequence_list = list(sequence)
                            
                            for a, nucleotide in enumerate(sequence_list):
                                if nucleotide == "A":
                                    nucleotide_composition[0]["A"][a] += 1
                                
                                elif nucleotide == "T":
                                    nucleotide_composition[1]["T"][a] += 1
                                    
                                elif nucleotide == "G":
                                    nucleotide_composition[2]["G"][a] += 1
                                    
                                elif nucleotide == "C":
                                    nucleotide_composition[3]["C"][a] += 1
                                    
                                else:
                                    Messenger("Programming error... \n").error_message()
                                    
                                       
                        for l, position in enumerate(total):
                            
                            for nucleotide_dictionary in nucleotide_composition:
                            
                                for key, item in nucleotide_dictionary.iteritems():
                                    
                                    if item[l] > 0:
                                        item[l] = float(float(item[l])/float(n_sequences))
                                    
                                    else:
                                        continue
                        
                        
                        with open(output_file, 'w') as output:
                            
                            for x in range(0, gene_size):
                                consensus = ''
                                ids = ['A','T','G','C']
                                val = []
                                for nucleotide_dictionary in nucleotide_composition:
                                    
                                    for key, item in nucleotide_dictionary.iteritems():        
                                        val.append(item[x])
                                
                                if all(i == 0 for i in val):
                                    final_sequence.append('-')
                                    liklihood.append('-')
                                
                                else:
                                    final_sequence.append(ids[val.index(max(val))])
                                    liklihood.append(str(max(val)*100))
                            
                            output.write('>'+tax_name+'\n')
                            output.write(''.join(final_sequence)+'\n')
                            output.write(','.join(liklihood)+'\n')
        
        
        cmd = 'rm tmp*'
        subprocess.check_call(cmd, shell = True)
        
        print Messenger("Finished assembling\n").date_n_message()
        print '========================================================='
                        
                        
                        
                    



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='''--- 16S assembler %s --- Groups and assembles 16S reads from a guppy file by taxonomic rank and classification.''' % __version__
                                , epilog='==============================================================================')
    parser.add_argument('--guppy', metavar='guppy file produced by graftM', help='guppy file produced by graftM', required=True)
    parser.add_argument('--output', metavar='output directory', help='Directory within which to place assembled sequences', required=True)
    parser.add_argument('--reverse_reads', metavar='reverse reads', help='The reverse reads to extract from if necessary')
    parser.add_argument('--reads_file', metavar='reads file',  help='unaligned reads file produced by graftM', required=True)
    parser.add_argument('--kmer', metavar='k-mer',  help='k-mer for assembly with felvet')
    parser.add_argument('--assembly_type', metavar='type of assembly', help='cap3 or native hmmalignment assembly',choices = ['cap3','hmm', 'phrap', 'felvet'], required=True)
    args = parser.parse_args()
    
    if args.reverse_reads:
        if args.kmer:
            pass
        else:
            print Messenger("No k-mer specified for felvet assembly \n").error_message()
            exit(1)
            

    print intro
    
    split_guppy = TaxoGroup().guppy_splitter(args.guppy)
    grouped_sequences, r_reads = TaxoGroup().rank_grouper(split_guppy, args.reads_file)
    
    if args.reverse_reads:
        TaxoGroup().all_read_extractor(args.reverse_reads, r_reads)
        
    TaxoGroup().assembler(grouped_sequences, args.output, args.assembly_type)
            
    exit(1)    





    

    