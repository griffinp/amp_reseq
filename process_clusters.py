#!/usr/bin/env python

__author__ = "Philippa Griffin"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Philippa Griffin"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Philippa Griffin"
__email__ = "pip.griffin@gmail.com"

import argparse
import re
import os
import heapq
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
alphabet = Gapped(IUPAC.ambiguous_dna)
#from qiime.util import create_dir

parser = argparse.ArgumentParser(description = 'Process OTU Map, \
BLAST Match and .fasta files.')

parser.add_argument('-m', '--otu_map_dir', action='store', 
type=str, help='Directory containing OTU map files', required=True)

parser.add_argument('-f', '--fasta_dir', action='store', 
type=str, help='Directory containing .fasta files: one file for each individual', required=True)

parser.add_argument('-b', '--blast_match_dir', action='store', 
type=str, help='Directory containing blast match files', required=True)

parser.add_argument('-o', '--output_dir', nargs='?', 
default='Process_OTUs_Output/',action='store', type=str, 
help='Optional: output directory.')

parser.add_argument('-n', '--top_n_otus', nargs='?', default=10,
action='store', type=int, help='Choosing the top "n" most abundant OTUs for analysis per amplicon per individual. Should be a positive integer.')

parser.add_argument('-a', '--amplicon_name', type=str, action='append',
help='Amplicon name. Must be written exactly as in the BLAST reference.')

args = parser.parse_args()
base_directory = os.getcwd()
otu_map_dir = base_directory + '/' + args.otu_map_dir
#check whether called 'individual fasta dir' or 'fasta dir'
fasta_dir = base_directory + '/' + args.fasta_dir
blast_match_dir = base_directory + '/' + args.blast_match_dir
output_dir = base_directory + '/' + args.output_dir
top_n_otus = args.top_n_otus


class OtusPerInd:
    """Data class that collates blast matches and sequences.
    
    """
    def __init__(self, filepath, otu_map_dir, fasta_dir, blast_match_dir):
        self.name = filepath.split('/')[-1]
        self.info = {}
        self.info['ind_name'] = str(self.name).split('_')[0]
        self.info['otu_file'] = otu_map_dir + '/' + str(self.name)
        self.info['seq_file'] = fasta_dir + '/' + self.info['ind_name'] \
                                + '.fasta'
        self.info['blast_file'] = blast_match_dir + '/' + self.info['ind_name'] \
                                  + '_blast.txt'
        self.info['blast_match_dict'] = self.make_blast_match_dictionary()
    
    def __str__(self):
        return 'OtusPerInd object ' + self.get_ind_name()
    
    def __repr__(self):
        return str(self.info)
    
    def get_ind_name(self):
        """ Gets the individual name from the OTU map file name
        """
        return self.info['ind_name']        
    def get_otu_fp(self):
        """ Gets the OTU map file
        """
        return self.info['otu_file']        
    def get_seq_fp(self):
        """ Gets the fasta sequence file
        """
        return self.info['seq_file']        
    def get_blast_fp(self):
        """Gets the BLAST match file
        """
        return self.info['blast_file']                         
    
    def make_blast_match_dictionary(self):
        """Makes a dictionary entry of blast match details.
        
        This function makes a dictionary for this individual
        where each key is a cluster name and the corresponding value 
        is a list, containing the blast match (i.e. the amplicon 
        matched by that cluster), the match direction 
        (forward or reverse), the no. reads in that cluster, 
        the identities of the reads in that cluster (as a list) 
        i.e. the format is {cluster name : [blast_match, match_dirn, 
        no_reads, [read_identities]]}
        It outputs the entire dictionary.
        """
        blast_match_dict = {}
        with open(self.get_blast_fp(), "r") as opened_blast_file:
            for line in opened_blast_file:
                split_line = re.split(r"[\s\t:]+", line)
                if len(split_line) == 7:
                    #nb counts the \n as a split point
                    otu_name = split_line[0]
                    blast_match = split_line[1]
                    blast_match_score = float(split_line[2])
                    if blast_match_score <= 1e-5:
                        blast_match_dict[otu_name] = [blast_match]
                        #Could include a test here to check that 
                        #the values of blast_match in blast_match_dict 
                        #are included in the amplicon dictionary created 
                        # from the command line input
                        for_or_rev = int(split_line[-2]) - int(split_line[-3])
                        if for_or_rev > 0:
                            direction = "forward"
                        else:
                            direction = "reverse"
                        blast_match_dict[otu_name].append(direction)
        with open(self.get_otu_fp(), "r") as opened_otu_file:
            for line in opened_otu_file:
                otu_name = re.split(r"[\s\t]+", line)[0]
                line_length = len(re.split(r"[\s\t]+", line)[1:-1])
                read_names = re.split(r"[\s\t]+", line)[1:-1]
                if otu_name in blast_match_dict:
                    blast_match_dict[otu_name].append(line_length)
                    blast_match_dict[otu_name].append(read_names)
                else:
                    blast_match_dict[otu_name] = ['no match', 'no direction', 
                                                  line_length, read_names]
        return blast_match_dict

def list_otu_map_files(otu_map_dir):
    """Makes a list of all OTU map files in the OTU Maps directory

    i.e. files that end in '.txt'
    """
    assert os.path.isdir(otu_map_dir), "OTU map directory path does not exist"
    otu_table_list = []
    for file in os.listdir(otu_map_dir):
        if file.endswith(".txt"):
            otu_table_list.append(otu_map_dir +'/' +  file)
    assert otu_table_list, "Files in OTU map directory must end in '.txt' \
                            Was the directory path correctly specified?"
    try:
        len(otu_table_list) == len(set(otu_table_list))
    except(ValueError):
        print(len(otu_table_list), " files in OTU Maps Directory but only \
              ", len(set(otu_table_list)), " unique file names.")
    return otu_table_list

def identify_top_n_otu_counts(x, amplicon_list, top_n_otus):
    """Takes an object of the class OtusPerInd and outputs a dictionary
    
    Output dictionary entries are of the form 
    {amplicon_name : [read count of the top 'n'th otu, total read count]} 
    One entry is made for each amplicon for the given individual.
    """
    assert isinstance(x, OtusPerInd), "Was the OtusPerInd object made \
    incorrectly?"
    # Use the amplicon names ( which were input on the command line) 
    # as keys for a new dictionary.
    assert amplicon_list, "Amplicon list is empty. Were the amplicon \
    names entered in the input command?"
    amplicon_top_n_dictionary = dict.fromkeys(amplicon_list, [])
    blast_match_dict = x.info['blast_match_dict']
     
    # For each amplicon - ie each key in the amplicon_count_directory - 
    # find all entries in the blast_match_dict that have this amplicon 
    # in their value (first element in the list)
    ###### NB need to include action if the amplicon doesn't appear in the 
    # blast_match_dict
    # could make this a new 'subset' dictionary
    # then, within this subset dictionary, make a new 'top_n' dictionary
    # that contains only entries with the top n numbers 
    # (second element in the value list)
 
    for i in amplicon_top_n_dictionary.keys():                             
        amp_match_dict = {otu_name: otu_value_list for otu_name, 
                         otu_value_list in blast_match_dict.items() 
                         if otu_value_list[0] == i}
        if amp_match_dict:
            read_count_list = []
            for amp_value_list in amp_match_dict.values():
                read_count = amp_value_list[2]
                read_count_list.append(read_count)
            min_count = min(heapq.nlargest(top_n_otus, read_count_list))
            read_sum = sum(read_count_list)
            amplicon_top_n_dictionary[i] = [min_count, read_sum]
            
    return amplicon_top_n_dictionary

def write_top_n_otus_as_fasta(x, amplicon_list, top_n_otus, output_dir):
    """Outputs fasta file of raw sequences for each of the top N OTUs.
    
    This function takes an object of the class OtusPerInd as input 
    and outputs the sequences in .fasta format for each of the top n OTUs 
    for that individual.
    It requires the functions 'identify_top_n_otu_counts' 
    and 'list_blast_matches'
    """

    temp_blast_match_dict = x.make_blast_match_dictionary()
    temp_ind_name = x.get_ind_name()
    amplicon_top_n_dictionary = identify_top_n_otu_counts(x, amplicon_list, top_n_otus)
                                   
    for amplicon_name in amplicon_top_n_dictionary.keys():
        # Only proceed for the otus with read counts in the top n per amplicon
        read_value = amplicon_top_n_dictionary[amplicon_name]
        if read_value:
            read_sum = read_value[1]
            amplicon_top_reads_dictionary = {otu_name: otu_value_list 
                                             for otu_name, otu_value_list 
                                             in temp_blast_match_dict.items() 
                                             if len(otu_value_list) == 4 
                                             and otu_value_list[0] == amplicon_name 
                                             and otu_value_list[2] >= 
            #not sure how to do spacing on following line
                             amplicon_top_n_dictionary[amplicon_name][0]
                                             }
            for otu_name in amplicon_top_reads_dictionary.keys():
                direction = amplicon_top_reads_dictionary[otu_name][1]
                percent_total_reads = \
                amplicon_top_reads_dictionary[otu_name][2]/float(read_sum)*100
                pc_total_reads = str(round(percent_total_reads, 2))
                number_reads = str(amplicon_top_reads_dictionary[otu_name][2])
                reads_for_output = amplicon_top_reads_dictionary[otu_name][3]
                
                #name the output fasta file
                output_fasta = str(output_dir + '/' + x.info['ind_name'] 
                                  + '_' + amplicon_name + '_' + otu_name 
                                  + '_' + number_reads + '_' + pc_total_reads 
                                  + '.fasta'
                                  )
                
            # finally write all sequences from fasta_sequences that match
            # an id occurring in line_as_list to an output file
                
                with open(output_fasta, "w") as output_file:
                    fasta_sequences = SeqIO.parse(open(x.info['seq_file']), 
                                                 'fasta')
                    for record in fasta_sequences:
                        # check if sequence needs reverse complementing, 
                        # and do so if needed
                        if record.id in reads_for_output:
                            if direction == 'reverse':
                                #temp_dna = Seq(record.seq, generic_dna)
                                temp_dna = record.seq
                                record.seq = temp_dna.reverse_complement()
                            #print ind_name, otu_name, seq.id
                            SeqIO.write([record], output_file, "fasta")
    #Would be nice to have this string show only the amplicons 
    #actually found in that individual   
    print "Top " + str(top_n_otus) + " clusters per amplicon written to directory " \
    + output_dir + " for each of the amplicons " + \
    str(amplicon_list) + " for individual " + temp_ind_name

def main():
    """Processes OTUs, outputs a consensus sequence for the top N OTUs.

    For each file in the directory otu_map_dir, it will find the 
    matching BLAST_matches file (in blast_match_dir) and .fasta file 
    (in fasta_dir). Then it will output (in output_dir) one file for 
    each of the top n OTUs by read count, for each amplicon.

    """
    #args = parser.parse_args()
    #base_directory = os.getcwd()
    #otu_map_dir = base_directory + '/' + args.otu_map_dir
    ##check whether called 'individual fasta dir' or 'fasta dir'
    #fasta_dir = base_directory + '/' + args.fasta_dir
    #blast_match_dir = base_directory + '/' + args.blast_match_dir
    #output_dir = base_directory + '/' + args.output_dir
    #top_n_otus = args.top_n_otus
    
    # Input tests
    if not os.path.isdir(otu_map_dir):
        raise IOError("Specified OTU Map Directory" +
                      "%s does not exist" % otu_map_dir)
    if not os.path.isdir(fasta_dir):
        raise IOError("Specified Fasta Directory" +
                      "%s does not exist" % fasta_dir)
    if not os.path.isdir(blast_match_dir):
        raise IOError("Specified BLAST Match Directory" +
                      "%s does not exist" % blast_match_dir)   
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    try:
        isinstance(top_n_otus, int)
    except ValueError:
        print("'top_n_otus' must be a positive integer" +
              "(%i was specified)" % top_n_otus)
    try:
        top_n_otus >= 1
    except ValueError:
        print("'top_n_otus' must be a positive integer" +
              "(%i was specified)" % top_n_otus)

    number_amplicons = len(args.amplicon_name)
    amplicon_list = args.amplicon_name   
    otu_file_list = list_otu_map_files(args.otu_map_dir)
    for otu_file in otu_file_list:
        otus_per_ind = OtusPerInd(otu_file, otu_map_dir, fasta_dir, blast_match_dir)
        print 'Processing ' + str(otus_per_ind)
        write_top_n_otus_as_fasta(otus_per_ind, amplicon_list, top_n_otus, output_dir)

if __name__ == '__main__':
    main()

        
#testing

#otu_file_list = ['/Users/pgriffin/Documents/Students/Liz James/Raw_Data_II/cdhit_picked_otus/_otus.txt']

###### testing ########

#otutest = OtusPerInd(otu_map_dir + '/EJP01_otus.txt')
