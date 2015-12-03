#!/usr/bin/env python

__author__ = "Philippa Griffin"
__copyright__ = ""
__credits__ = ["Philippa Griffin"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Philippa Griffin"
__email__ = "pip.griffin@gmail.com"

import argparse
import re
import os
import subprocess
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
alphabet = Gapped(IUPAC.ambiguous_dna)
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
from Bio import SeqRecord
#from qiime.util import create_dir

parser = argparse.ArgumentParser(description = 'Create an alignment \
and consensus sequence from each OTU fasta file')

parser.add_argument('-f', '--otu_fasta_dir', action='store', 
type=str, help='Path to directory containing OTU fasta files', required=True)

parser.add_argument('-a', '--align_output_dir', nargs='?', 
default='MAFFT_Alignments',action='store', type=str, 
help='Path to output directory for alignments')

parser.add_argument('-c', '--align_cons_dir', nargs='?',
default='MAFFT_Consensus_Seq', action='store', type=str, 
help='Path to output directory for consensus sequences')

parser.add_argument('-s', '--split_large_fasta', type=str, action='store',
choices=['True', 'False'], default='True',
help='Turn off splitting of large fasta files. Default is on')

parser.add_argument('-t', '--similarity_threshold', type=float, action='store',
nargs='?', default=0.2,
help='Similarity threshold for sequence consensus estimation. Default is 0.2')


# Testing
#args = parser.parse_args(['-f, './Process_OTUs_Output', '-a', './MAFFT_Alignments', '-c', './MAFFT_Consensus_Seq', '-s', 'True', '-t', '0.2'])


def file_len(file):
    with open(file) as f:
        no_fasta = 0
        for line in f:
            no_fasta = no_fasta + line.count('>')
    return no_fasta + 1

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    
    from http://biopython.org/wiki/Split_large_file
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def input_seq_iterator(fasta_file):
    """Splits fasta files containing over 5000 sequences
    """
    input_seq_iterator = SeqIO.parse(open(fasta_file), "fasta")
    for i, batch in enumerate(batch_iterator(input_seq_iterator, 5000)):
        filename = fasta_file + '_subset_%i.fasta' % (i + 1)
        output_handle = open(filename, "w")
        count = SeqIO.write(batch, output_handle, "fasta")
        output_handle.close()
    print "Wrote %i records to %s" % (count, filename)

def fasta_subset(fasta_directory, base_directory):
    """Identifies fasta files containing over 5000 sequences
    """
    if os.path.exists("./Large_OTU_Fasta_Files"):
        pass
    else:
        os.mkdir("./Large_OTU_Fasta_Files")
    os.chdir(fasta_directory)
    for file in os.listdir("."):
        if file_len(file) < 5001:
            pass
        else:
            input_seq_iterator(file)
            os.rename(file, str(base_directory + '/Large_OTU_Fasta_Files/' + file))
            
            
def list_alignment_files(alignment_dir):
    """Makes a list of all alignment files in the MAFFT Output directory
    i.e. files that end in '_aligned.fasta' """
    aligned_list = []
    for file in os.listdir(alignment_dir):
        if file.endswith("_aligned.fasta"):
            aligned_list.append(file)
    return aligned_list
       
def output_consensus(y, threshold_value, consensus_output_dir):
    """Takes as input an alignment file
    and outputs a consensus sequence in fasta format"""
    file_name = os.path.basename(y)
    fasta_name = file_name.split('_align')[0]
    alignment = AlignIO.read(open(y), "fasta")
    summary_align = SummaryInfo(alignment)
    consensus = summary_align.gap_consensus(threshold = threshold_value, 
                                            ambiguous = 'N', 
                                            consensus_alpha = alphabet, 
                                            require_multiple = 1)
    consensus_seq = SeqRecord.SeqRecord(consensus,id=fasta_name+"_consensus")
    output_file_name = str(consensus_output_dir+'/'+fasta_name+"_cons.fasta")
    output_handle = open(output_file_name, "w")
    print "Writing consensus sequence for " + fasta_name
    SeqIO.write(consensus_seq, output_handle, "fasta")
    output_handle.close()



def main():
    args = parser.parse_args()
    base_directory = os.getcwd()
    otu_fasta_dir = base_directory + '/' + args.otu_fasta_dir
    align_output_dir = base_directory + '/' + args.align_output_dir
    consensus_output_dir = base_directory + '/' + args.align_cons_dir
    if not os.path.exists(consensus_output_dir):
        os.makedirs(consensus_output_dir)
    if args.split_large_fasta == 'True':
        print "Splitting fasta files with > 5000 sequences"
        fasta_subset(otu_fasta_dir, base_directory)
        
    output_command_file = str(base_directory + '/mafft_alignment_commands.sh')
    
    print "Writing MAFFT alignment command file"
    with open(output_command_file, "w+") as command_file:
        for input_file in os.listdir(otu_fasta_dir):
            ind_name = re.split(r"\.fa", str(input_file))[0]
            out_file = str(align_output_dir+'/'+ind_name+'_aligned.fasta')
            mafft_command = str('mafft --thread -1 --quiet '+otu_fasta_dir+'/'+input_file+' > '+out_file)
            command_file.write(mafft_command+"\n")
    
    os.chmod(output_command_file, 0755)
    print "Running MAFFT alignment commands..."
    mafft_process = subprocess.Popen(output_command_file, shell=True)
    mafft_process.wait()

    print "Removing empty alignment files (produced from only 1 input sequence)"
    find_command = str('find ' + align_output_dir + ' -size 0 -delete -type f')
    find_empty_files = subprocess.Popen(find_command, shell=True)
    find_empty_files.wait()

    aligned_list = list_alignment_files(align_output_dir)
    for file in aligned_list[0:-1]:
        output_consensus(str(align_output_dir + '/' + file), args.similarity_threshold, consensus_output_dir)
        
if __name__ == '__main__':
    main()
