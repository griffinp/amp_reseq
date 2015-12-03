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
#from qiime.util import create_dir

parser = argparse.ArgumentParser(description = 'Create an assembly \
file for each amplicon for each individual')

parser.add_argument('-r', '--ref_seq_fp', action='store', 
type=str, help='Reference sequence file (.fasta format)', required=True)

parser.add_argument('-c', '--consensus_dir', nargs='?', 
default='MAFFT_Consensus_Seq',action='store', type=str, 
help='Path to directory containing consensus sequences')

parser.add_argument('-l', '--collated_cons_dir', nargs='?',
default='Collated_Consensus_Seq', action='store', type=str, 
help='Path to output directory for collating consensus sequences \
per amplicon per individual')

parser.add_argument('-a', '--assembly_dir', nargs='?',
default='MAFFT_Assemblies', action='store', type=str, 
help='Path to output directory for assembled consensus sequences \
per amplicon per individual')


    
    
def add_reference_to_consensus_dir(reference_file, consensus_dir):
    """Creates a .fasta file for each reference sequence
    
    Located in the consensus directory and to be used as a backbone
    for assembly.
    """
    with open(reference_file, "r"):
        for record in SeqIO.parse(reference_file, "fasta"):
            output_name = record.id
            output_handle = open(consensus_dir+'/'+output_name+'.fasta', "w")
            SeqIO.write(record, output_handle, "fasta")
            output_handle.close()      
        
def make_ind_plus_amp_dictionary(consensus_dir):
    """Make a dictionary for each individual/amplicon combination.
    
    Key = individual/amplicon combination
    Value = list of consensus sequence files for that ind/amp combination
    example file name: TD316_cp4_506_82_3.4_cons.fasta
    """
    assembly_dict = {}
    for input_file in os.listdir(consensus_dir):
        split_name = re.split(r"_", input_file)
        if len(split_name) > 1:
            ind_name = split_name[0]
            amp_name = split_name[1]
            assembly_key = ind_name+'_'+amp_name
            #if the key doesn't exist yet, add the key:value pair
            if assembly_dict.get(assembly_key, None) == None:
                assembly_dict[assembly_key] = [consensus_dir+'/'+input_file]
            else: 
                assembly_dict[assembly_key].append(consensus_dir+'/'+input_file)
    return assembly_dict

def collate_consensus_sequences(ind_plus_amp_dict, collated_cons_dir):
    """Collates consensus sequences for a given individual/amplicon combination.
    
    Outputs collated sequences in a single .fasta file
    """
    for ind_amp in ind_plus_amp_dict.items():
        collated_file = ind_amp[0]+'_collated.fasta'
        with open(collated_cons_dir+'/'+collated_file, "w+") as output_handle:
            for consensus_file in ind_amp[1]:
                #record_name = re.split(r"/", consensus_file)[-1].strip('.fasta')
                with open(consensus_file, "r"):
                    for record in SeqIO.parse(consensus_file, "fasta"):
                        #record.id = record_name
                        SeqIO.write(record, output_handle, "fasta")

def remove_reference(assembly_file):
    """Removes the reference sequence from an alignment."""
    assembly_file_name = os.path.basename(assembly_file)
    assembly_directory = os.path.dirname(assembly_file)
    output_handle = assembly_directory + '/' + assembly_file_name[13:]
    amp_name = re.split(r"_", assembly_file_name)[3]
    with open(assembly_file, "r") as input, open(output_handle, "w+") as out:
        for record in SeqIO.parse(input, "fasta"):
            if record.id != amp_name:
                SeqIO.write(record, out, "fasta")

def main():
    args = parser.parse_args()
    base_directory = os.getcwd()
    reference_file = args.ref_seq_fp
    consensus_dir = base_directory + '/' + args.consensus_dir
    collated_cons_dir = base_directory + '/' + args.collated_cons_dir
    assembly_dir = base_directory + '/' + args.assembly_dir
    assembly_command_file = base_directory + '/' + 'mafft_assembly_commands.sh'
    add_reference_to_consensus_dir(reference_file, consensus_dir)
    ind_plus_amp_dict = make_ind_plus_amp_dictionary(consensus_dir)
    print "Collating consensus sequences from "
    print consensus_dir
    print "into "
    print collated_cons_dir
    collate_consensus_sequences(ind_plus_amp_dict, collated_cons_dir)

    print "Writing MAFFT 'assembly' command file"
    with open(assembly_command_file, "w+") as command_file:
        for input_file in os.listdir(collated_cons_dir):
            ind_name = re.split(r"_", input_file)[0]
            amp_name = re.split(r"_", input_file)[1]
            inp = str(collated_cons_dir+'/'+input_file)
            out_file = str(assembly_dir+'/'+'ref_included_' + ind_name + '_' + amp_name + '_assembled.fasta')
            ref_file = str(consensus_dir+'/'+amp_name+'.fasta')
            mafft_command = str('mafft --thread -1 --maxiterate 1000 --legacygappenalty --reorder --quiet --addfragments '+inp+' '+ref_file+' > '+out_file)
            command_file.write(mafft_command+"\n")
    
    print "Running MAFFT assembly commands..."
    mafft_process = subprocess.Popen('sh mafft_assembly_commands.sh', shell=True)
    mafft_process.wait()
    
    print "Removing reference sequence from each assembly"
    for assembly_file in os.listdir(assembly_dir):
        if assembly_file[:12] == 'ref_included':
            assembly_path = assembly_dir + '/' + assembly_file
            remove_reference(assembly_path)
    
    remove_files_process = subprocess.Popen(str('rm ' + assembly_dir + '/' + 'ref_included*'), shell=True)
    remove_files_process.wait()
    
if __name__ == '__main__':
    main()

