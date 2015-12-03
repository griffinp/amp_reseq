#!/usr/bin/env python

__author__ = "Philippa Griffin"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Philippa Griffin"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Philippa Griffin"
__email__ = "pip.griffin@gmail.com"

import os
from shutil import rmtree

from unittest import TestCase, main
import tempfile

#from qiime.util import create_dir
from qiime.process_clusters import (OtusPerInd, list_otu_map_files,
                                    identify_top_n_otu_counts,
                                    write_top_n_otus_as_fasta)


class ProcessClustersTests(TestCase):
    
    def setUp(self):
        """First defining some test data"""
        self._paths_to_clean_up = []
        self._dirs_to_clean_up = []
#        self._files_to_clean_up = []   
        
        self.amplicon_list = ['amplicon1','amplicon2','amplicon3']
        self.top_n_otus = 2
        self.otu_file_list = ['ind1_otus.txt', 'ind2_otus.txt']
        self.otu_file_list_with_repeat = ['ind1_otus.txt', 'ind2_otus.txt', 'ind1_otus.txt']
        self.blast_file_list = ['ind1_blast.txt', 'ind2_blast.txt']
        self.otu_file_list_with_repeat = ['ind1_blast.txt', 'ind2_blast.txt', 'ind1_blast.txt']
        self.fasta_file_list = ['ind1.fasta', 'ind2.fasta']
        self.fasta_file_list_with_repeat = ['ind1.fasta', 'ind2.fasta', 'ind1.fasta']
                
        self.dir_path = "/tmp/qiimetestfiles/"
        
        try:
            mkdir(self.dir_path)
        except OSError:
            pass
        
        self.otu_map_dir = dir_path + '/' + 'cdhit_picked_otus'
        self.fasta_dir = dir_path + '/' + 'Individual_Relabelled_Fasta'
        self.blast_match_dir = dir_path + '/' + 'Blast_Matches'
        self.output_dir = dir_path + '/' + 'Process_Clusters_Output'

        try:
            mkdir(self.otu_output_dir)
        except OSError:
            pass

        # define directory to clean up
        self._dirs_to_clean_up = ["/tmp/qiimetestfiles/Process_Clusters_Output"]
                    
    def tearDown(self):
        map(remove, self._paths_to_clean_up)
        map(removedirs, self._dirs_to_clean_up)
                
                
#Plan for tests to include here

    def test_list_otu_map_files_duplicate_filename(self):
        """Tests that a duplicate OTU file name raises a ValueError"""
        
        for each_file in self.otu_file_list_with_repeat:
            each_file.mkstemp(dir=self.otu_map_dir)
        self.assertRaises(ValueError, list_otu_map_files, self.otu_map_dir)




#test that expected output matches actual output
#
#test that an appropriate error is raised if the same sample name occurs twice in the Qiime_Fasta directory
#
#test that an appropriate error is raised if the same sample name occurs twice in the BLAST_Matches directory
#
#test that an appropriate error is raised if the same sample name occurs twice in the OTU_Files directory
#
#test that an appropriate error is raised if a sample name doesn't occur in all three of the input directories
#
#test that the Qiime_Fasta input format is valid
#
#test that the BLAST_Matches input format is valid
#
#test that the OTU_Files input format is valid
#
#test that an appropriate error is raised if an invalid sample ID is specified (e.g. underscore in wrong place)
#
#test that an informative error is raised if BLAST match score is missing from BLAST_Matches input (i.e. alert user that they need to use the appropriate version of blast_wrapper.py)
#
#test that the appropriate files are detected as input
#
#test that the expected output is successfully written  


blast_matches = """:     
593: amplicon1 2e-94 96.41 392 198
690: amplicon1 1e-21 86.11 211 318
340: amplicon2 9e-122 98.25 712 484
"""

otu_clusters = """593	ind1_108240
690	ind1_111111	ind1_111657	\
ind1_112681
340	ind1_114743	ind1_101893	\
ind1_102736	ind1_113925	\
"""

individual_fasta = """>ind1_101893 W30DH:00253:00731
CATAAATTTGCTTGAAAAGCCCTGGGAAAGTGCTCTCTTTTTCCATAGAAATAAATTCGTTCAATAAGGGCTCCAGA\
AGATGTTGATCGTAAGTGAGAAGATTGGTTACGGAGAAAGAGGAAGCCGGATTCATATTCACATACATGAAAAGTAT\
ATAGGAAGAAGAATAATCTCTG
>ind1_102736 W30DH:00313:00974
CATAAATTTGCTTGAAAAGCCCTGGGAAAGTGCTCTCTTTTTCCATAGAAATAAATTCGTTCAATAAGGGCTCCAGA\
AGATGTTGATCGTAAGTGAGAAGATTGGTTACGGAGAAAGAGGAAGCCGGATTCATATTCACATACATGAAAAGTATA
>ind1_108240 W30DH:00672:00406
TCCAAAGGAATGATTGGTGGGAGCGACGTGATGCGTGATGCCCAGGCGGACGTGCCCTCGGCCGAATGGCTACGGGC\
GCAACTTGCGTTCAAAAACTCGATGGTTTACAGGATTCTGCAATTCACACCAAGTATCGCATTTTGCTACGTTCTTC\
ATCGATGCGAGAGCCGAGATATCCGTTGCCGAGAGTCGTTTGTG
>ind1_111111 W30DH:00849:00888
CACTTATTAAACCCCGTTTTTAAACAACGTCTTATCAAATGTAAAAATTTGAAACAACTTTTGACAACGGATCTCTT\
GGCTCTCGCAACGATGAAGAACGCAGCGAAATGCGATACGTAGTGTGAATTGCAGAACCGTGAATCATCGAATCTTT\
GAACGCATATTGCGCCTTTGGGTATTCCCT
>ind1_111657 W30DH:00884:00390
CTTGGCTCTCGCAACGATGAAGAACGCAGCGAAATGCGATACGTAGTGTGAATTGCAGAACCGTGAATCATCGAATC\
TTTGAACGCATATTGCGCCTTTGGGTATCC
>ind1_112681 W30DH:00945:00622
CGTTTTTAAACAACGTCTTATCAAATGTAAAAATTTGAAACAACTTTTGACAACGGATCTCTTGGCTCTCGCAACGA\
TGAAGAACGCAGCGAAATGCGATACGTAGTGTGAATTGCAGAACCGTGAATCATCGAATCTTTGAACGCATATTGCG\
CCTTTGGGTATCC
>ind1_113925 W30DH:00988:01147
CATAAATTTGCTTGAAAAGCCCTGGGAAAGTGCTCTCTTTTTCCATAGAAATAAATTCGTTCAATAAGGGCTCCAGA\
AGATGTTGATCGTAAGTGAGAAGATTGGTTACGGAGAAAGAGGAAGCCGGATTCATATTCACATACATGAAAAGTAT\
ATAGGAAGAAGAATAATCTCTGATTTCTTTTTGAAAAAGAAGAACTGGCTTTCTTTGAATTTGAAGTAATAA
>ind1_114743 W30DH:01081:00473
CATAAATTTGCTTGAAAAGCCCTGGGAAAGTGCTCTCTTTTTCCATAGAAATAAATTCGTTCAATAAGGGCTCCAGA\
AGATGTTGATCGTAAGTGAGAAGATTGGTTACGGAGAAAGAGGAAGCCGGATTCATATTCACATACATGAAAAGTAT\
ATAGGAAGAAGAATAATCTCTGATTTCTTTTTGAAAAAGAAGAACTGGCTTTCTTTGAAT
"""

# run tests if called from command line
if __name__ == "__main__":
    main()