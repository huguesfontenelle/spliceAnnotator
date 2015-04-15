"""
GeneSplicer

python wrapper for the GeneSplicer of Mihaela Pertea
See Copyright and README notice of the original software below.

Implementation of:
Pertea et al. "GeneSplicer: a new computational method for splice site 
prediction." Nucleic Acids Res (2001) vol. 29 (5) pp. 1185-90

@author: Hugues Fontenelle, 2014

------------------------------------------------
Copyright @ 2001, The Institute for Genomic Research (TIGR). All rights
reserved.

This system is licensed as open source under the terms in the file LICENSE.


GeneSplicer is a program trained to find splice sites in Arabidopsis 
Thaliana. 
USAGE: genesplicer <fasta-file> <specific-genome-training-directory> [options]
Options:
-f file_name     Write the results in file_name
-a t             Choose t as a threshold for the acceptor sites
-d t             Choose t as a threshold for the donor sites
-e n             The maximum acceptor score within n bp is chosen
-i n             The maximum donor score within n bp is chosen
-h               Display these options of the program
Output:
End5 End3 Score "confidence" splice_site_type
E.g. :
202 203 2.994010 Medium donor

Mihaela Pertea
mpertea@tigr.org
------------------------------------------------
"""

import os.path, inspect

from splice.splice_base import SpliceBase

REL_PATH_TO_3RDPARTY = '../thirdparty/genesplicer/'


class GeneSplicer( SpliceBase ):
    
    # ------------------------------------------------
    def __init__(self, seq="", offset=0):
        super(GeneSplicer, self).__init__(seq, offset)
    
    # ------------------------------------------------
    def find_all_sites(self, splice_type='all', strand='all'):  
        ori_path = os.getcwd()
        path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        os.chdir(path)
        path_to_binaries = path + '/' + REL_PATH_TO_3RDPARTY + "bin/mac/"
        path_to_specific_genome_training_directory = path + '/' + REL_PATH_TO_3RDPARTY + "human/"
        splicer_software = "genesplicer"
        
        input_filename = path + "/" + "genesplicer.fasta"
        output_filename = path + "/" + "genesplicer.output"

        with open(input_filename, "w") as f:
            f.write('>seq\n')
            f.write(str(self.seq))

        # $ genesplicer
        # USAGE:  genesplicer <fasta-file> <specific-genome-training-directory> [options] 
        cmd = [path_to_binaries+splicer_software, input_filename, 
               path_to_specific_genome_training_directory, "-f", output_filename]
        os.system(" ".join(cmd))
        
        sites = {}
        splice_type_d = {'donor':'Donor', 'acceptor':'Acceptor'}
        offset_d = {('Acceptor','+') : 2 + self.offset,
                    ('Acceptor','-') : 0 + self.offset,
                    ('Donor','+') : 0 + self.offset,
                    ('Donor','-') : 2 + self.offset}
                    
        with open(output_filename, "r") as splice_file:
            for line in splice_file:
                l = line.split()
                mi = min(int(l[0]),int(l[1]))
                #ma = max(int(l[0]),int(l[1]))
                if int(l[0])<int(l[1]): # forward case
                    local_strand = "+"
                else: # reverse case
                    local_strand = "-" 
                local_splice_type = splice_type_d[l[4]]
                score = float(l[2])
                sites.update({ (mi + offset_d[(local_splice_type, local_strand)], local_splice_type, local_strand) : score })

        
        os.chdir(ori_path)
        
        if splice_type in ['Donor', 'Acceptor', 'all', None]:
            sites = self.filter_sites(sites, splice_type=splice_type)
        if strand in ['+', '-', 'all', None]:
            sites = self.filter_sites(sites, strand=strand)        

        return sites
    
    # ------------------------------------------------
    def export_to_bed(self, filename, geneID, sites):
        super(GeneSplicer, self).export_to_bed(filename, geneID, sites, "GeneSplicer")

        
# ================================================
# main
# ================================================
if  __name__ == "__main__":
    # usage example
    s1 = GeneSplicer()
    seq1 = "\
ggtatgtaaaaagtaatataacagaaaaataaatatcttgtctgaatcagagagagagctctcaacaaccaaagctttggactctcaggccaggaagact\
atgggattgggcggcgatcaatcctttgttcccgtcatggattccggccaagtccgcctcaaggagctcggctacaagcaagagctcaagcgcgatctct\
cgtctgtctccgccactctttcctttttcttctttgatttgtgttagctaaacaattgttagatttcagggtcttctccaatttcgccatctccttctcc\
atcatatcggtgctcactggtatcaccaccacctacaacaccggcttaagattcggcggcactgtcactctggtctacggatggttcctggccggctcct\
tcacaatgtgcgttgggttatctatggccgagatctgctcctcttaccctacctccggtggtctctactactggagtgctatgctcgctggccctcgttg\
ggctcctcttgcctcttggatgactggctggtactcccttctctcacctttttcttaatttggaattggaagcaagagagcgtctctctctgactcattt\
tcacaaatcttttttcaggttcaacatcgttggtcaggtgcaacattcattcaagcagacttgttttcattagatatatctctcattggggctaaaagta\
aggttgtggtgtgtttgcagtgggcagtgacggccagcgttgacttctctctggcacagttgattcaggtgatcgtccttctctccaccggcggtagaaa\
cggcggcggttataaaggatcagactttgttgtgattggtatccatggtgggatcctcttcatccacgctcttctcaacagcctccccatctccgtcttg\
tctttcattggacagcttgctgctctttggaatctcctcggtacccaatctctcttttgtctctgtgtttcatgtcatgtctcaacaacaatcatcaatc\
aagcatctcttatattcaggggttttggtgctcatgattctgattcctttggtttctacggaaagagcaaccactaagtttgtctttaccaatttcaaca\
ctgataatggccttggcatcaccagctacgcttacatattcgttttgggactcctcatgagccagtacaccattacagggtatgatgcctctgcccacat\
ggtgagtgtaatggattctgtttttaatgattgattagcttttcttgaggattgaacaaaaagtaaatttgcagacagaagagacagtcgacgcagacaa\
gaacgggcccagaggaataatcagtgcaattggtatatcaattctgtttggatggggttatatattgggcataagctatgccgtcacagacataccttct\
cttctgagtgagaccaacaactctggtgggtatgccattgctgagatcttctacttagctttcaagaataggtttgggagcggtactggtggaatcgtgt\
gcttaggcgttgttgcggttgctgtgtttttctgtggcatgagctctgtcaccagcaattccaggtatatatatacattcatgtttggttaaaacatatt\
ctgtctgcatttaggacttgatcatatgtgtgtttggttaatcaggatggcgtatgcgttttcgagagatggagcgatgccaatgtcgccgttatggcac\
aaagtgaacagcagagaggtccccattaatgcggtttggctctctgctctcatatcattttgtatggccctcaccgtaagtccccttctctctttctgtg\
tatctctcttttgtgtgttgggattgctgattgatgagatgaaactggcagtcactggggagcatagtggcgttccaagcaatggtgtcgatcgcaacga\
ttggactgtacatagcatacgcaatcccgattatattgagagtgacgcttgcgcgcaacaccttcgtacctggaccattcagcctgggaaaatacggaat\
ggtagtcgggtgggtggcagtcctgtgggttgtaaccatatcagtcctcttctccttacccgtggcatatcccataacagcagagacactcaactacact\
ccggtggcagttgctggtttggtggccataaccctctcatattggcttttcagtgcccgccactggtttacgggtcccatctccaacattcttagctgaa\
acaacgtttcaggggttgatttatagtaaagccattgttagctagctgatatttacccttttatacattcaacaagcaacaacacatctgtcacatttta".upper()
    s1.seq = seq1
    sites = s1.find_all_sites()
    #s1.export_to_bed("test.bed","myID",sites)
