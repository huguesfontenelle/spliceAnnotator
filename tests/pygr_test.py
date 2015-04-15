# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:20:04 2015

@author: huguesfo
"""

REFSEQ = "/Users/huguesfo/Devel/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)


from pygr import seqdb
import unittest

class TestPygr(unittest.TestCase):

    def setUp(self):
        self.chrom = '7'
        self.start = 15100200
        self.end = 15100210

    def test_simple(self):
        sequenceDB = seqdb.SequenceFileDB(REFSEQ)
        fasta = str(sequenceDB[self.chrom][self.start-1:self.end])
        self.assertEqual(fasta, 'TATTCAGCATT')
        
#=========================================================
if  __name__ == "__main__":
    unittest.main()
    
