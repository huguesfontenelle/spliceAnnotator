# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:20:04 2015

@author: huguesfo
"""

REFSEQ = "/Users/huguesfo/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)


from pyfaidx import Fasta
import unittest


class TestPyfaidx(unittest.TestCase):

    def setUp(self):
        self.chrom = '13'
        self.start = 32912680
        self.end = 32912690
        self.genome = Fasta(REFSEQ)

    def test_simple(self):
        s = self.genome[self.chrom][self.start-1:self.end]
        self.assertEqual(str(s), 'AGAAGCATGTC')

    ''''
    def test_orientation(self):
        s = self.genome[self.chrom][self.start-1:self.end]
        sr = -s
        sr2 = s.reverse
        self.assertEqual(s.orientation, 1)
        self.assertEqual(sr.orientation, -1)
        self.assertEqual(sr2.orientation, -1)
    '''

    def test_complement(self):
        s = self.genome[self.chrom][self.start-1:self.end]
        scomp = s.complement
        self.assertEqual(scomp.orientation, -1) # Returns -1 (is this correct?)

    def test_start(self):
        s = self.genome[self.chrom][self.start-1:self.end]
        self.assertEqual(s.start, self.start)

    def test_end(self):
        s = self.genome[self.chrom][self.start-1:self.end]
        self.assertEqual(s.end, self.end)

#=========================================================
if  __name__ == "__main__":
    unittest.main()
