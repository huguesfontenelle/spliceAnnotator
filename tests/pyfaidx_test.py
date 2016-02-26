# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:20:04 2015

@author: huguesfo
"""

REFSEQ = "/Users/huguesfo/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)


from pyfaidx import Fasta
from unittest import TestCase


class TestPyfaidx(TestCase):


    def test_simple(self):
        chrom = '13'
        start, end = [32912680, 32912690]
        genome = Fasta(REFSEQ)
        fasta = genome[chrom][start:end]
        rc = fasta.reverse.complement
        
        assert fasta.seq == 'GAAGCATGTC'  
        assert fasta.orientation == 1
        assert fasta.start == start
        assert fasta.end == end
        
        assert rc.seq == 'GACATGCTTC'
        assert rc.orientation == -1

        