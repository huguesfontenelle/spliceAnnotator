# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:20:04 2015

@author: huguesfo
"""

from pyfaidx import Fasta
from unittest import TestCase
import os

DATA_DIR =  os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '../data/'
)

REFSEQ = os.path.join(DATA_DIR, '../data/human_g1k_v37.fasta')

class TestPyfaidx(TestCase):

    def test_simple(self):
        chrom = '13'
        start, end = [32912680, 32912690]
        genome = Fasta(REFSEQ)
        fasta = genome[chrom][start:end]
        rc = fasta.reverse.complement

        assert fasta.seq == 'GAAGCATGTC'
        assert fasta.orientation == 1
        assert fasta.start == start + 1
        assert fasta.end == end

        assert rc.seq == 'GACATGCTTC'
        assert rc.orientation == -1
