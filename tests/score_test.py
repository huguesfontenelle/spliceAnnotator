"""
Tests for the diverse splice finder algorithms.

@author: Hugues Fontenelle, 2014
"""

from unittest import TestCase
from splice import max_ent_scan as mes

class TestSpliceSiteFinder(TestCase):

    def test_MaxEntScan(self):
        # the following sequences and scores are computed by the website:
        # http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq_acc.html
        assert mes.score("cagGTaagt") == [10.86]
        assert mes.score("gagGTaagt") == [11.08]
        assert mes.score("taaATaagt") == [0.0]
        assert mes.score("ttccaaacgaacttttgtAGgga") == [2.89]
        assert mes.score("tgtctttttctgtgtggcAGtgg") == [8.19]
        assert mes.score("ttctctcttcagacttatAGcaa") == [0.0]
