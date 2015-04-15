"""
Tests for the diverse splice finder algorithms.

@author: Hugues Fontenelle, 2014
"""

import unittest

from splice.SSFL import SSFL
from splice.max_ent_scan import MaxEntScan

from Bio.Seq import Seq

class TestSpliceSiteFinder(unittest.TestCase):

    """
    def setUp(self):
        pass
    """

    def test_SSFL(self):
        seq1 = "\
AAGGTAAGTggggggggggggggggggggggggggggggggggggggggg\
TTTTTTTTTTTTCAGGTggggggggggggggggggggggggggggggggg\
ACTTACCTTggggggggggggggggggggggggggggggggggggggggg\
ACCTGAAAAAAAAAAAAggggggggggggggggggggggggggggggggg".upper()
        # For testing purpose, this sequence is manufactured to show 50 nt per line
        # line 1: Donor site, + strand     aag|GTaagt          0+ 3
        # line 2: Acceptor site, + strand  tttttttttttcAG|gt  50+15
        # line 3: Donor site, - strand     acttAC|ctt        100+ 6
        # line 4: Acceptor site, - strand  ac|CTgaaaaaaaaaaa 150+ 2
        s1 = SSFL()
        s1.seq = seq1
        sites = s1.find_all_sites()
        self.assertTrue((3, 'Donor', '+') in sites)
        self.assertEqual(sites[(3, 'Donor', '+')], 100.0)
        self.assertTrue((65, 'Acceptor', '+')  in sites)
        self.assertEqual(sites[(65, 'Acceptor', '+')], 100.0)
        self.assertTrue((106, 'Donor', '-') in sites)
        self.assertEqual(sites[(106, 'Donor', '-')], 100.0)
        self.assertTrue((152, 'Acceptor','-') in sites)
        self.assertEqual(sites[(152, 'Acceptor','-')], 100.0)
        r = s1.score_donor_by_position(3)
        self.assertEqual(r.items()[0][1], 100.0)
        r = s1.score_acceptor_by_position(65)
        self.assertEqual(r.items()[0][1], 100.0)

    def test_MaxEntScan(self):
        s = MaxEntScan()
        # the following sequences and scores are computed by the website:
        # http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq_acc.html
        self.assertEqual(10.86, s.score_donor_by_sequence("cagGTaagt"))
        self.assertEqual(11.08, s.score_donor_by_sequence("gagGTaagt"))
        self.assertEqual(-0.12, s.score_donor_by_sequence("taaATaagt"))
        self.assertEqual(2.89, s.score_acceptor_by_sequence("ttccaaacgaacttttgtAGgga"))
        self.assertEqual(8.19, s.score_acceptor_by_sequence("tgtctttttctgtgtggcAGtgg"))
        self.assertEqual(-0.08, s.score_acceptor_by_sequence("ttctctcttcagacttatAGcaa"))
        # usage example
        seq1 = "\
GAGGTAAGTggggggggggggggggggggggggggggggggggggggggg\
TTTCTTTTTTTTTTTGGCAGTGTggggggggggggggggggggggggggg\
ACTTACCTCggggggggggggggggggggggggggggggggggggggggg\
ACACTGCCAAAAAAAAAAAGAAAggggggggggggggggggggggggggg".upper()
        # For testing purpose, this sequence is manufactured to show 50 nt per line
        # line 1: Donor site, + strand     gag|GTaagt                 0+ 3
        # line 2: Acceptor site, + strand  tgtctttttctgtgtggcAG|tgg  50+20
        # line 3: Donor site, - strand     acttAC|ctc               100+ 6
        # line 4: Acceptor site, - strand  cca|CTgccacacagaaaaagaca 150+ 3
        s.seq = seq1
        sites = s.find_all_sites()
        self.assertTrue((3, "Donor", '+') in sites)
        self.assertEqual(sites[(3, 'Donor', '+')], 11.08)
        self.assertTrue((70, "Acceptor", '+') in sites)
        self.assertEqual(sites[(70, "Acceptor", '+')], 8.63)
        self.assertTrue((106, "Donor", '-') in sites)
        self.assertEqual(sites[(106, "Donor", '-')], 11.08)
        self.assertTrue((153, "Acceptor",'-') in sites)
        self.assertEqual(sites[(153, "Acceptor",'-')],  8.63)



# ============================================================
if  __name__ == "__main__":
    unittest.main()
