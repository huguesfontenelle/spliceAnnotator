"""
Tests for SpliceScore.

@author: Hugues Fontenelle, 2014
"""

import unittest
from splice.splice_score import SpliceScore

REFSEQGENE = "/Users/huguesfo/Documents/DATA/b37/refSeq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
GENEPANEL = "/Users/huguesfo/Devel/genevar/vcpipe-bundle/clinicalGenePanels/Bindevev_OUS_medGen_v02_b37/Bindevev_OUS_medGen_v02_b37.transcripts.csv"

class TestSpliceScore(unittest.TestCase):

    def setUp(self):
        self.data = [{'chrom': '4', 'pos':155490472, 'ref':'C', 'alt':'T'}]
        self.data_gp = [{'chrom': '17', 'pos': 41247865, 'ref': 'T', 'alt':'C'},
                        {'chrom': '13', 'pos': 32890599, 'ref': 'T', 'alt':'G'}]
        self.data_gp = [{'chrom': '17', 'pos': 48267950, 'ref': 'C', 'alt':'A'},
                        {'chrom': '14', 'pos': 74988732, 'ref': 'C', 'alt':'T'}]
                         
        self.data_outoftranscript = [{'ID': 'rs1800469', 'alt': u'G', 'chrom': u'19', 'pos': 41860296, 'ref': u'A'}]

    def test_refseqgene(self):
        '''
        Runs the scoring
        '''
        for record in self.data:
            s = SpliceScore(record)
            s.set_ref_seq(REFSEQ)
            s.set_ref_seq_gene(REFSEQGENE)
            s.use_algo(use_SSFL=True, use_MaxEntScan=True)
            s.score_splice_sites()
            self.assertTrue('wild' in s)
            self.assertTrue('mut' in s)

    def test_genepanel(self):
        for record in self.data_gp:
            s = SpliceScore(record)
            s.set_ref_seq(REFSEQ)
            s.set_gene_panel(GENEPANEL)
            s.use_algo(use_SSFL=True, use_MaxEntScan=True)
            s.score_splice_sites()
            self.assertTrue('wild' in s)
            self.assertTrue('mut' in s)
            
    def test_outoftranscript(self):
        for record in self.data_outoftranscript:
            s = SpliceScore(record)
            s.set_ref_seq(REFSEQ)
            s.set_gene_panel(GENEPANEL)
            s.use_algo(use_SSFL=True, use_MaxEntScan=True, use_GeneSplicer=True, use_NNSplice=True, use_HSF=True)
            s.score_splice_sites()
            self.assertTrue('wild' in s)
            self.assertTrue('mut' in s)
            self.assertTrue('authentic' in s)
            auth = s['authentic']
            self.assertTrue('pos' in auth)
            self.assertTrue('splice_type' in auth)
            self.assertTrue('strand' in auth)
            
# ============================================================                     
if  __name__ == "__main__":
    unittest.main()
