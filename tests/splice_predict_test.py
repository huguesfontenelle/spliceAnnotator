"""
Tests for SplicePredict.

@author: Hugues Fontenelle, 2014
"""

import unittest
from splice.splice_predict import SplicePredict

REFSEQGENE = "/Users/huguesfo/Documents/DATA/b37/refSeq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
GENEPANEL = "/Users/huguesfo/Devel/genevar/vcpipe-bundle/clinicalGenePanels/Bindevev_OUS_medGen_v02_b37/Bindevev_OUS_medGen_v02_b37.transcripts.csv"

      
class TestSpliceScore(unittest.TestCase):

    def setUp(self):
        self.strategies = ['Jian14nar', 'Houdayer', 'AMG-kreftgenetikk', 'AMG-diag']
        self.data = [{'chrom': '4', 'pos':155490472, 'ref':'C', 'alt':'T',
                      'strategy': 'Houdayer'}]
        self.data_gp = [{'chrom': '17', 'pos': 48267950, 'ref': 'C', 'alt':'A',
                         'strategy': 'Houdayer'},
                        {'chrom': '14', 'pos': 74988732, 'ref': 'C', 'alt':'T',
                         'strategy': 'Houdayer'}]
        self.data_strat = [{'ID': 'rs397509306', 'chrom': '17', 'pos': 41247865, 'ref': 'T', 'alt':'C'}]
        
                 
    def test_refseqgene(self):  
        for record in self.data:
            p1 = SplicePredict(record)
            p1.set_ref_seq(REFSEQ)
            p1.set_ref_seq_gene(REFSEQGENE)
            p1.strategy = record['strategy']
            effect = p1.predict()
            self.assertTrue('wild' in p1)
            self.assertTrue('mut' in p1)    
            self.assertTrue('predict' in p1)
            self.assertTrue(record['strategy'] in p1['predict'])
            self.assertTrue('Effect' in p1['predict'][record['strategy']][0])
            self.assertEqual(effect[record['strategy']], p1['predict'][record['strategy']])

    def test_genepanel(self):
        for record in self.data_gp:
            p1 = SplicePredict(record)
            p1.set_ref_seq(REFSEQ)
            p1.set_gene_panel(GENEPANEL)
            p1.strategy = record['strategy']
            effect = p1.predict()
            self.assertTrue('wild' in p1)
            self.assertTrue('mut' in p1)  
            self.assertTrue('predict' in p1)
            self.assertTrue(record['strategy'] in p1['predict'])
            self.assertTrue('Effect' in p1['predict'][record['strategy']][0])
            self.assertEqual(effect[record['strategy']], p1['predict'][record['strategy']])       

    def test_strategies(self):
        for record in self.data_strat:
            for strategy in self.strategies:
                p1 = SplicePredict(record)
                p1.set_ref_seq(REFSEQ)
                p1.set_ref_seq_gene(REFSEQGENE)
                p1.strategy = strategy
                effect = p1.predict()
                self.assertTrue('wild' in p1)
                self.assertTrue('mut' in p1)  
                self.assertTrue('predict' in p1)
                self.assertTrue(strategy in p1['predict'])
                self.assertTrue('Effect' in p1['predict'][strategy][0])
                self.assertEqual(effect[strategy], p1['predict'][strategy])
                print p1
                print effect
    
# ============================================================                     
if  __name__ == "__main__":
    unittest.main()

'''
record={'chrom': '17', 'pos': 41247865, 'ref': 'T', 'alt':'C', 'strategy': 'AMG-kreftgenetikk'}
p1 = SplicePredict(record)
p1.set_ref_seq(REFSEQ)
p1.set_gene_panel(GENEPANEL)
p1.strategy = record['strategy']
effect = p1.predict()
'''