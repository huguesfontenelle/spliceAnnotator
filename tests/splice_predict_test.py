"""
Tests for SplicePredict.

@author: Hugues Fontenelle, 2014
"""

import unittest
from splice.splice_predict import SplicePredict

REFSEQGENE = "/Users/huguesfo/genevar/vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
GENEPANEL = "/Users/huguesfo/genevar/vcpipe-bundle/clinicalGenePanels/Bindevev_v02/Bindevev_v02.transcripts.csv"


class TestSplicePredict(unittest.TestCase):

    def setUp(self):
        self.strategies = ['Jian14nar', 'Houdayer', 'AMG-kreftgenetikk', 'AMG-diag']
        self.data = [{'chrom': '4', 'pos':155490472, 'ref':'C', 'alt':'T',
                      'strategy': 'Houdayer'}]
        self.data_gp = [{'chrom': '17', 'pos': 48267950, 'ref': 'C', 'alt':'A',
                         'strategy': 'Houdayer'},
                        {'chrom': '14', 'pos': 74988732, 'ref': 'C', 'alt':'T',
                         'strategy': 'Houdayer'}]
        self.data_strat = [{'ID': 'rs397509306', 'chrom': '17', 'pos': 41247865, 'ref': 'T', 'alt':'C'}]
        self.data_indels = [
            {'chrom': '3', 'pos': 30686236, 'ref': 'CA', 'alt': 'AT',
                'ID': '0', 'strategy': 'Jian14nar'},
            {'chrom': '3', 'pos': 30686405, 'ref': 'ATG', 'alt': 'TTC',
                'ID': '0', 'strategy': 'Jian14nar'}]
        self.data_denovo = [
            {'chrom': '3', 'pos': 30686242, 'ref': 'A', 'alt': 'G',
                'ID': '0', 'strategy': 'Jian14nar'},
            {'chrom': '3', 'pos': 30686413, 'ref': 'C', 'alt': 'T',
                'ID': '0', 'strategy': 'Jian14nar'},
            {'chrom': '3', 'pos': 30686258, 'ref': 'C', 'alt': 'T',
                'ID': '0', 'strategy': 'Jian14nar'}]

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

    def test_indels(self):
        for record in self.data_indels:
            p1 = SplicePredict(record)
            p1.set_ref_seq(REFSEQ)
            p1.set_gene_panel(GENEPANEL)
            p1.strategy = record['strategy']
            effect = p1.predict()
            self.assertTrue('wild' in p1)
            self.assertTrue('mut' in p1)
            self.assertTrue('predict' in p1)
            self.assertTrue('Effect' in p1['predict'][record['strategy']][0])
            self.assertEqual(effect[record['strategy']], p1['predict'][record['strategy']])
            print p1
            print effect

    def test_denovo(self):
        for record in self.data_denovo:
            p1 = SplicePredict(record)
            p1.set_ref_seq(REFSEQ)
            p1.set_gene_panel(GENEPANEL)
            p1.strategy = record['strategy']
            effect = p1.predict()
            self.assertTrue('wild' in p1)
            self.assertTrue('mut' in p1)
            self.assertTrue('predict' in p1)
            self.assertTrue('Effect' in p1['predict'][record['strategy']][0])
            self.assertEqual(effect[record['strategy']], p1['predict'][record['strategy']])
            print p1
            print effect

# ============================================================
if  __name__ == "__main__":
    unittest.main()
