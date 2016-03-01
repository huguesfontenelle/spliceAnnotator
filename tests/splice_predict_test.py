"""
Tests for SplicePredict.

@author: Hugues Fontenelle, 2014
"""

import unittest
from unittest import TestCase, SkipTest
from splice import splice_predict as p

REFSEQGENE = "/Users/huguesfo/genevar/vcpipe-bundle/funcAnnot/refseq/refGene_150311.tsv" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
GENEPANEL = "/Users/huguesfo/genevar/vcpipe-bundle/clinicalGenePanels/Bindevev_v02/Bindevev_v02.transcripts.csv"


class TestSplicePredict(TestCase):

    def setUp(self):
        self.data = [
            {'chrom': '4', 'pos':155490472, 'ref':'C', 'alt':'T'},
            {'chrom': '17', 'pos': 48267950, 'ref': 'C', 'alt':'A'},
            {'chrom': '14', 'pos': 74988732, 'ref': 'C', 'alt':'T'},
            {'chrom': '17', 'pos': 41247865, 'ref': 'T', 'alt':'C'},
            {'chrom': '3', 'pos': 30686236, 'ref': 'CA', 'alt': 'AT'},
            {'chrom': '3', 'pos': 30686405, 'ref': 'ATG', 'alt': 'TTC'},
            {'chrom': '3', 'pos': 30686242, 'ref': 'A', 'alt': 'G'},
            {'chrom': '3', 'pos': 30686413, 'ref': 'C', 'alt': 'T'},
            {'chrom': '3', 'pos': 30686258, 'ref': 'C', 'alt': 'T'}
            ]
        self.keys =  ['auth_pos',
                      'distance',
                      'effect_descr',
                      'mut_score',
                      'mut_seq',
                      'splice_type',
                      'strand',
                      'transcript',
                      'wild_score',
                      'wild_seq']
        self.sequences = [
            {'chrom': '3', 'pos':101520831, 'ref':'G', 'alt':'',
             'wild_seq': 'GAGGTATGT', 'mut_seq': 'AAGGTATGT'}, # del on left of 5'ss+
            {'chrom': '3', 'pos':101520836, 'ref':'A', 'alt':'',
             'wild_seq': 'GAGGTATGT', 'mut_seq': 'GAGGTTGTA'}, # del on right of 5'ss+
            {'chrom': '3', 'pos':101520833, 'ref':'GGT', 'alt':'',
             'wild_seq': 'GAGGTATGT', 'mut_seq': 'GAATGTACT'}, # del start on left, finishes on right of 5'ss+: in this case the site is destroyed..
            {'chrom': '3', 'pos':101520832, 'ref':'A', 'alt':'AA',
             'wild_seq': 'GAGGTATGT', 'mut_seq': 'AAGGTATGT'}, # ins on left of 5'ss+
            {'chrom': '3', 'pos':101520836, 'ref':'A', 'alt':'AAA',
             'wild_seq': 'GAGGTATGT', 'mut_seq': 'GAGGTAAAT'}, # ins on right of 5'ss+
            {'chrom': '3', 'pos':101520832, 'ref':'A', 'alt':'G',
             'wild_seq': 'GAGGTATGT', 'mut_seq': 'GGGGTATGT'}, # SNP on left of 5'ss+
            {'chrom': '3', 'pos':101520836, 'ref':'A', 'alt':'C',
             'wild_seq': 'GAGGTATGT', 'mut_seq': 'GAGGTCTGT'}, # SNP on right of 5'ss+
            {'chrom': '2', 'pos':162060105, 'ref':'A', 'alt':'',
            'wild_seq': 'AAGGTATGT', 'mut_seq': 'CAGGTATGT', 'mut_score':9.8} # GIN-564
        ]
        self.genepanel = {
            'genepanel': "/Users/huguesfo/genevar/vcpipe-bundle/clinicalGenePanels/Bindevev_v02/Bindevev_v02.transcripts.csv",
            'data': [
                {'chrom': '1', 'pos':12017900, 'ref':'T', 'alt':'A', 'transcript':'NM_000302.3'},
                {'chrom': '2', 'pos':189859061, 'ref':'A', 'alt':'AAG', 'transcript':'NM_000090.3'},
                {'chrom': '3', 'pos':30691759, 'ref':'A', 'alt':'', 'transcript':'NM_003242.5'},
                {'chrom': '5', 'pos':127730815, 'ref':'C', 'alt':'G', 'transcript':'NM_001999.3'},
            ]
        }
        self.denovo = [
            {'chrom': '2', 'pos':162060116, 'ref':'T', 'alt':'G'}, # re-inforces existing 5'
            {'chrom': '2', 'pos': 162041853, 'ref': 'T', 'alt':'G'}, # introduces 5'
        ]
        

    def test_refseqgene(self):
        for record in self.data:
            effect = p.predict(record['chrom'],
                               record['pos'],
                               record['ref'],
                               record['alt'],
                               refseq=REFSEQ,
                               refseqgene=REFSEQGENE)
            for key in self.keys:
                assert key in effect[0][0]
                

    def test_genepanel(self):
        for record in self.genepanel['data']:
            effect = p.predict(record['chrom'],
                               record['pos'],
                               record['ref'],
                               record['alt'],
                               refseq=REFSEQ,
                               genepanel=self.genepanel['genepanel'])
            for key in set(self.genepanel['data'][0].keys()) - set(['chrom', 'pos', 'ref', 'alt']):
                assert effect[0][0][key] == record[key]
                    

    def test_sequences(self):
        '''
        Check that the wild and mutated sequenced are correct.
        '''
        for record in self.sequences:
            effect = p.predict(record['chrom'],
                               record['pos'],
                               record['ref'],
                               record['alt'],
                               refseq=REFSEQ,
                               refseqgene=REFSEQGENE)
            for key in set(self.sequences[0].keys()) - set(['chrom', 'pos', 'ref', 'alt']):
                    assert effect[0][0][key] == record[key]

    def test_denovo(self):
        '''
        Check that the wild and mutated sequenced are correct.
        '''
        for record in self.denovo:
            effect = p.predict(record['chrom'],
                               record['pos'],
                               record['ref'],
                               record['alt'],
                               refseq=REFSEQ,
                               refseqgene=REFSEQGENE)
            assert effect[0][1]['effect_descr'] == 'de_novo'
                    
if __name__ == "__main__":
    unittest.main()