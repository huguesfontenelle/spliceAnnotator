"""
Tests for SplicePredict.

@author: Hugues Fontenelle, 2014
"""

from unittest import TestCase, SkipTest
from splice import splice_predict as p

REFSEQGENE = "/Users/huguesfo/genevar/vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab" # RefSeqGene definitions
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

    @SkipTest
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
        for record in self.data:
            effect = p.predict(record['chrom'],
                               record['pos'],
                               record['ref'],
                               record['alt'],
                               refseq=REFSEQ,
                               genepanel=GENEPANEL)
            for key in self.keys:
                if effect and effect[0] and effect[0][0]:
                    assert key in effect[0][0]
