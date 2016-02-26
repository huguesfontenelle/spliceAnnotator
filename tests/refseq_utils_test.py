# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 12:39:58 2014

@author: huguesfo
"""


from unittest import TestCase
from splice import refseq_utils as rf

REFSEQGENE = "/Users/huguesfo/genevar/vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
GENEPANEL = "/Users/huguesfo/genevar/vcpipe-bundle/clinicalGenePanels/Bindevev_v02/Bindevev_v02.transcripts.csv"

class TestRefseqUtils(TestCase):

    def setUp(self):
        self.data = [
            {
            'chrom': '7', 'pos': 5097400, 'start': 5097378, 'end': 5097453,
            'fasta': 'TGGGAGGTGAATTTCCATGTCAACATAGTCCAGGTAAGTTAGTAGCGTATCAAAGGTTAAAAAATGCTCATCCCAG',
            'authentic': {'pos': 5097410, 'splice_type': 'Donor', 'strand': '+'}
            },
            {
             'chrom': '11', 'pos': 33777508, 'start': 33777493, 'end': 33777548,
             'fasta': 'GATTAGACAGTGCCATGCTTCCCAATAACCTGAGGAGCATACAGACTCAAACATGT',
             'authentic': {'pos': 33777521, 'splice_type': 'Acceptor', 'strand': '-'}
            },
            {
             'chrom': '14', 'pos': 58055916, 'start': 58055911, 'end': 58055966,
             'fasta': 'TGGGCAGTTAACATGAAACCTACCTGAATTTTTTCATTGGAGATTGCTTTTCTTGA',
             'authentic': {'pos': 58055933, 'splice_type': 'Donor', 'strand': '-'}
            }]
        self.data_genepanel = [
            {'chrom': '5', 'pos': 127640800, 'authentic': {'pos': 127640774, 'splice_type': 'Acceptor', 'strand': '-'}},
            {'chrom': '5', 'pos': 127599100, 'authentic': {'pos': 127599116, 'splice_type': 'Donor', 'strand': '-'}},
            {'chrom': '3', 'pos': 30713100, 'authentic': {'pos': 30713129, 'splice_type': 'Acceptor', 'strand': '+'}},
            {'chrom': '3', 'pos': 30691900, 'authentic': {'pos': 30691952, 'splice_type': 'Donor', 'strand': '+'}}]
        self.data_auth = [
            {
             'chrom': '8', 'pos': 143958300,
             'authentic': {'pos': 143958301, 'splice_type': 'Acceptor', 'strand': '-'}
            },
            {
             'chrom': '8', 'pos': 143958100,
             'authentic': {'pos': 143958097, 'splice_type': 'Donor', 'strand': '-'}
            },
            {
             'chrom': '8', 'pos': 143922530,
             'authentic': {'pos': 143922533, 'splice_type': 'Acceptor', 'strand': '+'}
            },
            {
             'chrom': '8', 'pos': 143922630,
             'authentic': {'pos': 143922641, 'splice_type': 'Donor', 'strand': '+'}
            }]

    def test_fasta_refseq(self):
        for data in self.data:
            fasta = rf.get_fasta(data['chrom'], data['start'], data['end'], refseq=REFSEQ)
            self.assertEqual(fasta, data['fasta'])

    def test_fasta_NCBI(self):
        for data in self.data:
            fasta = rf.get_fasta(data['chrom'], data['start'], data['end'])
            self.assertEqual(fasta, data['fasta'])

    def test_auth_refseqgene(self):
        for data in self.data:
            auth = rf.get_closest_authentic(data['chrom'], data['pos'], refseqgene=REFSEQGENE)
            self.assertEqual({k:auth[k] for k in ['pos', 'splice_type', 'strand']}, data['authentic'])
        for data in self.data_auth:
            auth = rf.get_closest_authentic(data['chrom'], data['pos'], refseqgene=REFSEQGENE)
            self.assertEqual({k:auth[k] for k in ['pos', 'splice_type', 'strand']}, data['authentic'])

    def test_auth_NCBI(self):
        for data in self.data:
            auth = rf.get_closest_authentic(data['chrom'], data['pos'])
            self.assertEqual({k:auth[k] for k in ['pos', 'splice_type', 'strand']}, data['authentic'])
        for data in self.data_auth:
            auth = rf.get_closest_authentic(data['chrom'], data['pos'])
            self.assertEqual({k:auth[k] for k in ['pos', 'splice_type', 'strand']}, data['authentic'])

    def test_auth_genepanel(self):
        for data in self.data_genepanel:
            auth = rf.get_closest_authentic(data['chrom'], data['pos'], genepanel=GENEPANEL)
            self.assertEqual({k:auth[k] for k in ['pos', 'splice_type', 'strand']}, data['authentic'])

#=========================================================
if  __name__ == "__main__":
    unittest.main()
