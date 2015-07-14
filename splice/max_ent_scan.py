# -*- coding: utf-8 -*-
"""
MaxEntScan

Calls the PERL scripts from Yeo03maxEntScan

echo acggtaagt > t ; perl score5.pl t | cut -f2

Implementation of:
Gene Yeo and Christopher B. Burge "Maximum entropy modeling of short sequence
motifs with applications to RNA splicing signals."
J Comput Biol (2004) vol. 11 (2-3) pp. 377-94

@author: Hugues Fontenelle, 2014
"""

import os.path, inspect
import itertools
import tempfile
import subprocess

from splice.splice_base import SpliceBase, reverse_complement

REL_PATH_TO_3RDPARTY = '../thirdparty/maxentscan/'

class MaxEntScan( SpliceBase ):

    # ------------------------------------------------
    def __init__(self, seq="", offset=0):
        super(MaxEntScan, self).__init__(seq, offset)

        # hard-coded thresholds
        self.donorThreshold = 0.0 # [0-16]
        self.acceptorThreshold = 0.0 # [0-16]

    # ------------------------------------------------
    def score_donor_by_position(self, pos):
        # for donors, pos 0 refers to the G of GT right after the splice junction
        # (ie pos 0 is the first nt in the intron)
        #  5' A A G|G T A A G T 3'
        #    -3-2-1|0 1 2 3 4 5
        candidate = self.seq[pos-3:pos+6]
        score = self.score_donor_by_sequence(candidate)
        return {(pos, 'Donor', '+'): score}

    # ------------------------------------------------
    def score_acceptor_by_position(self, pos):
        # for acceptors, pos 0 refers to the nt right after the splice junction
        # (ie pos 0 is the first nt in the exon)
        # 5' T T T T T T T T T T T T C A G|G T 3'
        #  -15 . . . -10  . . . -5-4-3-2-1|0 1
        candidate = self.seq[pos-20:pos+3]
        score = self.score_acceptor_by_sequence(candidate)
        return {(pos, 'Acceptor', '+'): score}

    # ------------------------------------------------
    def score_donor_by_sequence(self, candidate):
        # Each sequence must be 9 bases long. [3 bases in exon][6 bases in intron]
        s = self.find_splice_sites(seq=candidate, splice_type="Donor", reverse=False)
        return s.values()[0]

    # ------------------------------------------------
    def score_acceptor_by_sequence(self, candidate):
        # Each sequence must be 23 bases long. [20 bases in the intron][3 base in the exon]
        s = self.find_splice_sites(seq=candidate, splice_type="Acceptor", reverse=False)
        return s.values()[0]

    # ------------------------------------------------
    def find_splice_sites(self, seq=None, splice_type="Donor", reverse=False):
        '''
        Find and score splice sites from sequence
        Parameters:
         - seq: FASTA string, possibly containing several candidate sites.
                If empty, uses the class variable self.seq
                Can be as short as one candidate splice site sequence
         - splice_type: 'Donor' (default) or 'Acceptor'
         - reverse: if True, score the reverse strand (default False)
        Output:
         - returns a dictionary of {(pos, splice_type, strand): score}
                pos is w.r.t. the sequence,
                but an offset can be added (from class variable self.offset)
        '''
        if not seq:
            seq=self.seq

        if splice_type in ["Donor", "donor"]:
            pattern = "GT"
            ind_left = 3
            ind_right = 6
            ind_off = 0
            splice_type = "Donor"
            program = "score5.pl"
            threshold = self.donorThreshold
        elif splice_type in ["Acceptor", "acceptor"]:
            pattern = "AG"
            ind_off = 2
            ind_left = 20 - ind_off
            ind_right = 3 + ind_off
            splice_type = "Acceptor"
            program = "score3.pl"
            threshold = self.acceptorThreshold
        else:
            raise "Invalid splice type (", splice_type, ")\nMust be Donor|Acceptor."

        strand = "+"
        if reverse==True:
            seq = reverse_complement(seq)
            strand = "-"
        s = {}
        ind = 0

        path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        path = os.path.join(path, REL_PATH_TO_3RDPARTY)

        candidate_list = []

        outfd, outsock_path = tempfile.mkstemp(suffix='.tmp', prefix='MaxEntScan', dir=None, text=True)
        with open(outsock_path, 'w') as f:
            while ind >= 0:
                ind = seq.find(pattern, ind+1)
                candidate = seq[ind-ind_left:ind+ind_right]
                if len(candidate) == ind_left+ind_right:
                    f.write("> candidate "+str(ind+ind_off)+"\n")
                    f.write(str(candidate))
                    candidate_list.append([ind+ind_off, candidate])

        args = ['perl', program, outsock_path]
        proc = subprocess.Popen(args, cwd=path, stdout=subprocess.PIPE)
        score = proc.stdout.read().split('\n')
        proc.terminate()

        os.close(outfd)
        os.remove(outsock_path)

        score = filter(None, score)
        score = [float(line.split('\t')[-1]) for line in score]

        selectors = [x > threshold  for x in score]
        if reverse==True:
            for key1, key2 in itertools.compress(zip(score, candidate_list), selectors):
                s.update( {(len(seq)-key2[0] + self.offset, splice_type, strand) : key1 } )
        else:
            for key1, key2 in itertools.compress(zip(score, candidate_list), selectors):
                s.update( {(key2[0] + self.offset, splice_type, strand): key1 } )

        return s

    # ------------------------------------------------
    def find_donors(self, reverse=False):
        return self.find_splice_sites(splice_type="Donor", reverse=reverse)

    # ------------------------------------------------
    def find_acceptors(self, reverse=False):
        return self.find_splice_sites(splice_type="Acceptor", reverse=reverse)

    # ------------------------------------------------
    def find_all_sites(self, splice_type='all', strand='all'):
        sites = {}

        if splice_type in ['Donor', 'all', None] and strand in ['+', 'all', None]:
            sites.update( self.find_donors() )
        if splice_type in ['Acceptor', 'all', None] and strand in ['+', 'all', None]:
            sites.update( self.find_acceptors() )
        if splice_type in ['Donor', 'all', None] and strand in ['-', 'all', None]:
            sites.update( self.find_donors(reverse=True) )
        if splice_type in ['Acceptor', 'all', None] and strand in ['-', 'all', None]:
            sites.update( self.find_acceptors(reverse=True) )

        return sites

    # ------------------------------------------------
    def export_to_bed(self, filename, geneID, sites):
        super(MaxEntScan, self).export_to_bed(filename, geneID, sites, "MaxEntScan")


# ============================================================
if  __name__ == "__main__":
    # usage example
    seq1 = "\
GAGGTAAGTggggggggggggggggggggggggggggggggggggggggg\
TTTCTTTTTTTTTTTGGCAGTGTggggggggggggggggggggggggggg\
ACTTACCTCggggggggggggggggggggggggggggggggggggggggg\
ACACTGCCAAAAAAAAAAAGAAAggggggggggggggggggggggggggg".upper()
    # For testing purpose, this sequence is manufactured to show 50 nt per line
    # line 1: Donor site, + strand     GAG|GTAAGT                 0+ 3
    # line 2: Acceptor site, + strand  TTTCTTTTTTTTTTTGGCAG|TGT  50+20
    # line 3: Donor site, - strand     ACTTAC|CTC               100+ 6
    # line 4: Acceptor site, - strand  ACA|CTGCCAAAAAAAAAAAGAAA 150+ 3
    s1 = MaxEntScan()
    s1.seq = seq1
    sites = s1.find_all_sites()
    r = s1.score_acceptor_by_position(70) # Note: scores only + strand
    print "score= ", r.items()[0][1]
    print "score= ", s1.score_donor_by_sequence("AAGGTAAGT") # Note: scores only + strand
    print "score= ", s1.score_acceptor_by_sequence("TTTCTTTTTTTTTTTGGCAGTGT")
    #s1.export_to_bed("test.bed","myID",sites)
