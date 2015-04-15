# -*- coding: utf-8 -*-
"""
SSFL

Implementation of:
Shapiro, M.B., Senapathy, P. "RNA splice junctions of different classes of eukaryotes:
Sequence statistics and functional implications in gene expression."
Nucleic Acids Research, Volume 15, Issue 17, 11 September 1987, Pages 7155-7174

The Position Weight Matrix (PWM) is updated after Zhang98

Notes: the algorithm scans for GT (AG) for candidate donor (acceptor) sites, 
before scoring the 9 (17) nucleotides long candidate sequence. 

TODO use Bio.motifs
Check at http://ibis.tau.ac.il/ssat/SpliceSiteFrame.htm

@author: Hugues Fontenelle, 2014
"""

from splice.splice_base import SpliceBase, reverse_complement

class SSFL( SpliceBase ):
    
    # ------------------------------------------------
    def __init__(self, seq="", offset=0):
        super(SSFL, self).__init__(seq, offset)
        
        # hard-coded thresholds 
        self.donorThreshold = 70 # [0-100]
        self.acceptorThreshold = 70 # [0-100]
               
        # Zhang98
        # Position Weight Matrix for human donor 5' splice site GT
        #                       -3  -2  -1   0   1   2   3   4   5
        #                        *   *   *   G   T   *   *   *   *
        # candidate consensus :  A   A   G   G   T   A   A   G   T 
        # donor AAGGTAAGT
        self.donorPWM = {'A' : [38, 62, 12,  0,  0, 71, 73, 11, 21], 
                         'C' : [31, 10,  4,  0,  0,  2,  6,  6, 10], 
                         'G' : [18, 12, 77,100,  0, 24,  8, 75, 14], 
                         'T' : [13, 16,  7,  0,100,  3, 13,  8, 55]} 
        self._maxt =  38+62+77+100+100+71+73+75+55
        self._mint =  13+10 +4  +0  +0 +2 +6 +6+10  
        # Position Weight Matrix for human acceptor 3' splice site AG
        #                         -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1   0   1  2  3
        #                         <------ pyridimides --------------->  N  c   A   G  g  t   
        # candidate consesus :      T   T   T   T  T  T  T  T  T  T  T  T  C   A   G  G  T
        # acceptor TTTTTTTTTTTTCAGGT
        self.acceptorPWM = {'A' : [15, 14, 13, 11,10,10,11,12,13,11,10,26, 7,100,  0,26,24],
                            'C' : [24, 21, 20, 22,21,22,25,28,28,25,22,25,55,  0,  0,11,15],
                            'G' : [10, 12, 10,  9,10, 9,10,10, 8, 5, 5,15, 1,  0,100,50,20],
                            'T' : [51, 53, 57, 58,59,59,54,50,51,59,63,33,37,  0,  0,13,41]}
        self._l1 = 9+10+9+10 +10+8+5+5
        self._h1 = 53+57+58+59 +59+54+59+63
        self._l2 = 1+0+0+11
        self._h2 = 55+100+100+50
    
    # ------------------------------------------------
    def score_donor_by_position(self, pos):
        # for donors, pos 0 refers to the G of GT right after the splice junction
        # (ie pos 0 is the first nt in the intron)
        #  5' A A G|G T A A G T 3'
        #    -3-2-1|0 1 2 3 4 5
        candidate = self.seq[pos-3:pos+6]
        score = self.score_donor_by_sequence(candidate)
        return {(pos, 'Donor', '+'): score} # candidate.tostring()]

    # ------------------------------------------------
    def score_acceptor_by_position(self, pos):
        # for acceptors, pos 0 refers to the nt right after the splice junction 
        # (ie pos 0 is the first nt in the exon)
        # 5' T T T T T T T T T T T T C A G|G T 3'
        #  -15 . . . -10  . . . -5-4-3-2-1|0 1
        candidate = self.seq[pos-15:pos+2]
        score = self.score_acceptor_by_sequence(candidate)
        return {(pos, 'Acceptor', '+'): score} #candidate.tostring()]
 
    # ------------------------------------------------
    def score_donor_by_sequence(self, candidate):
        if len(candidate) != 9:
            return 0
        score = 0    
        ind = 0
        for c in candidate:
            score += self.donorPWM[c][ind]
            ind += 1
        score = 100*(score- self._mint)/float(self._maxt - self._mint)    
        return score 
        
    # ------------------------------------------------
    def score_acceptor_by_sequence(self, candidate):
        if len(candidate) != 17:
            return 0
  
        #  if AG is found in pyrimidine part then candidate is rejected
        if candidate[0:11].find("AG") != -1: 
            return 0
  
        t1 = []
        ind = 0
        for c in candidate[0:11]:
            t1.append( self.acceptorPWM[c][ind] )
            ind += 1
        t1.sort()
        t1.reverse()
        t1 = sum(t1[0:8])
        t2 = 0    
        ind = 12
        for c in candidate[ind:ind+4]:
            t2 += self.acceptorPWM[c][ind]
            ind += 1
        score1 = 100*((t1-self._l1)/float(self._h1-self._l1))
        score2 = 100*((t2-self._l2)/float(self._h2-self._l2))
        score = 0.5 * (score1 + score2) 
        return score 

    # ------------------------------------------------
    def find_donors(self, reverse=False):
        seq = self.seq
        if reverse==True:
            seq = reverse_complement(seq)
        s = {}
        ind = 0
        while ind >= 0:            
            ind = seq.find("GT", ind+1)
            candidate = seq[ind-3:ind+6]
            score = self.score_donor_by_sequence(candidate)
            if score > self.donorThreshold:
                if reverse==True:
                    s.update( {( len(seq)-ind+self.offset, 'Donor', '-'): score} ) #candidate.reverse_complement().tostring()   
                else:
                    s.update( {( ind+self.offset, 'Donor', '+'): score } ) #candidate.tostring()                 
        return s
    
    # ------------------------------------------------
    def find_acceptors(self, reverse=False):
        seq = self.seq
        if reverse==True:
            seq = reverse_complement(seq)
        s = {}
        ind = 0
        while ind >= 0:            
            ind = seq.find("AG", ind+1)        
            candidate = seq[ind-13:ind+4]
            score = self.score_acceptor_by_sequence(candidate)        
            if score > self.acceptorThreshold:
                if reverse==True:
                    # the index is shifted by 2, so that 0 now correspond to the splice junction itself
                    s.update( {(len(seq)-(ind+2)+self.offset, 'Acceptor', '-'): score} ) #candidate.reverse_complement().tostring()
                else:
                    # the index is shifted by 2, so that 0 now correspond to the splice junction itself
                    s.update( {(ind+2+self.offset, 'Acceptor', '+'): score} ) #candidate.tostring()           
        return s
    
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
        super(SSFL, self).export_to_bed(filename, geneID, sites, "SSFL")


# ============================================================
if  __name__ == "__main__":
    # usage example
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
    r = s1.score_acceptor_by_position(65) # Note: scores only + strand
    print "score= ", r.items()[0][1]     
    print "score= ", s1.score_donor_by_sequence("AAGGTAAGT") # Note: scores only + strand
    #s1.export_to_bed("test.bed","myID",sites) 