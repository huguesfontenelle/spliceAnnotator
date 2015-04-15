# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
NNSplice

Implementation of:
Reese MG, Eeckman, FH, Kulp, D, Haussler, D, 1997.
"Improved Splice Site Detection in Genie". J Comp Biol 4(3), 311-23.

The web version is available at:
http://www.fruitfly.org/seq_tools/splice.html

@author: Hugues Fontenelle, 2014
"""

from splice.splice_base import SpliceBase

from bs4 import BeautifulSoup
import urllib
import re

class NNSplice( SpliceBase ):
    # ------------------------------------------------
    def __init__(self, seq="", offset=0):
        super(NNSplice, self).__init__(seq, offset)
        # hard-coded thresholds 
        self.donorThreshold = 0.0 # [0-1.0]
        self.acceptorThreshold = 0.0 # [0-1.0]    

    # ------------------------------------------------
    def submit(self):    
        url = "http://www.fruitfly.org/cgi-bin/seq_tools/splice.pl"
       
        post_params = {
            'organism' : 'human',     # ['human', 'drosophila']
            'which' : 'both',         # ['acceptors', 'donors', 'both']
            'reverse' : 'yes',        # ['yes', 'no']
            'min_donor' : str(self.donorThreshold),
            'min_acc' : str(self.acceptorThreshold),
            'text' : str(self.seq)
        }
        post_args = urllib.urlencode(post_params)

        fp = urllib.urlopen(url, post_args)
        soup = BeautifulSoup(fp)
        
        prog = re.compile('\d+ +\d+ +[-+]?[0-9]*\.?[0-9]+')
        
        soup.find_all('h3')[2].find_next('pre')
        
        record_sites = [] 
        for h3 in soup.find_all('h3'):
            if 'Donor' in str(h3):
                splice_type = 'Donor'
                strand = '+'
            elif 'Acceptor' in str(h3):
                splice_type = 'Acceptor'
                strand = '+'   
            elif 'reverse strand' in str(h3):
                strand = '-'
            pre = h3.find_next('pre')
            for br in pre.find_all('br'):
                ss = br.next
                
                match = re.search(prog,str(ss))
                if match:
                    start, end, score = match.group().split()
                    record_sites.append( [start, end, score, splice_type, strand] )

        # reformat record_sites[] to sites{}
        nn_offset = {('Donor','+') : 6,
                     ('Donor','-') : -7,
                     ('Acceptor','+') : 20,
                     ('Acceptor','-') : -21}
        sites = {}
        for rec in record_sites:
            start, end, score, splice_type, strand = rec
            pos = int(start) + nn_offset[(splice_type, strand)]
            sites[(pos + self.offset, splice_type, strand)] = float(score)
            
        return sites

     # ------------------------------------------------
    def find_all_sites(self, splice_type='all', strand='all'):
        sites = {} 
        
        if len(self.seq) > 100000:
            raise RuntimeError("The length of the sequence is too long. NNSplice does not handle sequences longer than 100,000 bp.") 

        sites = self.submit()
        
        if splice_type in ['Donor', 'Acceptor']:
            sites = self.filter_sites(sites, splice_type=splice_type)
        if strand in ['+', '-']:
            sites = self.filter_sites(sites, strand=strand)
    
        return sites        
        
    # ------------------------------------------------
    def export_to_bed(self, filename, geneID, sites):
        super(NNSplice, self).export_to_bed(filename, geneID, sites, "NNSplice")       

# ============================================================
if  __name__ == "__main__":
    # usage example
    seq1 = "\
gggggggggggggggggggggggggggggggggggggggggggggggggg\
AAGGTAAGTggggggggggggggggggggggggggggggggggggggggg\
TTTTTTTTTTTTCAGGTggggggggggggggggggggggggggggggggg\
ACTTACCTTggggggggggggggggggggggggggggggggggggggggg\
ACCTGAAAAAAAAAAAAggggggggggggggggggggggggggggggggg".upper()
    # For testing purpose, this sequence is manufactured to show 50 nt per line
    # line 1: Donor site, + strand     aag|GTaagt          0+ 3
    # line 2: Acceptor site, + strand  tttttttttttcAG|gt  50+15
    # line 3: Donor site, - strand     acttAC|ctt        100+ 6
    # line 4: Acceptor site, - strand  ac|CTgaaaaaaaaaaa 150+ 2
    s1 = NNSplice()
    s1.seq = seq1
    sites = s1.find_all_sites()
    #s1.export_to_bed("test.bed","myID",sites)         