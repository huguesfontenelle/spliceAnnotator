# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 12:13:16 2014

@author: Hugues Fontenelle, 2014
"""

# ============================================================
def reverse_complement(seq):
    """Return reverse complement of DNA strand"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq_rev_cpl = "".join(complement.get(base, base) for base in reversed(seq))
    return seq_rev_cpl

# ============================================================
# base class SpliceBase
# all splice site detection algorithms derive from this class, in separate files
# ============================================================
class SpliceBase( object ):
    
    def __init__(self, seq="", offset=0):
        self.seq = seq
        self.offset = offset
        
    # ------------------------------------------------
    def export_to_bed(self, filename, geneID, sites, trackname):
        # color = {'+':"0,0,255", '-':"0,255,0"}
        color = {'Donor':"0,0,255", 'Acceptor':"0,255,0"}
        # TODO color by strand or by splice_type ?
        bed_file = open(filename,"w")
        bed_file.write("track name="+trackname+" description=\"Splice Sites\"\n")
        for site,score in sites.iteritems():
            (pos, splice_type, strand) = site
            bedline = [geneID, str(pos), str(pos+1), splice_type, "{:.2f}".format(score),\
            strand, str(pos), str(pos+1), color[splice_type], "\n" ]
            bed_file.write(" ".join(bedline))
        bed_file.close()
        return 0

    # ------------------------------------------------   
    def filter_sites(self, sites, *args, **kwargs):
        '''
        Filter sites by splice_type, strand and/or greater than threshold
        
        Usage examples:
        $ filtered = filter_sites(sites, splice_type='Acceptor', strand='-', threshold=60.0)
        $ filtered = filter_sites(sites, splice_type='Donor')
        $ filtered = filter_sites(sites, strand='+')
        $ filtered = filter_sites(sites, threshold=60.0)
        
        Notes:
        The 'sites' data structure is a dictionary where the keys are tuples of 
            (position, splice_type, strand), and the value is the score.
        $ sites = {(3, 'Donor', '+'): 100.0,
                 (65, 'Acceptor', '+'): 100.0}

        '''        
        strand = kwargs.get('strand', None)
        splice_type = kwargs.get('splice_type', None)
        threshold = kwargs.get('threshold', None)
        
        if strand not in ['+','-', 'all', None]:
            raise ValueError("%r is not a proper strand orientation. Specify '+' or '-'." % strand)
        if splice_type not in ['Donor','Acceptor', 'all', None]:
            raise ValueError("%r is not a proper splice type. Specify 'Donor' or 'Acceptor'." % strand)
        
        filtered = sites
        
        if splice_type not in ['all', None]:                   
            filtered = {k: v for k, v in filtered.iteritems() if k[1] == splice_type }
        if strand not in ['all', None]:                   
            filtered = {k: v for k, v in filtered.iteritems() if k[2] == strand }
        if threshold is not None:        
            filtered = {k: v for k, v in filtered.iteritems() if v >= threshold }  
    
        return filtered

    # ------------------------------------------------   
    def count_decoy(self, splice_type='Donor'):
        '''
        Count potential Donor (GT) -or Acceptor (AG)- sites on both strands.
        Does not substract the actual authentic sites (unknown at this point)
        '''
        if splice_type == 'Donor':
            pattern = 'GT'
        elif splice_type == 'Acceptor':
            pattern = 'AG' 
        
        rseq = reverse_complement(self.seq)
        
        counter = -1;
        ind = 0
        while ind >= 0:            
            ind = self.seq.find(pattern, ind+1)
            counter += 1
        ind = 0
        while ind >= 0:            
            ind = rseq.find(pattern, ind+1)
            counter += 1    
            
        return counter-1 
            
            
            
            
            
            
