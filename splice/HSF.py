# -*- coding: utf-8 -*-
'''
Human Splicing Finder (HSF)

Implementation of:
Desmet, F.-O., Hamroun, D., Lalande, M., Collod-Béroud, G., Claustres, M.,
& Béroud, C. (2009). Human Splicing Finder: an online bioinformatics tool to
predict splicing signals. Nucleic Acids Research, 37(9), e67.
doi:10.1093/nar/gkp215

The web version is available at:
http://www.umd.be/HSF/

@author: Hugues Fontenelle, 2014
'''


from splice.splice_base import SpliceBase, reverse_complement


from bs4 import BeautifulSoup
import urllib


class HSF(SpliceBase):
    '''
    Human Splicing Finder (HSF)

    Splice Site prediction algorithm.
    It is not implemented but instead makes use of the web version at:
    http://www.umd.be/HSF/
    '''

    # ------------------------------------------------
    def __init__(self, seq="", offset=0):
        SpliceBase.__init__(self, seq, offset)

    # ------------------------------------------------
    def submit(self, reverse=False):
        '''
        Using BeautifulSoup4, POST the parameters to webpage, get the results,
        and process to store the HSF splice sites.
        Convert coordinates.
        '''
        if reverse:
            seq = reverse_complement(self.seq)
            strand = '-'
        else:
            seq = self.seq
            strand = '+'

        url = "http://www.umd.be/HSF/4DACTION/input_SSF"
        post_params = {
            'choix_analyse' : 'ssf_sequence',
            'autoselect' : 'yes',
            'snp_select' : 'no',
            'matrice_3' : 'no',
            'showonly' : 'yes',
            'nuclposition' : '',
            'nuclposition5' : '',
            'nuclposition3' : '',
            'choix_bdd' : 'paste',
            'texte' : '',
            'exon_number' : '',
            'intron_number' : '',
            'remLentext' : '4800',
            'sequence' : seq,
            'remLentextmut' : '5000',
            'sequencemut' : '',
            'Firstnucleotide' : '',
            'Lastnucleotide' : '',
            'batch' : '',
            'remLentextBP' : '100',
            'sequenceBP' : '',
            'MDE_sequences' : '',
            'paramfulltables' : 'onlyvariants',
            'fenetreintron' : 'yes',
            'fenetretaille' : '24',
            'paramimages' : 'no',
            'Matrice' : 'PSS',
            'seuil_maxent5' : '0',
            'seuil_maxent3' : '0',
            'seuil_sf2' : '72.9',
            'seuil_sf2_esef' : '1.956',
            'seuil_sf2ig' : '70.51',
            'seuil_sf2ig_esef' : '1.867',
            'seuil_sc35' : '75.05',
            'seuil_sc35_esef' : '2.383',
            'seuil_srp40' : '78.08',
            'seuil_srp40_esef' : '2.67',
            'seuil_srp55' : '73.86',
            'seuil_srp55_esef' : '2.676',
            'seuil_9g8' : '59.245',
            'seuil_tra2' : '75.964',
            'seuil_sironi1' : '60',
            'seuil_sironi2' : '60',
            'seuil_sironi3' : '60',
            'seuil_hnrnpa1' : '65.476'
        }

        post_args = urllib.urlencode(post_params)

        fp = urllib.urlopen(url, post_args)
        soup = BeautifulSoup(fp)

        tr = soup.findAll('tr', {'class':'voir'})
        tr.extend(soup.findAll('tr', {'class':'vu'}))
        record_sites = {}

        for tr1 in tr:
            td = tr1.findAll('td')
            pos = int(td[0].string) #Sequence Position
            #td[1].string #cDNA Position
            splice_type = td[2].string.encode('utf-8') #Splice site type
            #td[3].string #Motif
            #td[4].string #New potential splice site
            score = float(td[5].string) #Consensus value (0-100)
            record_sites.update({(pos, splice_type, strand) : score})

        # reformat record_sites[] to sites{}
        nn_offset = {('Donor', '+') : 2,
                     ('Donor', '-') : -2,
                     ('Acceptor', '+') : 11,
                     ('Acceptor', '-') : -11}
        sites = {}
        if reverse:
            for k, v in record_sites.iteritems():
                sites.update({(len(seq) - k[0] + nn_offset[(k[1], k[2])] + \
                    self.offset, k[1], k[2]): v})
        else:
            for k, v in record_sites.iteritems():
                sites.update({(k[0] + nn_offset[(k[1], k[2])] + \
                    self.offset, k[1], k[2]): v})


        return sites

     # ------------------------------------------------
    def find_all_sites(self, splice_type='all', strand='all'):
        '''
        Maim method.
        Calls submit twice, once for forward, once for reverse strand.
        '''
        sites = {}

        if len(self.seq) > 100000:
            raise RuntimeError("The length of the sequence is too long. \
HSF does not handle sequences longer than 100,000 bp.")

        if strand in ['+', 'all', None]:
            sites.update(self.submit())
        if strand in ['-', 'all', None]:
            sites.update(self.submit(reverse=True))

        if splice_type in ['Donor', 'Acceptor', 'all', None]:
            sites = self.filter_sites(sites, splice_type=splice_type)
            
        return sites

    # ------------------------------------------------
    def export_to_bed(self, filename, geneID, sites):
        '''
        Export function
        '''
        SpliceBase.export_to_bed(self, filename, geneID, sites, "HSF")

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
    s1 = HSF()
    s1.seq = seq1
    sites = s1.find_all_sites()
