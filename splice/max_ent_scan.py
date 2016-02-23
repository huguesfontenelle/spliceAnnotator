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

import os.path
import inspect
import tempfile
import subprocess


REL_PATH_TO_3RDPARTY = '../thirdparty/maxentscan/'

def reverse_complement(seq):
    '''
    Return reverse complement of DNA strand
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq_rev_cpl = ''.join(complement.get(base, base) for base in reversed(seq))
    return seq_rev_cpl


def score(seq, strand='+'):
    '''
    Score each each candidate sequence with MaxEntScan
    '''
    plus_strand = ['+', 'forward', '1', 1]
    minus_strand = ['-', 'reverse', '-1', -1]
    
    if type(seq) is str:
        seq = [seq]
        
    assert strand in plus_strand or minus_strand
    for s in seq:
        assert len(s) in [9, 23]
        assert set(s.upper()) <= set('ATCG')
    
    # guess splice type from sequence length
    if len(seq[0]) == 9:
        program = "score5.pl"
    elif len(seq[0]) == 23:
        program = "score3.pl"
    
    if strand in minus_strand:
        seq = [reverse_complement(s) for s in seq]
    
    path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    path = os.path.join(path, REL_PATH_TO_3RDPARTY)

    outfd, outsock_path = tempfile.mkstemp(suffix='.tmp', prefix='MaxEntScan', dir=None, text=True)
    with open(outsock_path, 'w') as f:
        for s in seq:
                f.write("> candidate " + s + "\n")
                f.write(s)

    args = ['perl', program, outsock_path]
    proc = subprocess.Popen(args, cwd=path, stdout=subprocess.PIPE)
    score = proc.stdout.read().split('\n')
    proc.terminate()

    os.close(outfd)
    os.remove(outsock_path)

    score = filter(None, score)
    score = [float(line.split('\t')[-1]) for line in score]
        
    return score