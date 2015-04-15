# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 14:06:36 2014

@author: huguesfo
"""

from Bio import Entrez
from Bio import SeqIO
from annotation.splice.SSFL import SSFL
from annotation.splice.max_ent_scan import MaxEntScan
from annotation.splice.gene_splicer import GeneSplicer
from annotation.splice.NNSplice import NNSplice

import pylab as pl
from sklearn.metrics import auc
import numpy as np
import os.path, inspect

genepanel_name = "trusight_tumor.genepanel"
path_to_data = "/Users/huguesfo/Documents/DATA/trusight_tumor/"
Entrez.email = "hugues.fontenelle@medisin.uio.no"

algos = {'SSFL':SSFL, 'MaxEntScan':MaxEntScan, 'GeneSplicer':GeneSplicer, 'NNSplice':NNSplice}

#%%
# ----------------------------------------
# analyze_splice_output: do the stats
# inputs: auth_sites: authentic splice sites in the RefSeqGene
#         crypt_sites : cryptic splice sites found by an algorithm
# outputs: ?
# ----------------------------------------
def analyze_splice_output(decoys, auth_sites, crypt_sites, thresholds):
    '''
    Does some stats.
    '''
    #FPR = []
    #TPR = []
    P = len(auth_sites)
    N = decoys - P   
    sa = set(auth_sites)
    TP = []
    FP = []
    for threshold in thresholds:
        crypt_sites_t = {k: v for k, v in crypt_sites.iteritems() if v >= threshold }
    
        sc = set(crypt_sites_t)       
        # http://en.wikipedia.org/wiki/Confusion_matrix
        TP.append( len(sc & sa) ) # true positives
        FP.append( len(sc - sa) ) # false positives
        #FN = len(sa - sc) # false negatives
        #TN = N- FP        # true negatives    

    return (TP, FP, P, N)

#%%
# ----------------------------------------
# Analyze gene
# ----------------------------------------
def analyze_gene(algo_name, splice_type, gene_name):
    '''
    Inputs:
        algo_name is in ['SSFL', 'MaxEntScan', 'GeneSplicer' ]
        splice_type is in ['Donor', 'Acceptor']
        gene_name is, for example, 'AKT1' or 'BRAF'
    Outputs:
        auth_sites: authentic splice sites, defined as intron-exon junctions in RefSeqGene        
        crypt_sites: cryptic splice sites found by the prediciton algorithm
        decoys: number of splice sites candidates, ie counts any GT (or AG)
    '''
    output_gb = gene_name+ ".gb"
    h_gb = open(path_to_data+output_gb, "rU")
    records = SeqIO.parse(h_gb, "gb")
    record=records.next()
    h_gb.close()
       
    # find mRNA exons start:end in the subsequence; open the GenBank file
    exons = []
    for f in record.features:
        if f.type == 'mRNA':
            for sf in f.sub_features:
                exons.append((sf.location.start.position, # start of exon (acceptor)
                              sf.location.end.position,   # end of exon (donor)
                              sf.location.strand) ) # append to list of tuples          
    
    # open the FASTA as well
    seq = None
    output_fasta = gene_name+ ".fasta"
    h_fasta = open(path_to_data+output_fasta, "rU")
    records = SeqIO.parse(h_fasta, "fasta")
    record=records.next()
    seq = record.seq
    h_fasta.close()
    
    # convert exons to splice sites
    auth_sites = {}
    strand = {1:'+',-1:'-'}
    for exon in exons:
        ss = (exon[0], 'Acceptor',strand[exon[2]])
        auth_sites[ss] = 100
        ss = (exon[1], 'Donor',strand[exon[2]])
        auth_sites[ss] = 100
    
    # run prediction algorithms
    pred = algos[algo_name]()
    pred.seq = seq
    crypt_sites = pred.find_all_sites()
    decoys = pred.count_decoy(splice_type = splice_type)
    
    # filtering 
    auth_sites = pred.filter_sites(auth_sites, splice_type = splice_type)
    crypt_sites = pred.filter_sites(crypt_sites, splice_type = splice_type)
    
    return auth_sites, crypt_sites, decoys

#%%
# ----------------------------------------
# ROC analysis for each algorithm and each splice_type
# ----------------------------------------
def roc_analysis(algo_name, splice_type):
    '''
    Read the gene panel, and iterate on each gene to run the analysis   
    '''
    # read the genepanel file
    genepanel = []
    path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    os.chdir(path)
    h_panel = open(path+'/'+genepanel_name, "rU")  
    for line in h_panel:
        items = line.split()
        if len(items)>1:
            genepanel.append(items)
    h_panel.close()
    
    if algo_name == 'SSFL':
        thresholds = np.arange(50, 100, 5)
    elif algo_name == 'MaxEntScan':
        thresholds = np.arange(0, 16, 1)
    elif algo_name == 'GeneSplicer':
        thresholds = np.arange(0, 16, 1)
    elif algo_name == 'NNSplice':
        thresholds = np.arange(0, 1, 0.1)
        
    T_tp = np.zeros(np.size(thresholds))
    T_fp = np.zeros(np.size(thresholds))
    T_p = 0
    T_n = 0
    for gene in genepanel:
        gene_name = gene[0]
        
        print "Analyzing " + gene_name + " for '" + splice_type + "' with '" + algo_name + "':"
        try:
            auth_sites, crypt_sites, decoys = analyze_gene(algo_name, splice_type, gene_name)
        except RuntimeError: # Skip this one
            print "Skipped!"
            continue
        print "Found " + str(len(auth_sites)) + " authentic sites and " + str(len(crypt_sites)) + " cryptic sites."
        # do the stats
        (TP, FP, P, N) = analyze_splice_output(decoys, auth_sites, crypt_sites, thresholds)
        print "TP_max=" + str(max(TP)) + ", FP_max=" + str(max(FP)) + ", P=" + str(P) + ", N=" + str(N) + "."
        T_tp += TP
        T_fp += FP
        T_p += P
        T_n += N
        
        
    TPR = T_tp / float(T_p)  # true positive rate = sensitivity
    FPR = T_fp / float(T_n) # false positive rate = fall-out
   
    return (FPR, TPR)

#%%
# ============================================================
# Main
# ============================================================
if  __name__ == "__main__":
    
    splice_type = 'Acceptor'
    (FPR1, TPR1) = roc_analysis('SSFL', splice_type)
    (FPR2, TPR2) = roc_analysis('MaxEntScan', splice_type)
    (FPR3, TPR3) = roc_analysis('GeneSplicer', splice_type)
    (FPR4, TPR4) = roc_analysis('NNSplice', splice_type)
    
    # ROC   
    roc_auc1 = auc(FPR1, TPR1)
    roc_auc2 = auc(FPR2, TPR2)
    roc_auc3 = auc(FPR3, TPR3)
    roc_auc4 = auc(FPR3, TPR3)
    print "Area under the ROC curve : \nSSFL %f\nMaxEntScan %f\nGeneSplicer %f\nNNSplice %f" % (roc_auc1, roc_auc2, roc_auc3, roc_auc4)
    
    # Plot ROC curve
    pl.clf()
    pl.plot(FPR1, TPR1, '-b', label='SSFL' )
    pl.plot(FPR2, TPR2, '--c', label='MaxEntScan' )
    pl.plot(FPR3, TPR3, ':m', label='GeneSplicer' )
    pl.plot(FPR4, TPR4, '.-r', label='NNSplice' )
    pl.plot([0, 1], [0, 1], 'k--')
    pl.xlim([0.0, 0.2])
    pl.ylim([0.0, 1.0])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title( splice_type )
    pl.legend(loc="lower right")
    pl.show()
