# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 12:39:25 2014

@author: huguesfo
"""

import os
import csv
from pyfaidx import Fasta

email = "hugues.fontenelle@medisin.uio.no"

# Human Genome Assembly Data
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
RefSeqAccession = {'hg19': {
    '1': "NC_000001.10",    '2': "NC_000002.11",    '3': "NC_000003.11",
    '4': "NC_000004.11",    '5': "NC_000005.9",    '6': "NC_000006.11",
    '7': "NC_000007.13",    '8': "NC_000008.10",    '9': "NC_000009.11",
    '10': "NC_000010.10",    '11': "NC_000011.9",    '12': "NC_000012.11",
    '13': "NC_000013.10",    '14': "NC_000014.8",    '15': "NC_000015.9",
    '16': "NC_000016.9",    '17': "NC_000017.10",    '18': "NC_000018.9",
    '19': "NC_000019.9",    '20': "NC_000020.10",    '21': "NC_000021.8",
    '22': "NC_000022.10",    'X': "NC_000023.10",    'Y': "NC_000024.9",
},
'hg38': {
    '1': "NC_000001.11",    '2': "NC_000002.12",    '3': "NC_000003.12",
    '4': "NC_000004.12",    '5': "NC_000005.10",    '6': "NC_000006.12",
    '7': "NC_000007.14",    '8': "NC_000008.11",    '9': "NC_000009.12",
    '10': "NC_000010.11",    '11': "NC_000011.10",    '12': "NC_000012.12",
    '13': "NC_000013.11",    '14': "NC_000014.9",    '15': "NC_000015.10",
    '16': "NC_000016.10",    '17': "NC_000017.11",    '18': "NC_000018.10",
    '19': "NC_000019.10",    '20': "NC_000020.11",    '21': "NC_000021.9",
    '22': "NC_000022.11",    'X': "NC_000023.11",    'Y': "NC_000024.10",
}
}

try:
    import cPickle as pickle
except:
    import pickle

__all__ = ["memoize"]

def memoize(function, limit=None):
    if isinstance(function, int):
        def memoize_wrapper(f):
            return memoize(f, function)

        return memoize_wrapper

    dict = {}
    list = []
    def memoize_wrapper(*args, **kwargs):
        key = pickle.dumps((args, kwargs))
        try:
            list.append(list.pop(list.index(key)))
        except ValueError:
            dict[key] = function(*args, **kwargs)
            list.append(key)
            if limit is not None and len(list) > limit:
                del dict[list.pop(0)]

        return dict[key]

    memoize_wrapper._memoize_dict = dict
    memoize_wrapper._memoize_list = list
    memoize_wrapper._memoize_limit = limit
    memoize_wrapper._memoize_origfunc = function
    memoize_wrapper.func_name = function.func_name
    return memoize_wrapper


# ------------------------------------------------------------
def chr_to_refseq(chrom, ref='hg19'):
    '''Returns the NCBI Reference SequenceID from chromosome number'''
    return RefSeqAccession[ref].get(chrom, "Unkwnown chromosome number")


# ------------------------------------------------------------
def get_fasta(chrom, start, end, refseq=None, ref='hg19'):
    '''Retrieves FASTA sequence'''
    fasta = ''
    chrom = str(chrom)
    if chrom.startswith("chr"): chrom = chrom[3:]
    if refseq:
        '''Retrieves FASTA sequence from local refSeq file'''
        if os.path.isfile(refseq):
            genome = Fasta(refseq)
            fasta = str(genome[chrom][start-1:end])
    else:
        '''Retrieves FASTA sequence from NCBI (internet connection)'''

        from Bio import Entrez, SeqIO

        refseq_id = chr_to_refseq(chrom, ref)
        Entrez.email = email
        handle = Entrez.efetch(db="nucleotide",
                               id=refseq_id,
                               rettype="fasta",
                               strand=1,
                               seq_start=start,
                               seq_stop=end)
        entrez_record = SeqIO.read(handle, "fasta")
        handle.close()
        fasta = str(entrez_record.seq)

    return fasta


# ------------------------------------------------------------
def get_auth_from_refseqgene(chrom, pos, refseqgene=None,  ref='hg19'):
    '''Retrieves exon definitions from local refSeqGene file'''
    refSeq = []
    try:
        with open(refseqgene, 'r') as f:
            refSeqList = csv.reader(f, delimiter='\t')
            for row in refSeqList:
                if row[0].strip().startswith('#'):
                    continue
                if row[2][3:] == chrom and int(row[4]) <= pos <= int(row[5]):
                    refSeq = row
                    break
    except Exception as ex:
        raise ex

    if refSeq == []:
        print 'chr%s:%s Not in a RefSeq transcript\n' % (chrom, str(pos))
        return None

    exonList = [(int(x), int(y)) for (x,y) in zip(refSeq[9].rstrip(',').split(','), refSeq[10].rstrip(',').split(','))]
    strand = refSeq[3]
    transcript = refSeq[1]

    rel_exons = [(abs(start-pos), abs(end-pos)) for (start, end) in exonList]
    a = [(min(e), e.index(min(e))) for i, e in enumerate(rel_exons)]
    id0 = a.index(min(a))
    id1 = min(a)[1]
    closest_splice_site = exonList[id0][id1]

    if strand == '+':
        splice_type = {0: "Acceptor", 1: "Donor"}[id1]
    else:  # the convention is opposite for - strand!!
        splice_type = {1: "Acceptor", 0: "Donor"}[id1]

    return {'chrom': chrom,
            'pos': closest_splice_site, 'splice_type': splice_type,
            'strand': strand, 'transcript': transcript}


# ------------------------------------------------------------
def get_auth_from_genepanel(chrom, pos, genepanel=None, ref='hg19'):
    '''Retrieves exon definitions from local genepanel file'''
    gp = []
    try:
        with open(genepanel, 'r') as f:
            gpList = csv.reader(f, delimiter='\t')
            for row in gpList:
                if row[0].strip().startswith('#'):
                    continue
                txStart, txEnd = int(row[1]), int(row[2])
                if row[0] == chrom and txStart <= pos <= txEnd:
                    gp = row
                    break
    except Exception as ex:
        raise ex

    if gp == []:
        print 'chr%s:%s Not in a genepanel transcript\n' % (chrom, str(pos))
        return None

    exonList = [(int(x), int(y)) for (x,y) in zip(gp[12].rstrip(',').split(','), gp[13].rstrip(',').split(','))]
    strand = gp[5]
    transcript = gp[3]

    rel_exons = [(abs(start-pos), abs(end-pos)) for (start, end) in exonList]
    a = [(min(e), e.index(min(e))) for i, e in enumerate(rel_exons)]
    id0 = a.index(min(a))
    id1 = min(a)[1]
    closest_splice_site = exonList[id0][id1]

    if strand == '+':
        splice_type = {0: "Acceptor", 1: "Donor"}[id1]
    else:  # the convention is opposite for - strand!!
        splice_type = {1: "Acceptor", 0: "Donor"}[id1]

    return {'chrom': chrom, 'pos':closest_splice_site, 'splice_type': splice_type,
            'strand': strand, 'transcript': transcript}


# ------------------------------------------------------------
def get_auth_from_NCBI(chrom, pos, ref='hg19'):
    '''Retrieves exon definitions from from NCBI (internet connection)'''
    from Bio import Entrez, SeqIO

    refseq_id = chr_to_refseq(chrom, ref)
    startpos = pos - 5000
    endpos = pos + 5000
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide",
                           id=refseq_id,
                           rettype="gb",
                           retmode="text",
                           seq_start=startpos,
                           seq_stop=endpos)
    entrez_record = SeqIO.read(handle, "gb")
    handle.close()

    # find mRNA exons start:end in the subsequence
    # exons = {}
    exonList = []
    for feature in entrez_record.features:
        if feature.type == 'mRNA':
            if not feature.sub_features:
                exonList.append((feature.location.start.position + startpos,
                                 feature.location.end.position + startpos,
                                 feature.location.strand))
            else:
                for sub_feature in feature.sub_features:
                    exonList.append((sub_feature.location.start.position + startpos,
                                     sub_feature.location.end.position + startpos,
                                     sub_feature.location.strand))
    if not exonList:
        print 'chr%s:%s No authentic splice site nearby\n' % (chrom, str(pos))
        return None

    rel_exons = [(abs(start-pos), abs(end-pos)) for (start, end, s) in exonList]
    a = [(min(e), e.index(min(e))) for i, e in enumerate(rel_exons)]
    id0 = a.index(min(a))
    id1 = min(a)[1]
    closest_splice_site = exonList[id0][id1]

    strand = {1: "+", -1: "-"}[exonList[id0][2]]
    if strand == '+':
        splice_type = {0: "Acceptor", 1: "Donor"}[id1]
    else:  # the convention is opposite for - strand!!
        splice_type = {1: "Acceptor", 0: "Donor"}[id1]

    return {'chrom': chrom, 'pos': closest_splice_site-1,
            'splice_type': splice_type,
            'strand': strand, 'transcript': refseq_id}


# ------------------------------------------------------------
@memoize(10)
def get_closest_authentic(chrom, pos, refseqgene=None, genepanel=None, ref='hg19', get_sequence=False, refseq=None, seq_size=50):
    '''
    Find the closest intron-exon junction in the RefSeq w.r.t. the mutation
    position. Returns a dictionnary with the splice site position, splice type,
    strand for a transcript.
    Eg:
    {'pos': 32944538, 'splice_type': 'Acceptor', 'strand': '+', 'transcript': 'NM_000059'}
    '''
    if chrom.startswith("chr"): chrom = chrom[3:]
    if refseqgene:
        auth = get_auth_from_refseqgene(chrom, pos, refseqgene=refseqgene, ref='hg19')
    elif genepanel:
        auth = get_auth_from_genepanel(chrom, pos, genepanel=genepanel, ref='hg19')
    else:
        auth = get_auth_from_NCBI(chrom, pos, ref='hg19')

    if auth and get_sequence:
        auth['fasta'] = get_fasta(chrom=auth['chrom'], start=auth['pos']-seq_size, end=auth['pos']+seq_size, refseq=refseq)

    return auth


# ============================================================
if __name__ == "__main__":
    pass
