'''
SplicePredict

@author: Hugues Fontenelle, 2014
'''


import os.path
import re
from splice import max_ent_scan as mes
from splice import refseq_utils as rf

THRESHOLD_LOST = -0.216  # threshold at which a drop in MES score is considered "LOST"
SEQ_HALF_SIZE = 50  # half-length of the FASTA sequence retrieved around the splice site
EFFECT_KEYWORD = {
    'LOST': 'predicted_lost',                  # the splice site is predicted to be lost
    'CONSERVED': 'predicted_conserved',        # the splice site is predicted to be conserved
    'NO_EFFECT': 'consensus_not_affected',     # the variant does not impact the splice site (therefore it is conserved)
    'DE_NOVO': 'de_novo',                      # the variant creates or re-inforces a 'de novo' cryptic splice site
    'NOT_TRANSCRIBED': 'not_transcribed',      # the variant is not in the unspliced transcript (pre-mRNA) region
    'NOT_AVAILABLE': 'NA',                     # the prediction failed and is therefore not available
    'SPLICE_DONOR_VARIANT': 'splice_donor_variant',  # (VEP-style annotation of) a splice variant that changes the 2 base pair region at the 5' end of an intron.
    'SPLICE_ACCEPTOR_VARIANT': 'splice_acceptor_variant',  # (VEP-style annotation of) a splice variant that changes the 2 base region at the 3' end of an intron.
}


# ------------------------------------------------
def predict(chrom, pos, ref, alt, refseq=None, refseqgene=None, genepanel=None, denovo=True):
    '''
    Predicts
    '''
    assert os.path.isfile(refseq)
    if refseqgene:
        assert os.path.isfile(refseqgene)
    if genepanel:
        assert os.path.isfile(genepanel)

    if type(alt) is str:
        alt = [alt]

    effects = list()
    for alt1 in alt:
        effects.append(predict_one(chrom, pos, ref, str(alt1), refseq=refseq, refseqgene=refseqgene, genepanel=genepanel, denovo=denovo))

    return effects


# ------------------------------------------------
def predict_one(chrom, pos, ref, alt, refseq=None, refseqgene=None, genepanel=None, denovo=True):
    '''
    Predicts
    '''
    effects = predict_lost_auth(chrom, pos, ref, alt, refseq=refseq,
                                refseqgene=refseqgene, genepanel=genepanel)
    if denovo:
        effects += predict_de_novo(chrom, pos, ref, alt, refseq=refseq,
                                   refseqgene=refseqgene, genepanel=genepanel)
    effects += annotate_2bp_conserved(chrom, pos, ref, alt,
                                      refseq=refseq,
                                      refseqgene=refseqgene,
                                      genepanel=genepanel)

    return effects


# ------------------------------------------------
def predict_de_novo(chrom, pos, ref, alt, refseq=None, refseqgene=None, genepanel=None):

    auth = rf.get_closest_authentic(chrom=chrom, pos=pos, refseqgene=refseqgene,
                                    genepanel=genepanel, refseq=refseq,
                                    get_sequence=True, seq_size=SEQ_HALF_SIZE)
    if not auth:
        return []

    dist = pos - auth['pos']

    consensus_size = {
        ('Donor', '+'): [-2, 7],
        ('Acceptor', '+'): [-19, 4],
        ('Donor', '-'): [-3, 6],
        ('Acceptor', '-'): [-20, 3],
    }
    dimer = {'Donor': 'GT', 'Acceptor': 'AG'}
    s, e = consensus_size[(auth['splice_type'], auth['strand'])]

    fasta_atpos = rf.get_fasta(chrom=chrom, start=pos-SEQ_HALF_SIZE, end=pos+SEQ_HALF_SIZE, refseq=refseq)
    if auth['strand'] == '+':
        fasta_atpos_mut = fasta_atpos[:SEQ_HALF_SIZE] + alt + fasta_atpos[SEQ_HALF_SIZE+len(ref):]
    elif auth['strand'] == '-':
        fasta_atpos = mes.reverse_complement(fasta_atpos)
        fasta_atpos_mut = fasta_atpos[:SEQ_HALF_SIZE] + mes.reverse_complement(alt) + fasta_atpos[SEQ_HALF_SIZE+len(ref):]

    wild = auth['fasta'][SEQ_HALF_SIZE+s:SEQ_HALF_SIZE+e]
    try:
        wild_score = mes.score(wild)
    except Exception:
        return [{'effect_descr': EFFECT_KEYWORD['NOT_AVAILABLE']}]

    ff = lambda x: -s+1 < x < 2*SEQ_HALF_SIZE-e-1
    idx_mut = [m.start()-1 for m in re.finditer(dimer[auth['splice_type']], fasta_atpos_mut)]
    idx_mut = filter(ff, idx_mut)
    idx_wild = [idx if idx<=SEQ_HALF_SIZE else idx+len(ref)-len(alt) for idx in idx_mut]

    dn_wild = [fasta_atpos[idx+s:idx+e] for idx in idx_wild]
    dn_mut = [fasta_atpos_mut[idx+s:idx+e] for idx in idx_mut]

    denovo_seqs = [list(dn) for dn in zip(dn_wild, dn_mut, idx_wild)
                   if dn[0] != dn[1]
                   and len(dn[0]) == len(dn[1]) in [9, 23]]
    if not denovo_seqs:
        return []
    try:
        denovo_scores_w = mes.score([dm[0] for dm in denovo_seqs])
        denovo_scores_m = mes.score([dm[1] for dm in denovo_seqs])
    except AssertionError:
        print 'AssertionError when predicting de novo splice site'
        print '%s:%i%s>%s (wild "%s" --> mut "%s")' % (chrom, pos, ref, alt, dm[0], dm[1])
        return [{'effect_descr': EFFECT_KEYWORD['NOT_AVAILABLE']}]

    effects = []
    for seq, score_w, score_m in zip(denovo_seqs, denovo_scores_w, denovo_scores_m):
        effect = {'effect_descr': EFFECT_KEYWORD['DE_NOVO'],
                  'distance': dist,  # distance btw clostest auth and de novo ss
                  'de_novo_pos': seq[2]-SEQ_HALF_SIZE+pos,  # pos of de novo ss
                  'auth_score': wild_score[0],
                  'wild_score': score_w,
                  'mut_score': score_m,
                  'wild_seq': seq[0],
                  'mut_seq': seq[1],
                  'auth_seq': wild,
                  'auth_pos': auth['pos'],
                  'splice_type': auth['splice_type'],
                  'strand': auth['strand'],
                  'transcript': auth['transcript']}
        effect['de_novo_dist'] = effect['de_novo_pos'] - effect['auth_pos']
        effects.append(effect)

    ff_denovo = lambda x: ((x['wild_score'] <= 0) and (x['mut_score'] >= 4)) or \
                          ((x['wild_score'] > 0) and (x['mut_score'] >= 0) and
                           (x['mut_score'] / x['wild_score'] - 1 >= 0.25))

    return filter(ff_denovo, effects)


# ------------------------------------------------
def predict_lost_auth(chrom, pos, ref, alt, refseq=None, refseqgene=None, genepanel=None):
    '''
    Predicts whether the authentic splice site is conserved, lost,
    or the mutation has no effect on its nucleotide sequence.
    '''

    auth = rf.get_closest_authentic(chrom=chrom, pos=pos, refseqgene=refseqgene,
                                    genepanel=genepanel, refseq=refseq,
                                    get_sequence=True, seq_size=SEQ_HALF_SIZE)
    if not auth:
        return [{'effect_descr': EFFECT_KEYWORD['NOT_TRANSCRIBED']}]

    dist = pos - auth['pos']

    consensus_size = {
        ('Donor', '+'): [-2, 7],
        ('Acceptor', '+'): [-19, 4],
        ('Donor', '-'): [-3, 6],
        ('Acceptor', '-'): [-20, 3],
    }
    s, e = consensus_size[(auth['splice_type'], auth['strand'])]

    if auth['strand'] == '+':
        fasta = auth['fasta']
        fasta_mut = fasta[:SEQ_HALF_SIZE+dist] + alt + fasta[SEQ_HALF_SIZE+dist+len(ref):]
    elif auth['strand'] == '-':
        fasta = mes.reverse_complement(auth['fasta'])
        fasta_mut = fasta[:SEQ_HALF_SIZE-dist] + mes.reverse_complement(alt) + fasta[SEQ_HALF_SIZE-dist+len(ref):]

    wild = fasta[SEQ_HALF_SIZE+s:SEQ_HALF_SIZE+e]
    if 0 <= dist < SEQ_HALF_SIZE:
        mut = fasta_mut[SEQ_HALF_SIZE+s:SEQ_HALF_SIZE+e]
    elif -SEQ_HALF_SIZE < dist < 0:
        shift = -len(ref) + len(alt)
        mut = fasta_mut[SEQ_HALF_SIZE+shift+s:SEQ_HALF_SIZE+shift+e]
    else:
        mut = wild

    try:
        wild_score = mes.score(wild)
        mut_score = mes.score(mut)
    except Exception:
        return [{'effect_descr': EFFECT_KEYWORD['NOT_AVAILABLE']}]

    try:
        ratio = mut_score[0] / wild_score[0] - 1
    except ZeroDivisionError:
        if mut_score == 0.0:
            ratio = 0.0
        else:
            ratio = 100.0

    if wild == mut:
        effect_descr = EFFECT_KEYWORD['NO_EFFECT']  # the variant does not impact the site
    elif ratio <= THRESHOLD_LOST:
        effect_descr = EFFECT_KEYWORD['LOST']
    else:
        effect_descr = EFFECT_KEYWORD['CONSERVED']

    effect = {'effect_descr': effect_descr,
              'distance': dist,
              'wild_score': wild_score[0],
              'mut_score': mut_score[0],
              'wild_seq': wild,
              'mut_seq': mut,
              'auth_pos': auth['pos'],
              'splice_type': auth['splice_type'],
              'strand': auth['strand'],
              'transcript': auth['transcript']}

    return [effect]

# ------------------------------------------------
def annotate_2bp_conserved(chrom, pos, ref, alt, refseq=None, refseqgene=None, genepanel=None):
    '''
    Annotate if the conserved 2-base-pairs around the splice site are affected.
    This is similar to VEP where `splice_donor_variant` and
    `splice_acceptor_variant` is used for variants affecting GT/AG.
    Note that prediction is still handled by predict_lost_auth()
    '''

    auth = rf.get_closest_authentic(chrom=chrom, pos=pos, refseqgene=refseqgene,
                                    genepanel=genepanel, refseq=refseq,
                                    get_sequence=True, seq_size=SEQ_HALF_SIZE)
    if not auth:
        return []

    dist = pos - auth['pos']

    conserved_2bp = {
        ('Donor', '+'): [1, 2],
        ('Acceptor', '+'): [-1, 0],
        ('Donor', '-'): [-1, 0],
        ('Acceptor', '-'): [1, 2],
    }

    # import pdb; pdb.set_trace()

    if dist in conserved_2bp[(auth['splice_type'], auth['strand'])]:
        if auth['splice_type'] == 'Donor':
            effect_descr = EFFECT_KEYWORD['SPLICE_DONOR_VARIANT']
        else:
            effect_descr = EFFECT_KEYWORD['SPLICE_ACCEPTOR_VARIANT']

        effect = {'effect_descr': effect_descr,
                  'distance': dist,
                  'auth_pos': auth['pos'],
                  'splice_type': auth['splice_type'],
                  'strand': auth['strand'],
                  'transcript': auth['transcript']}
        return [effect]

    return []

# ------------------------------------------------
def print_vcf(effects):
    '''
    Here comes the VCF formatting
    '''

    p = list()
    for allele_effect in effects:
        s = list()
        for single_effect in allele_effect:
            if single_effect['effect_descr'] == EFFECT_KEYWORD['NOT_TRANSCRIBED']:
                s += [single_effect['effect_descr']]
            elif single_effect['effect_descr'] in [EFFECT_KEYWORD['NO_EFFECT'], EFFECT_KEYWORD['CONSERVED'], EFFECT_KEYWORD['LOST']]:
                s += ['|'.join([single_effect['transcript'],
                               single_effect['effect_descr'],
                               str(single_effect['wild_score']),
                               str(single_effect['mut_score'])
                              ])]

            elif single_effect['effect_descr'] == EFFECT_KEYWORD['DE_NOVO']:
                s += ['|'.join([single_effect['transcript'],
                                single_effect['effect_descr'],
                                str(single_effect['wild_score']),
                                str(single_effect['mut_score']),
                                str(single_effect['auth_score']),
                                str(single_effect['distance']),
                                ])]
            elif single_effect['effect_descr'] in [ EFFECT_KEYWORD['SPLICE_DONOR_VARIANT'], EFFECT_KEYWORD['SPLICE_ACCEPTOR_VARIANT'] ]:
                s += ['|'.join([single_effect['transcript'],
                                single_effect['effect_descr']
                                ])]
            else:
                s += ['NOT_IMPLEMENTED']
        p += ['&'.join(s)]

    return ','.join(p)


# ============================================================
def main():
    '''
    Testing
    '''
    from settings import settings
    REFSEQGENE = settings.concatBundlePath('funcAnnot/refseq/refGene_131119.tab')
    REFSEQ = settings.concatBundlePath(settings.getBundle()['reference']['fasta'])
    GENEPANEL = settings.concatBundlePath(settings.getBundle()['clinicalGenePanels']['Bindevev_v02']['transcripts'])

    records = [
        ('2', 162060108, 'T', 'A'),  # SNP after junction
        ('2', 162060108, 'T', 'AT'),  # indel after junction
        ('2', 162060108, '', 'A'),  # ins after junction
        ('2', 162060108, 'T', ''),  # del after junction
        ('2', 162060105, 'A', 'G'),  # SNP before junction
        ('2', 162060105, 'A', 'TT'), # indel before junction
        ('2', 162060105, '', 'G'),  # ins before junction
        ('2', 162060105, 'A', ''),  # del before junction
        ('2', 162060108, 'T', ['A', 'G']),  # SNP multiple alleles
        ('2', 162060228, 'T', 'A'),  # too far
        ('17', 41222943, 'A', 'C'),  # minus strand, SNP Donor after junction BRCA1:exon_15+2T>G
        ('17', 41222946, 'A', 'T'),  # minus strand, SNP Donor before junction BRCA1:exon_c15-2T>A
        ('17', 41222600, 'T', 'C'),  # too far from junction: no_effect
        ('2', 162060126, 'T', 'G'),  # intronic denovo donor
        ('2', 162060089, 'T', 'G'),  # exonic denovo donor
        ('2', 162060089, 'T', 'G'),  # exonic denovo donor
        ('2', 162060086, '', 'G'),  # exonic denovo donor
    ]
    for record in records:
        chrom, pos, ref, alt = record
        print record
        effects = predict(chrom, pos, ref, alt, refseq=refseq, refseqgene=refseqgene)
        print print_vcf(effects)

# ============================================================
if __name__ == "__main__":
    main()
