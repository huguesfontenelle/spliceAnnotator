'''
SplicePredict

@author: Hugues Fontenelle, 2014
'''


import os.path
from splice import max_ent_scan as mes
from splice import refseq_utils as rf

THRESHOLD_LOST = -0.216 # threshold at which a drop in MES score is considered "LOST"
SEQ_HALF_SIZE = 50 # half-length of the FASTA sequence retrieved around the splice site

# ------------------------------------------------
def predict(chrom, pos, ref, alt, refseq=None, refseqgene=None, genepanel=None):
    '''
    Predicts
    '''
    assert os.path.isfile(refseq)
    assert os.path.isfile(refseqgene) or os.path.isfile(genepanel)

    if type(alt) is str:
        alt = [alt]

    effects = list()        
    for alt1 in alt:
        effect = predict_one(chrom, pos, ref, str(alt1), refseq=refseq, refseqgene=refseqgene, genepanel=genepanel)
        effects += [print_vcf(effect)]
    
    return ','.join(effects)
            
        

# ------------------------------------------------
def predict_one(chrom, pos, ref, alt, refseq=None, refseqgene=None, genepanel=None):
    '''
    Predicts
    '''
    assert os.path.isfile(refseq)
    assert os.path.isfile(refseqgene) or os.path.isfile(genepanel)
    
    auth = rf.get_closest_authentic(chrom=chrom, pos=pos, refseqgene=refseqgene, genepanel=genepanel, refseq=refseq, get_sequence=True, seq_size=SEQ_HALF_SIZE)
    
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
    else:
        return {'effect_descr': 'not_in_transcript'}
    
    wild = fasta[SEQ_HALF_SIZE+s:SEQ_HALF_SIZE+e]
    mut = fasta_mut[SEQ_HALF_SIZE+s:SEQ_HALF_SIZE+e]

    wild_score = mes.score(wild)
    mut_score = mes.score(mut)
    
    ratio = mut_score[0] / wild_score[0] - 1
    
    if wild == mut:
        effect_descr = 'no_effect'  # the variant does not impact the site
    elif ratio <= THRESHOLD_LOST:
        effect_descr = 'lost_site'
    else:
        effect_descr  = 'conserved_site'

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

    return effect

# ------------------------------------------------
def print_vcf(effect):
    '''
    Here comes the VCF formatting
    '''
    if effect['effect_descr'] == 'not_in_transcript':
        return effect['effect_descr']
        
    s = '|'.join([effect['transcript'], effect['effect_descr'], str(effect['wild_score']), str(effect['mut_score'])])
    
    return s


# ============================================================
def main():
    '''
    Testing
    '''
    refseqgene = "/Users/huguesfo/genevar/vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab" # RefSeqGene definitions
    refseq = "/Users/huguesfo/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
    genepanel = "/Users/huguesfo/genevar/vcpipe-bundle/clinicalGenePanels/Bindevev_v02/Bindevev_v02.transcripts.csv" # gene panel transcript file

    records = [
               ('2', 162060108, 'T', 'A'), #SNP after junction
               ('2', 162060108, 'T', 'AT'), #indel after junction
               ('2', 162060108, '', 'A'), #ins after junction
               ('2', 162060108, 'T', ''), #del after junction
               ('2', 162060105, 'A', 'G'), #SNP before junction
               ('2', 162060105, 'A', 'TT'), #indel before junction
               ('2', 162060105, '', 'G'), #ins before junction
               ('2', 162060105, 'A', ''), #del before junction
               ('2', 162060108, 'T', ['A', 'G']), #SNP multiple alleles
               ('2', 162060228, 'T', 'A'), #too far
               ('17', 41222943, 'A', 'C'), # minus strand, SNP Donor after junction BRCA1:exon_15+2T>G    
               ('17', 41222946, 'A', 'T'), # minus strand, SNP Donor before junction BRCA1:exon_c15-2T>A   
               ('17', 41222600, 'T', 'C'), # too far from junction: no_effect
               ]
    for record in records:
        chrom, pos, ref, alt = record
        print record
        print predict(chrom, pos, ref, alt, refseq=refseq, refseqgene=refseqgene)


# ============================================================
if __name__ == "__main__":
    main()
