# -*- coding: utf-8 -*-
'''
SpliceScore

@author: Hugues Fontenelle, 2014
'''

from Bio.Seq import Seq
from splice.SSFL import SSFL
from splice.max_ent_scan import MaxEntScan
from splice.gene_splicer import GeneSplicer
from splice.NNSplice import NNSplice
from splice.HSF import HSF
from splice.refseq_utils import *
import vcf, json
import argparse, os, sys


class SpliceScore(dict):
    '''
    This class fetches the FASTA sequence around the position of interest,
    mutate the sequence, and score both sequences with all the splice scoring
    algorithms. The scores are stored in JSON format.
    The closest authentic splice site (type, position, strand) is added.

    $ python splice_score.py -h
    usage: splice_score.py [-h] -i INPUT [-o OUTPUT] [--all] [--ssfl] [--mes]
                           [--gs] [--nn] [--hsf] [-v]

    Score cryptic splice sites around mutation position for each variant in VCF
    file. Outputs in a JSON file.

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Input VCF file
      -o OUTPUT, --output OUTPUT
                            Output JSON file (defaults to <input__name>.json
      --all                 Score with all splice site scoring algorithms
      --ssfl                Score with Splice Site Finder-Like (default) [Shapiro
                            et al. Nucleic Acids Research 15(17):7155-74, 1987]
      --mes                 Score with MaxEntScan (default) [Yeo and Burge. J
                            Comput Biol 11(2-3):377-94, 2004]
      --gs                  Score with GeneSplicer [Pertea et al. Nucleic Acids
                            Res 29(5):1185-90, 2001]
      --nn                  Score with NNSplice [Reese et al. J Comp Biol
                            4(3):311-23, 1997]
      --hsf                 Score with Human Splicing Finder [Desmet et al.
                            Nucleic Acids Research 37(9), 2009]
      -v, --verbose

    Example usage:
    $ python splice_score.py --vcf data/case_study0.vcf --mes --ssfl
    $ python splice_score.py --vcf data/case_study0.vcf --json data/case_study1.json --all
    '''

    _sub_seq_size = [-100, 99]

    # ------------------------------------------------
    def __init__(self, record={}):
        '''
        Loads a dictionary containing the variant
        Example:
        record = {'chrom':17, 'pos':41246536, 'ref':'T', 'alt':'A', 'ID':'CM057535'}
        set_dict(record)
        '''
        super(SpliceScore, self).__init__(record)

        self.use_SSFL = True
        self.use_MaxEntScan = True
        self.use_GeneSplicer = False
        self.use_NNSplice = False
        self.use_HSF = False

        self.wild = None # wild type sequence
        self.mut = None # mutated sequence
        self.wild_ss = {} # scored ss wild seq
        self.mut_ss = {} # scored ss mut seq
        
        self.gene_panel = None
        self.ref_seq_gene = None
        self.ref_seq = None

    # ------------------------------------------------
    def set_variant(self, chrom, pos, ref, alt, ID=None):
        '''
        Loads a variant by parameters
        Example:
        set_variant('17', 41246536, 'T', 'A', ID='CM057535')
        '''
        self['chrom'] = chrom
        self['pos'] = pos
        self['ref'] = ref
        self['alt'] = str(alt).strip('[]').split(',')[0]
        self['ID'] = ID
        
    # ------------------------------------------------
    def __str__(self):
        return '(%s) chr%s:%d%s>%s' % (self['ID'], self['chrom'], self['pos'], self['ref'], self['alt'])
        
    # ------------------------------------------------
    def set_gene_panel(self, genepanel):
        '''
        Set the path to the GenePanel file
        eg: PATH/TO/Bindevev_OUS_medGen_v01_b37.transcripts.csv
        '''
        self.gene_panel = genepanel
    
    # ------------------------------------------------
    def set_ref_seq_gene(self, refseqgene):
        '''
        Set the path to the RefSeqGene file
        eg: PATH/TO/refGene_131119.tab
        '''
        self.ref_seq_gene = refseqgene
        
    # ------------------------------------------------
    def set_ref_seq(self, refseq):
        '''
        Set the path to RefSeq file
        eq: PATH/TO/b37/human_g1k_v37_decoy.fasta
        '''
        self.ref_seq = refseq
                    
    # ------------------------------------------------
    def use_algo(self, use_SSFL=True, use_MaxEntScan=True,
                 use_GeneSplicer=False, use_NNSplice=False,
                 use_HSF=False):
        '''
        Sets which algorithm(s) to use, by setting the algorithm's name
        to True.
        By default, only SSFL and MaxEntScan are set to True.
        '''
        self.use_SSFL = use_SSFL
        self.use_MaxEntScan = use_MaxEntScan
        self.use_GeneSplicer = use_GeneSplicer
        self.use_NNSplice = use_NNSplice
        self.use_HSF = use_HSF

    # ------------------------------------------------
    def get_closest_authentic(self):
        if self.gene_panel:
            auth = get_closest_authentic(self['chrom'], self['pos'], genepanel=self.gene_panel)
        elif self.ref_seq_gene:
            auth = get_closest_authentic(self['chrom'], self['pos'], refseqgene=self.ref_seq_gene)
        else:
            auth = get_closest_authentic(self['chrom'], self['pos']) # NCBI internet connection     
        self['authentic'] = auth
        return auth

    # ------------------------------------------------------------
    def load_sequence(self):
        startpos = self['pos'] + self._sub_seq_size[0]
        endpos = self['pos'] + self._sub_seq_size[1]
        fasta = get_fasta(self['chrom'], startpos, endpos, refseq=self.ref_seq)
        self.wild = Seq(fasta)

    # ------------------------------------------------------------
    def mutate(self):
        '''Mutate wild seq to mut seq by replacing REF with ALT'''
        relpos = abs(self._sub_seq_size[0])
        self.mut = self.wild[0:relpos] + self['alt'] \
            + self.wild[relpos+len(self['ref']):]

    # ------------------------------------------------------------
    def append_algo_name(self, sites, algo_name):
        '''
        Appends the name of the algorithm to the splice site dictionary
        Input dict {(pos, splice_type, strand): score}
        Output dict {(pos, splice_type, strand, algo_name): score}
        '''
        out = {}
        for key, value in sites.iteritems():
            out[(key[0], key[1], key[2], algo_name)] = value
        return out

    # ------------------------------------------------------------
    def score_splice_sites(self):
        '''
        Score both the wild and mut sequences with all selected algorithms
        '''
        self.load_sequence()
        self.mutate()
        self.get_closest_authentic()

        splice_type = self['authentic']['splice_type']
        strand = self['authentic']['strand']
        
        offset = self['pos'] + self._sub_seq_size[0] - 1
        self.wild_ss = {}
        self.mut_ss = {}
        if self.use_SSFL:
            wild_ss = SSFL(self.wild, offset).find_all_sites(splice_type=splice_type, strand=strand)
            mut_ss = SSFL(self.mut, offset).find_all_sites(splice_type=splice_type, strand=strand)
            self.wild_ss.update(self.append_algo_name(wild_ss, 'SSFL'))
            self.mut_ss.update(self.append_algo_name(mut_ss, 'SSFL'))
        if self.use_MaxEntScan:
            wild_ss = MaxEntScan(self.wild, offset).find_all_sites(splice_type=splice_type, strand=strand)
            mut_ss = MaxEntScan(self.mut, offset).find_all_sites(splice_type=splice_type, strand=strand)
            self.wild_ss.update(self.append_algo_name(wild_ss, 'MaxEntScan'))
            self.mut_ss.update(self.append_algo_name(mut_ss, 'MaxEntScan'))
        if self.use_GeneSplicer:
            wild_ss = GeneSplicer(self.wild, offset).find_all_sites(splice_type=splice_type, strand=strand)
            mut_ss = GeneSplicer(self.mut, offset).find_all_sites(splice_type=splice_type, strand=strand)
            self.wild_ss.update(self.append_algo_name(wild_ss, 'GeneSplicer'))
            self.mut_ss.update(self.append_algo_name(mut_ss, 'GeneSplicer'))
        if self.use_NNSplice:
            wild_ss = NNSplice(self.wild, offset).find_all_sites(splice_type=splice_type, strand=strand)
            mut_ss = NNSplice(self.mut, offset).find_all_sites(splice_type=splice_type, strand=strand)
            self.wild_ss.update(self.append_algo_name(wild_ss, 'NNSplice'))
            self.mut_ss.update(self.append_algo_name(mut_ss, 'NNSplice'))
        if self.use_HSF:
            wild_ss = HSF(self.wild, offset).find_all_sites(splice_type=splice_type, strand=strand)
            mut_ss = HSF(self.mut, offset).find_all_sites(splice_type=splice_type, strand=strand)
            self.wild_ss.update(self.append_algo_name(wild_ss, 'HSF'))
            self.mut_ss.update(self.append_algo_name(mut_ss, 'HSF'))

        # store wild scores in dict
        wild = {}
        for k, v in self.wild_ss.iteritems():
            key = (k[0], k[1], k[2])
            value = {k[3] : v}
            if key not in wild:
                wild[key] = value
            else:
                wild[key].update(value)
        self['wild'] = [{"pos":k[0], "scores":v} for k, v in wild.iteritems()]
        # store mut scores in dict
        mut = {}
        for k, v in self.mut_ss.iteritems():
            key = (k[0], k[1], k[2])
            value = {k[3] : v}
            if key not in mut:
                mut[key] = value
            else:
                mut[key].update(value)
        self['mut'] = [{"pos":k[0], "scores":v} for k, v in mut.iteritems()]

# ============================================================
def main():
    '''
    Main.
    Usage:
    $ python splice_score.py --vcf data/case_study0.vcf --o data/case_study0.json
    '''
    parser = argparse.ArgumentParser(description=
        "Score cryptic splice sites around mutation \
        position for each variant in VCF file.\
        Outputs in a JSON file.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input VCF file')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='Output JSON file (defaults to <input__name>.json')
    parser.add_argument('--all', action='store_true',
                        default=False,
                        help='Score with all splice site scoring algorithms')
    parser.add_argument('--ssfl', action='store_true',
                        default=True,
                        help='Score with Splice Site Finder-Like (default) \
                        [Shapiro et al. Nucleic Acids Research 15(17):7155-74, 1987]')
    parser.add_argument('--mes', action='store_true',
                        default=True,
                        help='Score with MaxEntScan (default) \
                        [Yeo and Burge. J Comput Biol 11(2-3):377-94, 2004]')
    parser.add_argument('--gs', action='store_true',
                        default=False,
                        help='Score with GeneSplicer \
                        [Pertea et al. Nucleic Acids Res 29(5):1185-90, 2001]')
    parser.add_argument('--nn', action='store_true',
                        default=False,
                        help='Score with NNSplice \
                        [Reese et al. J Comp Biol 4(3):311-23, 1997]')
    parser.add_argument('--hsf', action='store_true',
                        default=False,
                        help='Score with Human Splicing Finder \
                        [Desmet et al. Nucleic Acids Research 37(9), 2009]')
    parser.add_argument('--genepanel', type=str, required=False,
                        help='Filepath for genepanel')
    parser.add_argument('--refseqgene', type=str, required=False,
                        help='Filepath for refSeqGene')
    parser.add_argument('--refseq', type=str, required=False,
                        help='Filepath for refSeq reference FASTA sequences')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)

    args = parser.parse_args()
    vcf_filename = args.input
    if not args.output:
        json_filename = os.path.splitext(vcf_filename)[0] + '.json'
    else:
        json_filename = args.output
    verbose = args.verbose
    use_SSFL = args.ssfl
    use_MaxEntScan = args.mes
    use_GeneSplicer = args.gs
    use_NNSplice = args.nn
    use_HSF = args.hsf
    if args.all:
        use_SSFL = True
        use_MaxEntScan = True
        use_GeneSplicer = True
        use_NNSplice = True
        use_HSF = True

    ori_path = os.getcwd()

    vcf_handle = open(vcf_filename, 'r')
    vcf_reader = vcf.Reader(vcf_handle)
    records = list()
    for vcf_record in vcf_reader:
        if verbose:
            print "Processing %s: chr%s:%d%s>%s" % (vcf_record.ID,
                                                    vcf_record.CHROM,
                                                    vcf_record.POS,
                                                    vcf_record.REF,
                                                    vcf_record.ALT)
        splice = SpliceScore()
        if args.genepanel:
            splice.set_gene_panel(args.genepanel)
        if args.refseqgene:
            splice.set_ref_seq_gene(args.refseqgene)
        if args.refseq:
            splice.set_ref_seq(args.refseq)
        splice.set_variant(vcf_record.CHROM, vcf_record.POS,
                           vcf_record.REF, vcf_record.ALT, ID=vcf_record.ID)
        splice.use_algo(use_SSFL=use_SSFL,
                        use_MaxEntScan=use_MaxEntScan,
                        use_GeneSplicer=use_GeneSplicer,
                        use_NNSplice=use_NNSplice,
                        use_HSF=use_HSF)
        splice.score_splice_sites()
        records.append(splice)

    vcf_handle.close()

    json_handle = open(json_filename, 'w')
    json_content = json.dumps(records, indent = 4, ensure_ascii=False)
    json_handle.write(json_content)
    json_handle.close()

    os.chdir(ori_path)

# ============================================================
if  __name__ == "__main__":
    sys.exit(main())

'''
# Example code "inline" to copy in the console
from splice.splice_score import SpliceScore
REFSEQGENE = "/Users/huguesfo/Documents/DATA/b37/refSeq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/Documents/DATA/b37/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
record = {'chrom':'17', 'pos':41246536, 'ref':'T', 'alt':'A', 'ID':'CM057535'}
splice = SpliceScore(record)
splice.set_ref_seq(REFSEQ)
splice.set_ref_seq_gene(REFSEQGENE)
splice.use_algo(use_SSFL=True, use_MaxEntScan=True)
splice.score_splice_sites()
print splice
'''

