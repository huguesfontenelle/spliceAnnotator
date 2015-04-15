'''
SplicePredict

@author: Hugues Fontenelle, 2014
'''

import argparse, os, sys
import json
from splice.splice_score import SpliceScore

STRATEGIES = {
    "Jian14nar": {
        "description": "MaxEntScan -21.7%. Any new cryptic site",
        "reference": "Jian, X., Boerwinkle, E., Liu, X. \
In silico prediction of splice-altering single nucleotide variants in the human genome. \
Nucleic Acids Research, 42(22), 13534\342\200\22313544. doi:10.1093/nar/gku1206",
        "algorithms" : ["MaxEntScan"],
        "MaxEntScan": -.217, "consensus": 1},
    "Houdayer": {
        "description": "MaxEntScan -15% and SSFL -5%. Any new cryptic site",
        "reference": "Houdayer et al., Guidelines for Splicing Analysis in Molecular Diagnosis Derived from a Set of 327 Combined In Silico-In Vitro Studies on BRCA1 and BRCA2 Variants. Hum Mutat 33:1228-1238, 2012",
        "algorithms" : ["MaxEntScan", "SSFL"],
        "MaxEntScan": -.15, "SSFL": -.05, "consensus": 2},
    "AMG-diag": {
        "description": "3 algorithms among 5 (SSFL, MaxEntScan, NNSplice, GeneSplicer and HSF) agree that splice site strength is decreased (any decrease). Any new cryptic site agreed by 3 algorithms.",
        "algorithms": ["MaxEntScan", "SSFL", "NNSplice", "GeneSplicer", "HSF"],
        "MaxEntScan": -0.0,
        "SSFL": -0.0,
        "NNSplice": -0.0,
        "GeneSplicer": -0.0,
        "HSF": -0.0,
        "consensus": 3},
    "AMG-kreftgenetikk": {
        "description": "3 algorithms among 4 (SSFL, MaxEntScan, NNSplice and GeneSplicer) agree that splice site strength is decreased by at least -10%. Any new cryptic site agreed by 3 algorithms.",
        "reference": "Hellen, B., Splice site tools: A Comparative Analysis Report. NGRL, 2009",
        "algorithms": ["MaxEntScan", "SSFL", "NNSplice", "GeneSplicer"],
        "MaxEntScan": -0.10,
        "SSFL": -0.10,
        "NNSplice": -0.10,
        "GeneSplicer": -0.10,
        "consensus": 3}
    }

class SplicePredict(dict):
    '''
    Class for predicting the effect of a variant on splicing.

    $ python splice_predict.py --help
    usage: splice_predict.py [-h] -i INPUT -o OUTPUT [-m METHOD]

    Predict splice sites around mutation position for each variant stored in JSON.

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Input JSON file with precomputed scores
      -o OUTPUT, --output OUTPUT
                            Output JSON file for predictions (default=same as
                            input)
      -m METHOD, --method METHOD
                            Prediction method: choose among AMG-diag,
                            AMG-kreftgenetikk or Houdayer(default)

    Example usage:
    $ python splice_predict.py -i data/case_study0.json
    '''

    # ------------------------------------------------
    def __init__(self, record={}):
        '''
        Loads a dictionary containing the variant. The dictionary may contain the precomputed scores.
        Example:
        record = {'chrom':17, 'pos':41246536, 'ref':'T', 'alt':'A', 'ID':'CM057535'}
        set_dict(record)
        '''
        super(SplicePredict, self).__init__(record)

        self.strategy_default = 'Jian14nar'
        self.effect_dict = {'no_effect' : 'no_effect',
                            'lost_splice_site' : 'lost_site',
                            'new_cryptic_splice_site' : 'de_novo'}
        self.effect_default = 'no_effect'
        self.strategy = self.strategy_default

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
    def get_changes(self, change=0.0):
        '''
        Get the cryptic splice sites that differ between wild and mut.
        Filters by change % (optional; default=0 means that any change is
        reported)
        Output format is dict
        {(pos, splice_type, strand, algo_name) :
            [wild_score, mut_score, change]}
        '''
        wild_set = set(self.wild_ss) # set of wild splice sites
        mut_set = set(self.mut_ss) # set of mut splice sites
        out = {}
        common_ss = wild_set & mut_set # intersection: elements in both
        for splice_site in common_ss:
            wild_score = self.wild_ss[splice_site]
            mut_score = self.mut_ss[splice_site]
            ratio = (mut_score / float(wild_score)) - 1
            if abs(ratio) > change:
                out[splice_site] = [wild_score, mut_score, ratio]
        return out

    # ------------------------------------------------
    def get_new(self):
        '''
        Get the new cryptic splice sites that appear in mut but not in wild.
        Output format is dict
        {(pos, splice_type, strand, algo_name) : mut_score}
        '''
        wild_set = set(self.wild_ss) # set of wild splice sites
        mut_set = set(self.mut_ss) # set of mut splice sites
        out = {}
        new_ss = mut_set - wild_set # elements of mut not in wil: new ss
        for splice_site in new_ss:
            mut_score = self.mut_ss[splice_site]
            out[splice_site] = mut_score
        return out

    # ------------------------------------------------
    def get_lost(self):
        '''
        Get the lost cryptic splice sites that appear in wild but not in mut.
        Output format is dict
        {(pos, splice_type, strand, algo_name) : wild_score}
        '''
        wild_set = set(self.wild_ss) # set of wild splice sites
        mut_set = set(self.mut_ss) # set of mut splice sites
        out = {}
        lost_ss = wild_set - mut_set # elements of wild not in mut: lost ss
        for splice_site in lost_ss:
            wild_score = self.wild_ss[splice_site]
            out[splice_site] = wild_score
        return out

    # ------------------------------------------------
    def filter_sites(self, sites, *args, **kwargs):
        '''
        Filter sites by splice_type, strand and/or greater than threshold

        Usage examples:
        $ filtered = filter_sites(sites, splice_type='Acceptor', strand='-',
                                  threshold=60.0)
        $ filtered = filter_sites(sites, splice_type='Donor')
        $ filtered = filter_sites(sites, strand='+')
        $ filtered = filter_sites(sites, threshold=60.0)
        $ filtered = filter_sites(sites, algo_name='MaxEntScan')

        Notes:
        1) The 'sites' data structure is a dictionary where the keys are
           tuples of (position, splice_type, strand); the value is the score.
           $ sites = {(3, 'Donor', '+'): 100.0,
                      (65, 'Acceptor', '+'): 100.0}
        2) When filtering by algo_name, the returned tuple does not contain
           the algo_name anymore.
        '''
        strand = kwargs.get('strand', None)
        splice_type = kwargs.get('splice_type', None)
        threshold = kwargs.get('threshold', None)
        algo_name = kwargs.get('algo_name', None)
        algo_names = kwargs.get('algo_names', None)

        if strand not in ['+', '-', None]:
            raise ValueError("%r is not a proper strand orientation. \
                              Specify '+' or '-'." % strand)
        if splice_type not in ['Donor', 'Acceptor', None]:
            raise ValueError("%r is not a proper splice type. Specify 'Donor' \
                              or 'Acceptor'." % splice_type)
        if algo_name not in ['SSFL', 'MaxEntScan', 'GeneSplicer', 'NNSplice', None]:
            raise ValueError("%r is not a algorithm name for splice scoring.\
                              Specify 'SSFL', 'MaxEntScan', 'GeneSplicer', \
                              'NNSplice' or 'HSF'." % algo_name)

        filtered = sites

        if splice_type is not None:
            filtered = {k: v for k, v in filtered.iteritems()
                        if k[1] == splice_type}
        if strand is not None:
            filtered = {k: v for k, v in filtered.iteritems()
                        if k[2] == strand}
        if threshold is not None:
            filtered = {k: v for k, v in filtered.iteritems()
                        if v >= threshold}
        if algo_name is not None:
            filtered = {(k[0], k[1], k[2]): v for k, v in filtered.iteritems()
                        if k[3] == algo_name}
        if algo_names is not None:
            filtered = {(k[0], k[1], k[2], k[3]): v for k, v in filtered.iteritems()
                        if k[3] in algo_names}
        return filtered

    # ------------------------------------------------
    def predict(self):
        '''
        Filters the wild, mut, changed, new, lost ss
        with splice_type and strand info, as well as selected algorithm.
        Predicts
        '''
        if self.strategy not in STRATEGIES:
            raise RuntimeError('Unknown method.')

        # TODO: currently does SNP only
        if (len(self['ref']) != 1) or (len(self['alt']) != 1):
            return {self.strategy: [{'Effect':'N/A'}]}

        ####################
        # check if SCORES available
        ####################
        if 'authentic' not in self or 'wild' not in self or 'mut' not in self:
            scores = SpliceScore(self)
            scores.set_gene_panel(self.gene_panel)
            scores.set_ref_seq_gene(self.ref_seq_gene)
            scores.set_ref_seq(self.ref_seq)
            
            if self.strategy not in STRATEGIES:
                raise Exception("Cannot score: Unknown strategy %s" % self.strategy)
            use_strategy = STRATEGIES[self.strategy]['algorithms']
            scores.use_algo(use_SSFL='SSFL' in use_strategy,
                            use_MaxEntScan='MaxEntScan' in use_strategy,
                            use_NNSplice='NNSplice' in use_strategy,
                            use_GeneSplicer='GeneSplicer' in use_strategy,
                            use_HSF='HSF' in use_strategy)

            scores.score_splice_sites()
            self['authentic'] = scores['authentic']
            self['wild'] = scores['wild']
            self['mut'] = scores['mut']

        ####################
        # load AUTH
        ####################
        closest_real = self['authentic']['pos']
        splice_type = self['authentic']['splice_type']
        strand = self['authentic']['strand']
        self.auth_ss = (closest_real, splice_type, strand)

        ####################
        # convert SCORES
        ####################
        # get wild {(pos, splice_type, strand, algo): score}
        wild = {}
        for ss in self['wild']:
            for k, v in ss['scores'].iteritems():
                wild.update({(ss['pos'], ss.get('splice_type', splice_type), ss.get('strand', strand), k) : v})
        self.wild_ss = wild

        # get mut {(pos, splice_type, strand, algo): score}
        mut = {}
        for ss in self['mut']:
            for k, v in ss['scores'].iteritems():
                mut.update({(ss['pos'], ss.get('splice_type', splice_type), ss.get('strand', strand), k) : v})
        self.mut_ss = mut

        ####################
        # FILTER
        ####################
        changed_ss = self.get_changes()
        new_ss = self.get_new()
        lost_ss = self.get_lost()
        # 'lost' must be part of 'changed' because they are changed by -100%
        changed_ss.update({k:[v, 0.0, -1.0] for k, v in lost_ss.iteritems()})

        # CHANGED
        # filter by splice_type and strand
        changed_ss = self.filter_sites(changed_ss,
                                       splice_type=splice_type,
                                       strand=strand)
        # keep only relevant algorithms
        changed_ss = self.filter_sites(changed_ss,
                                       algo_names=STRATEGIES[self.strategy]['algorithms'])
        # filter by x% decrease:
        changed_ss = {k: v for k, v
                      in changed_ss.iteritems() if v[2] < STRATEGIES[self.strategy][k[3]]}
        # reformat the changed dictionnary to
        # unique keys (pos, splice_type, strand)
        # and values as tuples in a list
        # [(algo_name, score_wild, score_mut, change)]
        changed_list = {}
        for key, value in changed_ss.iteritems():
            (pos, splice_type, strand, algo_name) = key
            (score_wild, score_mut, change) = value
            site = (pos, splice_type, strand)
            if not site in changed_list:
                changed_list[site] = [(algo_name, score_wild,
                                       score_mut, change)]
            else:
                changed_list[site].append((algo_name, score_wild,
                                           score_mut, change))
        # now create dict when there is 'maj_no' or more value tuples
        self.changed_ss = {k:v for k, v in changed_list.iteritems()
                           if len(v) >= STRATEGIES[self.strategy]['consensus']}
        # NEW
        new_ss = self.filter_sites(new_ss, splice_type=splice_type,
                                   strand=strand)
        # keep MES, SSFL, NNSplice, GeneSplicer
        new_ss = self.filter_sites(new_ss,
                                   algo_names=STRATEGIES[self.strategy]['algorithms'])
        new_list = {}
        for key, score in new_ss.iteritems():
            (pos, splice_type, strand, algo_name) = key
            site = (pos, splice_type, strand)
            if not site in new_list:
                new_list[site] = [(algo_name, score)]
            else:
                new_list[site].append((algo_name, score))
        # now create dict when there is 'maj_no' or more value tuples
        self.new_ss = {k:v for k, v in new_list.iteritems()
                       if len(v) >= STRATEGIES[self.strategy]['consensus']}


        ####################
        # PREDICT
        ####################
        splice_effect = []

        # ............................................................
        # authentic splice site
        # is the authetic splice site in the wild ss?
        dist_d = {('Donor', '+'):[-3, 6],
                      ('Donor', '-'):[-6, 3],
                      ('Acceptor', '+'):[-13, 4],
                      ('Acceptor', '-'):[-4, 13]}
        try:
            dist = self.pos - self.auth_ss[0]

            dd = dist_d[(splice_type, strand)]
        except:
            dist = float('Inf')
            dd = dist_d[('Donor', '+')]

        if dist >= dd[0] or dist <= dd[1]:
            # is authentic ss changed?
            if self.auth_ss in self.changed_ss:
                # *** lost_site ***
                # the real splice site is significantly decreased
                effect = {'Effect': self.effect_dict['lost_splice_site']}
                scores = {}
                for key in self.changed_ss[self.auth_ss]:
                    algo_name, wild_real_score, mut_real_score, ratio = key
                    scores.update({algo_name: [wild_real_score, mut_real_score]})
                effect.update({'scores':scores})
                splice_effect.append(effect)

        # ............................................................
        # cryptic splice sites
        # change
        # (irrelevant)
        # new (can be useful)
        for (pos, splice_type, strand) in self.new_ss:
            # *** de_novo ***
            effect = {'Effect': self.effect_dict['new_cryptic_splice_site']}
            effect.update({'pos': pos})
            scores = {}
            for algo_name, score in self.new_ss[(pos, splice_type, strand)]:
                scores.update({algo_name: [0, score]})
            effect.update({'scores':scores})
            splice_effect.append(effect)
        # lost
        # (irrelevant)

        if not splice_effect:
            # *** no_effect ***
            effect={'Effect': self.effect_dict['no_effect']}
            # Also report scores that do not lead to lost_site
            if self.auth_ss[0] in [k[0] for k,v in self.wild_ss.iteritems()]:
                scores = {}
                for algo_name in STRATEGIES[self.strategy]['algorithms']:
                    wild_score = self.wild_ss.get( self.auth_ss+tuple([algo_name]), 0.0)
                    mut_score = self.mut_ss.get( self.auth_ss+tuple([algo_name]), 0.0)
                    scores.update({algo_name: [wild_score, mut_score]})
                effect.update({'scores':scores})
            splice_effect.append(effect)

        if 'predict' in self:
            self['predict'].update({self.strategy: splice_effect})
        else:
            self['predict'] = {self.strategy: splice_effect}

        return {self.strategy: splice_effect}

# ============================================================
def main():
    '''
    Main.
    Usage:
    $ python splice_predict.py -i data/case_study0.json
    '''
    ori_path = os.getcwd()

    parser = argparse.ArgumentParser(description="Predict splice sites \
around mutation position for each variant stored in JSON.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input JSON file with precomputed scores')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='Output JSON file for predictions (default=same as input)')
    parser.add_argument('-m', '--method', type=str, required=False,
                        help='Prediction method: choose among ''Jian14nar''(default), ''AMG-diag'', ''AMG-kreftgenetikk'' or ''Houdayer''')
    parser.add_argument('--genepanel', type=str, required=False,
                        help='Filepath for genepanel')
    parser.add_argument('--refseqgene', type=str, required=False,
                        help='Filepath for refSeqGene')
    parser.add_argument('--refseq', type=str, required=False,
                        help='Filepath for refSeq reference FASTA sequences')
    args = parser.parse_args()
    if not args.method:
        method = 'Jian14nar'
    else:
        method = args.method
    input_filename = args.input
    if not args.output:
        output_filename = input_filename
    else:
        output_filename = args.output

    records = list()
    with open(input_filename, 'r') as f:
        for record in json.load(f):
            predict = SplicePredict(record)
            if args.genepanel:
                predict.set_gene_panel(args.genepanel)
            if args.refseqgene:
                predict.set_ref_seq_gene(args.refseqgene)
            if args.refseq:
                predict.set_ref_seq(args.refseq)
            predict.strategy = method
            predict.predict()
            records.append(predict)

    with open(output_filename, 'w') as f:
        f.write(json.dumps(records, indent=4, ensure_ascii=False))

    os.chdir(ori_path)

# ============================================================
if  __name__ == "__main__":
    sys.exit(main())

'''
# ------------------------------------------------------------
#from splice.splice_predict import SplicePredict
REFSEQGENE = "/Users/huguesfo/Devel/genevar/vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab" # RefSeqGene definitions
REFSEQ = "/Users/huguesfo/Devel/genevar/vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta" # RefSeq FASTA sequences (hg19)
record = {'chrom':'17', 'pos':41246536, 'ref':'T', 'alt':'A', 'ID':'CM057535'}
p1 = SplicePredict(record)
p1.set_ref_seq(REFSEQ)
p1.set_ref_seq_gene(REFSEQGENE)
p1.strategy = 'Jian14nar'
p1.predict()
print p1

from splice.splice_score import SpliceScore
score = SpliceScore(record)
score.set_ref_seq(REFSEQ)
score.set_ref_seq_gene(REFSEQGENE)
score.use_algo(use_SSFL=True, use_MaxEntScan=True, use_NNSplice=True, use_GeneSplicer=True)
score.score_splice_sites()
p2 = SplicePredict(score)
p2.strategy = 'Jian14nar'
p2.predict()
print p2

assert p1 == p2
'''


