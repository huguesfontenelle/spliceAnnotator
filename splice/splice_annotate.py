# -*- coding: utf-8 -*-
'''
SpliceAnnotate

@author: Hugues Fontenelle, 2015
'''

__version__ = '0.1'

import argparse, sys, os
import vcf, json
from vcf.parser import _Info as VcfInfo
from splice.splice_predict import SplicePredict

STRATEGY = 'Jian14nar'

# ============================================================
def main():
        '''
        $ python ../../annotation/splice/splice_annotate.py --help
        usage: splice_annotate.py [-h] -i INPUT [-o OUTPUT] [--json JSON]
                                  [--genepanel GENEPANEL] [--refseqgene REFSEQGENE]
                                  [--refseq REFSEQ]

        Annotate a VCF with splice site effect prediction.

        optional arguments:
          -h, --help            show this help message and exit
          -i INPUT, --input INPUT
                                Input VCF file
          -o OUTPUT, --output OUTPUT
                                Output VCF file for annotation of predictions
                                (default=same as input)
          --json JSON           Output JSON file with scores and predictions
          --genepanel GENEPANEL
                                Filepath for genepanel
          --refseqgene REFSEQGENE
                                Filepath for refSeqGene
          --refseq REFSEQ       Filepath for refSeq reference FASTA sequences


        Example usage:
        $ python splice_annotate.py -i data/case_study0.vcf
        $ python splice_annotate.py -i data/case_study0.vcf -o data/case_study0_annot.vcf
        $ python splice_annotate.py -i data/case_study0.vcf --genepanel HBOC_OUS_medGen_v00_b37.transcripts.csv --refseq b37/human_g1k_v37_decoy.fasta
        $ python splice_annotate.py -i data/case_study0.vcf --refseqgene b37/refSeq/refGene_131119.tab --refseq b37/human_g1k_v37_decoy.fasta
        '''
        parser = argparse.ArgumentParser(description='Annotate a VCF with splice site effect prediction.')
        parser.add_argument('-i', '--input', type=str, required=True,
                            help='Input VCF file')
        parser.add_argument('-o', '--output', type=str, required=False,
                            help='Output VCF file for annotation of predictions (default=same as input)')
        parser.add_argument('--json', type=str, required=False,
                            help='Output JSON file with scores and predictions')
        parser.add_argument('--genepanel', type=str, required=False,
                            help='Filepath for genepanel')
        parser.add_argument('--refseqgene', type=str, required=False,
                            help='Filepath for refSeqGene')
        parser.add_argument('--refseq', type=str, required=False,
                            help='Filepath for refSeq reference FASTA sequences')

        args = parser.parse_args()
        input_filename = args.input
        if not args.output:
            output_filename = input_filename
        else:
            output_filename = args.output

        vcf_handle = open(input_filename, 'r')
        vcf_reader = vcf.Reader(vcf_handle)
        vcf_info_desc = 'Splice effect (method=%s). Format: Transcript|Effect|MaxEntScan-wild|MaxEntScan-mut' % STRATEGY
        vcf_reader.infos['splice'] = VcfInfo(id='splice',
                                             num=1,
                                             type='String',
                                             desc=vcf_info_desc,
                                             source='AMG-OUS-splice',
                                             version=__version__)

        output_file = open('tmp.vcf', 'w')
        vcf_writer = vcf.Writer(output_file, vcf_reader)

        records = list()
        for record in vcf_reader:
            predict = SplicePredict()
            if args.genepanel:
                predict.set_gene_panel(args.genepanel)
            if args.refseqgene:
                predict.set_ref_seq_gene(args.refseqgene)
            if args.refseq:
                predict.set_ref_seq(args.refseq)
            predict.set_variant(record.CHROM, record.POS, record.REF, record.ALT, ID=record.ID)
            predict.strategy = STRATEGY
            splice_effect = predict.predict()
            records.append(predict)

            # write in VCF
            splice_annot = []
            transcript = predict.get('authentic', {}).get('transcript', '')
            for effect in splice_effect[predict.strategy]:
                effect_name = effect.get('Effect', 'unknown')
                scores = effect.get('scores', {})
                MES = [round(x, 2) if x != '' else '' for x in scores.get('MaxEntScan', ['', ''])]
                #SSFL = [round(x, 2) if x != '' else '' for x in scores.get('SSFL', ['', ''])]
                if effect_name is not 'de_novo':
                    splice_annot.append('%s|%s|%s|%s' % (transcript, effect_name, MES[0], MES[1]))
                else:
                    pos = effect.get('pos', '?')
                    splice_annot.append('%s|%s@%s|%s|%s' % (transcript, effect_name, pos, MES[0], MES[1]))
            if record.INFO.has_key('splice'):
                record.INFO['splice'] = '&'.join(splice_annot)
            else:
                record.add_info('splice', '&'.join(splice_annot))
            vcf_writer.write_record(record)

        vcf_handle.close()
        vcf_writer.close()

        os.rename('tmp.vcf', output_filename)

        if args.json:
            with open(args.json, 'w') as f:
                f.write(json.dumps(records, indent=4, ensure_ascii=False))

# ============================================================
if  __name__ == "__main__":
    sys.exit(main())
