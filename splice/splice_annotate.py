# -*- coding: utf-8 -*-
'''
SpliceAnnotate

@author: Hugues Fontenelle, 2015
'''

__version__ = '0.2'

import argparse, sys, os
import vcf
from vcf.parser import _Info as VcfInfo
from splice import splice_predict

# ============================================================
def main():
        '''
        $ python splice_annotate.py --help
        usage: splice_annotate.py [-h] -i INPUT [-o OUTPUT] [--genepanel GENEPANEL]
                                  [--refseqgene REFSEQGENE] [--refseq REFSEQ]
        
        Annotate a VCF with splice site effect prediction.
        
        optional arguments:
          -h, --help            show this help message and exit
          -i INPUT, --input INPUT
                                Input VCF file
          -o OUTPUT, --output OUTPUT
                                Output VCF file for annotation of predictions
                                (default=same as input)
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
        vcf_info_desc = 'Splice effect. Format: Transcript|Effect|MaxEntScan-wild|MaxEntScan-mut'
        vcf_reader.infos['splice'] = VcfInfo(id='splice',
                                             num=1,
                                             type='String',
                                             desc=vcf_info_desc,
                                             source='AMG-OUS-splice',
                                             version=__version__)

        output_file = open('tmp.vcf', 'w')
        vcf_writer = vcf.Writer(output_file, vcf_reader)

        for record in vcf_reader:
            effect = splice_predict.predict(chrom=record.CHROM,
                                   pos=record.POS,
                                   ref=record.REF,
                                   alt=record.ALT,
                                   genepanel=args.genepanel,
                                   refseqgene=args.refseqgene,
                                   refseq=args.refseq)  

            if 'splice' in record.INFO:
                record.INFO['splice'] = effect
            else:
                record.add_info(effect)

            vcf_writer.write_record(record)

        vcf_handle.close()
        vcf_writer.close()

        os.rename('tmp.vcf', output_filename)


# ============================================================
if  __name__ == "__main__":
    sys.exit(main())
