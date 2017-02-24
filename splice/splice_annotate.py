# -*- coding: utf-8 -*-
'''
SpliceAnnotate

@author: Hugues Fontenelle, 2015
'''

__version__ = '0.4'

import argparse, sys, os
import vcf
import tempfile, shutil
from vcf.parser import _Info as VcfInfo
from splice import splice_predict
from multiprocessing import Pool
import StringIO

HEADER_INFO = '##INFO=<ID=splice,Number=1,Type=String,Description="Splice effect. Format: Transcript|Effect|MaxEntScan-wild|MaxEntScan-mut|MaxEntScan-closest|dist">'


# ============================================================
if __name__ == "__main__":
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

        with open(output_filename, 'w') as vcf_output:

            result_list = []

            def flush(result_list):
                for line in result_list:
                    vcf_output.write(line)

            def log_result(line):
                global result_list
                result_list.append(line)
                if len(result_list) >= 100:
                    flush(result_list)
                    result_list = []

            def process(line):
                elems = line.split('\t')
                effect = splice_predict.predict(chrom=elems[0],
                                                pos=int(elems[1]),
                                                ref=elems[3],
                                                alt=elems[4].split(','),
                                                genepanel=args.genepanel,
                                                refseqgene=args.refseqgene,
                                                refseq=args.refseq)
                INFO = elems[7]
                new_INFO = []
                if 'splice=' not in INFO:
                    new_INFO = [INFO, splice_predict.print_vcf(effect)]
                else:
                    for info_item in INFO.split(';'):
                        if info_item.startswith('splice='):
                            new_INFO.append('splice=' + splice_predict.print_vcf(effect))
                        else:
                            new_INFO.append(info_item)

                new_line= '\t'.join(elems[:6] + [';'.join(new_INFO)] + elems[8:])
                return new_line

            def write_header(vcf_input):
                vcf_reader = vcf.Reader(vcf_input)
                vcf_info_desc = 'Splice effect. Format: Transcript|Effect|MaxEntScan-wild|MaxEntScan-mut|MaxEntScan-closest|dist'
                vcf_reader.infos['splice'] = VcfInfo(id='splice',
                                                     num=1,
                                                     type='String',
                                                     desc=vcf_info_desc,
                                                     source='AMG-OUS-splice',
                                                     version=__version__)
                vcf.Writer(vcf_output, vcf_reader)
            
            with open(input_filename, 'r') as vcf_input:
                write_header(vcf_input)

                pool = Pool(processes=8)
                for line in vcf_input:
                    if line.startswith('#'):
                        pass
                    else:
                        pool.apply_async(process, args=(line, ), callback=log_result)
                pool.close()
                pool.join()
                flush(result_list)
