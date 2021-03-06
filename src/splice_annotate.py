# -*- coding: utf-8 -*-
'''
SpliceAnnotate

@author: Hugues Fontenelle, 2015
'''

__version__ = '0.5'

import argparse
import vcf
from vcf.parser import _Info as VcfInfo
import splice_predict
from multiprocessing import Pool
import tempfile, shutil
import subprocess
import os


def sort_vcf(fname):
    '''Due to multithreading, splice_annotate does not keep the VCF sorted'''
    outfd, outsock_path = tempfile.mkstemp(suffix='.vcf', prefix='sort', dir=None, text=True)
    with open(outsock_path, 'w') as output_file:
        # sorting VCF
        # grep ^# input.vcf > output.vcf
        subprocess.call(['grep', '^#', fname], stdout=output_file, stderr=subprocess.STDOUT)
        # grep -v ^# input.vcf | sort -V -k1 -k2 >> output.vcf
        p1 = subprocess.Popen(['grep', '-v', '^#', fname], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p2 = subprocess.Popen(['sort', '-V', '-k1', '-k2'], stdin=p1.stdout, stdout=output_file, stderr=subprocess.STDOUT)
        output = p2.communicate()
        shutil.copy(outsock_path, fname)
    os.close(outfd)
    os.remove(outsock_path)


# ============================================================
if __name__ == "__main__":
    '''
    $ python splice/splice_annotate.py --help
    usage: splice_annotate.py [-h] -i INPUT [-o OUTPUT] [-t GENEPANEL]
                              [--refGene REFGENE] [-R REFSEQ]

    Annotate a VCF with splice site effect prediction.

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Input VCF file
      -o OUTPUT, --output OUTPUT
                            Output VCF file for annotation of predictions
                            (default=same as input)
      -t GENEPANEL, --genepanel GENEPANEL
                            Filepath for genepanel
      --refGene REFGENE     Filepath for refGene
      -R REFSEQ, --refseq REFSEQ
                            Filepath for RefSeq (reference FASTA sequence)
      -np NP                Number of processes (default=8)
      --no-denovo           Do NOT make predictions for cryptic denovo splice
                            sites
    Remark: The VCF is sorted with UNIX bash sort, not using the contigs. This is
    fine for most purposes.

    Example usage:
    $ python splice_annotate.py -i data/case_study0.vcf -o data/case_study0_annot.vcf -t HBOC_OUS_medGen_v00_b37.transcripts.csv -R human_g1k_v37_decoy.fasta
    $ python splice_annotate.py -i data/case_study0.vcf -o data/case_study0_annot.vcf -refGene refGene_131119.tab -R human_g1k_v37_decoy.fasta
    '''
    parser = argparse.ArgumentParser(
        description='Annotate a VCF with splice site effect prediction.',
        epilog='Remark: The VCF is sorted with UNIX bash sort, not using the contigs. This is fine for most purposes.')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input VCF file')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Output VCF file for annotation of predictions')
    parser.add_argument('-t', '--genepanel', type=str, required=False,
                        help='Filepath for genepanel')
    parser.add_argument('--refGene', type=str, required=False,
                        help='Filepath for refGene')
    parser.add_argument('-R', '--refseq', type=str, required=True,
                        help='Filepath for RefSeq (reference FASTA sequence)')
    parser.add_argument('-np', type=int, required=False, default=8,
                        help='Number of processes (default=8)')
    denovo_parser = parser.add_mutually_exclusive_group(required=False)
    denovo_parser.add_argument('--denovo', dest='denovo', action='store_true',
        help='Make predictions for cryptic denovo splice sites (default)')
    denovo_parser.add_argument('--no-denovo', dest='denovo', action='store_false',
        help='Do NOT make predictions for cryptic denovo splice sites')
    parser.set_defaults(denovo=True)
    args = parser.parse_args()
    assert args.genepanel or args.refGene

    # Developper notes: everything is written within __main__ because of issues with multiprocessing.Pool:
    # - the callback of apply_async take only one param
    # - the process and callback must be picke-able wich (difficult in Python2 with class methods)
    # - result_list should be a nonlocal, but again this comes in Python3

    with open(args.output, 'w') as vcf_output, open(args.input, 'r') as vcf_input:
        result_list = []

        def flush(result_list):
            for line in result_list:
                vcf_output.write(line + '\n')

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
                                            refseqgene=args.refGene,
                                            refseq=args.refseq,
                                            denovo=args.denovo)
            INFO = elems[7]
            new_INFO = []
            for info_item in INFO.split(';'):
                info_item = info_item.strip()
                if not info_item.startswith('splice=') and info_item not in ['', '.']:
                    new_INFO.append(info_item.strip())
            new_INFO.append(splice_predict.print_vcf(effect))

            new_line = '\t'.join(elems[:7] + [';'.join(new_INFO)] + elems[8:])
            return new_line

        # write_header
        vcf_reader = vcf.Reader(vcf_input)
        vcf_info_desc = 'Splice effect. Format: Transcript|Effect|MaxEntScan-wild|MaxEntScan-mut|MaxEntScan-closest|dist'
        vcf_reader.infos['splice'] = VcfInfo(id='splice',
                                             num=1,
                                             type='String',
                                             desc=vcf_info_desc,
                                             source='spliceAnnotator',
                                             version=__version__)
        vcf.Writer(vcf_output, vcf_reader)

        # distribute annotation
        pool = Pool(processes=args.np)
        for line in vcf_input:
            if line.startswith('#'):
                pass
            else:
                pool.apply_async(process, args=(line, ), callback=log_result)
        pool.close()
        pool.join()
        flush(result_list)

    sort_vcf(args.output)
