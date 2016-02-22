#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

FILES="${DIR}/../data/case_study0.vcf"

REFSEQ="${DIR}/../../../../../vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta"
REFSEQGENE="${DIR}/../../../../../vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab"

export PYTHONPATH="${DIR}/.."

for f in $FILES
do
	python ${DIR}/../splice/splice_annotate.py -i $f  --refseqgene $REFSEQGENE --refseq $REFSEQ
done
