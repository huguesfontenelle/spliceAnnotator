#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

REFSEQ="${DIR}/../../../../../vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta"
REFSEQGENE="${DIR}/../../../../../vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab"

export PYTHONPATH="${DIR}/.."

python ${DIR}/../splice/splice_annotate.py -i $1 --refseqgene $REFSEQGENE --refseq $REFSEQ
