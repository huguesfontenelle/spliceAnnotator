#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

REFSEQ="${DIR}/../data/human_g1k_v37_decoy.fasta"
REFGENE="${DIR}/../data/refGene_180611.tsv"

export PYTHONPATH="${DIR}/../src/:${DIR}/../thirdparty/:${PYTHONPATH}"

python ${DIR}/../src/splice_annotate.py -i $1 --refGene $REFGENE --refseq $REFSEQ -o $2
