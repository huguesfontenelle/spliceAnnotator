#!/bin/bash

FILES="../data/case_study0.vcf
../data/case_study1.vcf
../data/case_study2.vcf
../data/case_study3.vcf
../data/de_novo.vcf"

REFSEQ="../../../../vcpipe-bundle/genomic/gatkBundle_2.5/human_g1k_v37_decoy.fasta"
REFSEQGENE="../../../../vcpipe-bundle/funcAnnot/refseq/refGene_131119.tab"

for f in $FILES
do
	# annotate the VCF file using Jian14nar (default) method
	python ../splice/splice_annotate.py -i $f --refseqgene $REFSEQGENE --refseq $REFSEQ
	# Optional: create JSON output for all methods:
	python ../splice/splice_score.py -i $f --all --refseqgene $REFSEQGENE --refseq $REFSEQ
	python ../splice/splice_predict.py -i ${f%.*}.json --method Jian14nar --refseqgene $REFSEQGENE --refseq $REFSEQ
	python ../splice/splice_predict.py -i ${f%.*}.json --method Houdayer --refseqgene $REFSEQGENE --refseq $REFSEQ
	python ../splice/splice_predict.py -i ${f%.*}.json --method AMG-diag --refseqgene $REFSEQGENE --refseq $REFSEQ
	python ../splice/splice_predict.py -i ${f%.*}.json --method AMG-kreftgenetikk --refseqgene $REFSEQGENE --refseq $REFSEQ
done

