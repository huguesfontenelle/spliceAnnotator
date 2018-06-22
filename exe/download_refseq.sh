#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ ! -f ${DIR}/../data/human_g1k_v37.fasta.fai ]; then
    echo "Downloading 'human_g1k_v37.fasta.fai' ..."
    wget -P ${DIR}/../data/ ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz
    wget -P ${DIR}/../data/ ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz.md5
    if [ $(md5 ${DIR}/../data/human_g1k_v37.fasta.fai.gz | cut -d' ' -f4) = $(cat ${DIR}/../data/human_g1k_v37.fasta.fai.gz.md5 | cut -d' ' -f1) ]; then
        gunzip -d ${DIR}/../data/human_g1k_v37.fasta.fai.gz
        rm ${DIR}/../data/human_g1k_v37.fasta.fai.gz.md5
    else
        echo "MD5 for 'human_g1k_v37.fasta.fai.gz' does not match!"
        exit 1
    fi
else
    echo "'human_g1k_v37.fasta.fai' alread downloaded!"
fi

if [ ! -f ${DIR}/../data/human_g1k_v37.fasta ]; then
    echo "Downloading 'human_g1k_v37.fasta' ..."
    wget -P ${DIR}/../data/ ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz
    wget -P ${DIR}/../data/ ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz.md5
    if [ $(md5 ${DIR}/../data/human_g1k_v37.fasta.gz | cut -d' ' -f4) = $(cat ${DIR}/../data/human_g1k_v37.fasta.gz.md5 | cut -d' ' -f1) ]; then
        gunzip -d ${DIR}/../data/human_g1k_v37.fasta.gz
        rm ${DIR}/../data/human_g1k_v37.fasta.md5
    else
        echo "MD5 for 'human_g1k_v37.fasta.gz' does not match!"
        exit 1
    fi
else
    echo "'human_g1k_v37.fasta' alread downloaded!"
fi
