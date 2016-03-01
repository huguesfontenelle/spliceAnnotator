# Splice site annotation and prediction

## Background

There has been countless approaches at predicting the deleteriousness of a variant on splicing. MaxEntScan [1] for example can preform this prediction when the score of the wild type is compared to the variant type. As shown by Jian in 2014 [2] and Fontenelle in 2015 [3], this algorithm is not only the best, but also better in isolation than combined with others.

[1] G. Yeo and C. B. Burge, “Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals.,” J. Comput. Biol., vol. 11, no. 2–3, pp. 377–94, Jan. 2004.
[2] X. Jian, E. Boerwinkle, and X. Liu, “In silico prediction of splice-altering single nucleotide variants in the human genome,” Nucleic Acids Res., vol. 42, no. 22, pp. 13534–13544, Nov. 2014.
[3] H. Fontenelle, M. C. Eike, T. Håndstad, E. Skorve, S. T. Koksrud Seljebotn, T. Grünfeld and D. E. Undlien, "Best practice for use of splice variant prediction tools.," 4th course in Next Generation Sequencing, Bertinoro, Italy, 12-16 May 2015 [PDF] (awarded Best Poster)

This library automates the annotation of variants with MaxEntScan [1]. It predicts the effect by calculating the ratio between the wild score and variant type score. The site is predicted to be lost when that ratio is below -21.4% [2].

## Annotation

### Effect on authentic splice site
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
3   30686236    id1 CA  AT  .   .   splice=NM_003242.5|predicted_lost|10.82|0.0
3    30686405    id2    ATG    TTC  .    .   splice=NM_003242.5|predicted_lost|9.26|5.42
3   30686258    id3 C   T   .   .   splice=NM_003242.5|no_effect|10.82|10.82
```

`splice=` is the INFO keyword for this annotation, it is followed by the transcript for which this prediction is valid, the predicted effect, the MaxEntScan score for the wild type and the score for the variant type. The fields are separated by `|` (pipes).
`predicted_lost` means the the variant score dropped significantly with respect to the wild type score.
`predicted_conserved` means the score did not drop significantly, or not at all
`no_effect` means the nucleotide sequence of the authentic splice site was not affected by the variant.
From Burge's paper for the MaxEntScan implementation, we know that "5' donors" are 9-mer with a consensus GT after the junction; and that "3' acceptors" are 23-mer with a consensus AG before the junction. Sites that do not conform to the canonical AG/GT are not taken into consideration (there exists no known algorithm that can).
In addition, a variant may be annotated with `not_in_transcript` if it falls outside a transcript; or with `NA` if something has failed.

### De novo: putative new splice Sites

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2   162060116    id4 T  G  .   .  splice=NM_004180.2|no_effect|9.79|9.79&NM_004180.2|de_novo|0.0|4.7|9.79|10
2    162041853    id5    T    G .    .   splice=NM_004180.2|no_effect|6.64|6.64&NM_004180.2|de_novo|0.0|10.86|6.64|5581
```

Several effects are joined with the `&` character. Here the variant are scored for their effect on the authentic splice site, but in addition "de novo" sites were found. Following the `de_novo` keyword is the MaxEntScan score of the wild sequence in that position, the score of the de novo site in that position, the score of the closest authentic site, and the distance (in nucleotides) to the closest authentic site.
The rules for marking a de novo are:
* If MES score for REF is 0, report "de novo" if score for ALT is ≥4
* If MES score for REF is >0 , report "de novo" if score for ALT is ≥ 25% increased compared to REF score
The score of the authentic closest site is not taken into consideration in those rules, but it is usually considered important [4] so we added it here.

[4] C. Houdayer, “In Silico Prediction of Splice-Affecting Nucleotide Variants,” in In Silico Tools for Gene Discovery, vol. 760, B. Yu and M. Hinchcliffe, Eds. Totowa, NJ: Humana Press, 2011, pp. 269–281.

Notes:
* De novo searches for Donors if the closest authentic site is a Donor. Conversely it searches for acceptors if the closest is an Acceptor.
* Only AG/GT are searched since the score of any other sequence is meaningless.
* Multi-allelic sites are scored/predicted as well, and separated by `,` (comma)

## Command-line tools

```
$ python splice_annotate.py -h
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
```
