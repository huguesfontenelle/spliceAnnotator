# Splice site annotation and prediction

## Annotation

`splice_score.py`: annotates a VCF file.

## Scoring

`splice_score.py`: scores a variant given in JSON input
Scoring algorithms are:  
  * `SSFL.py`: Product Weight Matrix-based algo, implemented in Python
  * `MaxEntScan`: Bayesian Decision, Python wrapper for Perl code 
  * `GeneSplicer`: Decision Tree, Python wrapper for compiled code
  * `NNSplice`: Neural Network, Python code for web scrapping
  * `HSF`: Human Splice Finder
  
## Prediction

* `splice_predict`: predicts the effect of a variant given in JSON format.

