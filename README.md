# BE_Endo_Smart
title: Deep learning models incorporating endogenous factors beyond DNA sequences improve the prediction accuracy of base editing outcomes

####set up running envirnoment within conda version 4.14.0
1) create running ENV (BE_Endo_Smart)
conda create -n BE_Endo_Smart python=3.9.4
2) activating ENV (BE_Endo_Smart)
conda activate BE_Endo_Smart
3)install package dependency
conda install -c conda-forge joblib=0.17.0
conda install -c conda-forge numpy=1.20.2
conda install -c conda-forge pandas=1.2.4
conda install -c anaconda scipy=1.6.2
conda install -c anaconda seaborn=0.11.2
conda install -c anaconda statsmodels=0.12.2
conda install tensorflow=2.4.1
conda install -c bioconda bedtools=2.29.1
conda install -c anaconda scikit-learn=0.24.2
4) add package path to environment
export PYTHONPATH=$PYTHONPATH:{PATH of BE_Endo_Smart Package}
example: export PYTHONPATH=$PYTHONPATH:~/BE_Endo_Smart

####Usage informations
python Prediction_main.py -h
usage: Prediction_main.py [-h] -T EDITOR -P PACKAGE -B BEDOUT -O OUTPREFIX Input_file
positional arguments:
  Input_file            40bp(10bp upstream + 20bp sgRNA+3bp PAM + 7bp downstream)(Tab seperate)chr:start-end
optional arguments:
  -h, --help            show this help message and exit
  -T EDITOR, --Editor EDITOR
                        BASE editor(ABE/CBE/abe/cbe)
  -P PACKAGE, --Package PACKAGE
                        Package path of BE_Endo_Smart
  -B BEDOUT, --BEDOUT BEDOUT
                        Output BED file for intersection
  -O OUTPREFIX, --Outprefix OUTPREFIX
                        Outfile prefix to store effiency and propotion result

example: 
python Prediction_main.py -T ABE -P /home/wull01 -B $PWD/test1.bed -O $PWD/test/test1_ test_input
python Prediction_main.py -T CBE -P /home/wull01 -B $PWD/test1.bed -O $PWD/test/test1_ test_input
Notice: here Package "BE_Endo_Smart" in the director "/home/wull01"

Input file detailed format 
TCCAGGCCAGGCAGTCAAGGAAGAAGCCCTGGGCTCCCAG        chr14:24431547-24431567
#input 40bp Composition
TCCAGGCCAG GCAGTCAAGGAAGAAGCCCT GGG  CTCCCAG        
|10bp up | |    protospacer   ||PAM| | 7bp |      


output

1) Predicted single factor eff/proportion table
2) Predicted all factor eff/proportion table

Error fix
PermissionError: [Errno 13] Permission denied: 'BD_intersect.sh'
chmod a+x BD_intersect.sh

Endo BED file download
In directory endo_factors, There is a "endo_factors.bed" file (endo_factors/endo_factors.bed).
The "endo_factors.bed" is shared via Pan.baidu: "https://pan.baidu.com/s/1wXlLXSgPXN4tdgTqgq8DeQ" with extraction code: "endo"
One downloaded, extracted from zip file and put "endo_factors.bed" in the director "endo_factors".
