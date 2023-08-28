# BE_Endo
title: Deep learning models incorporating endogenous factors beyond DNA sequences improve the prediction accuracy of base editing outcomes

##download 
git clone https://github.com/zhaodalv/BE_Endo.git (No further action is required!!)
Or download zip file at https://github.com/zhaodalv/BE_Endo
Once downlaod the zip file the "-master" suffix need to removed by this command: "mv BE_Endo-master BE_Endo", which is a necessary step to run the progam!!!

####set up running envirnoment within conda version 4.14.0
1) create running ENV (BE_Endo)
conda create -n BE_Endo python=3.9.4
2) activating ENV (BE_Endo)
conda activate BE_Endo
3)install package dependency
conda install -c conda-forge joblib=0.17.0
conda install -c conda-forge numpy=1.20.2
conda install -c conda-forge pandas=1.2.4
conda install -c anaconda scipy=1.6.2
conda install -c anaconda seaborn=0.11.2
conda install -c anaconda statsmodels=0.12.2
conda install tensorflow=2.4.1 // conda install -c conda-forge tensorflow=2.4.1
conda install -c bioconda bedtools=2.29.1
conda install -c anaconda scikit-learn=0.24.2
4) add package path to environment
export PYTHONPATH=$PYTHONPATH:{PATH of BE_Endo Package}
example: export PYTHONPATH=$PYTHONPATH:~/BE_Endo

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
                        Package path of BE_Endo
  -B BEDOUT, --BEDOUT BEDOUT
                        Output BED file for intersection
  -O OUTPREFIX, --Outprefix OUTPREFIX
                        Outfile prefix to store effiency and propotion result
  -E ENDOGENOUS, --Endogenous ENDOGENOUS
                        Choose one of endogenous factor name:expression methylation Dnase H3K27ac H3K4me3 PII H3K4me1 CTCF H3K36me3'

Main code "Prediction_main.py" offers two Prediction Mode:
First: Default Prediction Mode. This Mode includes prediction using models: BE_Seq and BE_Endo model.
example:
python Prediction_main.py -T CBE -P /home/wull01 -B $PWD/test1.bed -O $PWD/test/test1_ test_input

output:
1) Predicted BE_Seq eff/proportion table
2) Predicted BE_Endo eff/proportion table


Second: Customizing Prediction Mode. This Mode includes prediction using user provided endogenous factor selected from "expression methylation Dnase H3K27ac H3K4me3 PII H3K4me1 CTCF H3K36me3".
example:
python Prediction_main.py -T CBE -P /home/wull01 -B $PWD/test1.bed -O $PWD/test/CBE_custize_methylation -E methylation test_input

output:
1) Predicted single factor BE_Endo eff/proportion table

Notice: Here Package "BE_Endo" in the director "/home/wull01"

Input file detailed format 
TCCAGGCCAGGCAGTCAAGGAAGAAGCCCTGGGCTCCCAG        chr14:24431547-24431567
#input 40bp Composition
TCCAGGCCAG GCAGTCAAGGAAGAAGCCCT GGG  CTCCCAG        
|10bp up | |    protospacer   ||PAM| | 7bp |      

####Error fix
1)
PermissionError: [Errno 13] Permission denied: 'BD_intersect.sh'
chmod a+x BD_intersect.sh

2)
Traceback (most recent call last):
  File "XXX/BE_Endo/Main/Prediction_main.py", line 19, in <module>
    PIP_F.pip_file(args.Input_file,args.Editor,args.BEDOUT,args.Package,args.Outprefix)
  File "XXX/BE_Endo/Main/File_input_main.py", line 58, in pip_file
    Model_prepare_obj.read_ins("{}/BE_Endo/Main/Intersection_temp/{}".format(workdir,tmpname))
  File "/home/wull01/ACBE_endomodel_reanalysis/BE_Endo/BE_endo_Pack/Preprocess/Preprocessing.py", line 94, in read_ins
    result=pd.read_csv(intersected_path,sep="\t",header=None)
  File "/home/wull01/anaconda3/envs/py3.8/lib/python3.9/site-packages/pandas/io/parsers.py", line 610, in read_csv
    return _read(filepath_or_buffer, kwds)
  File "/home/wull01/anaconda3/envs/py3.8/lib/python3.9/site-packages/pandas/io/parsers.py", line 462, in _read
    parser = TextFileReader(filepath_or_buffer, **kwds)
  File "/home/wull01/anaconda3/envs/py3.8/lib/python3.9/site-packages/pandas/io/parsers.py", line 819, in __init__
    self._engine = self._make_engine(self.engine)
  File "/home/wull01/anaconda3/envs/py3.8/lib/python3.9/site-packages/pandas/io/parsers.py", line 1050, in _make_engine
    return mapping[engine](self.f, **self.options)  # type: ignore[call-arg]
  File "/home/wull01/anaconda3/envs/py3.8/lib/python3.9/site-packages/pandas/io/parsers.py", line 1898, in __init__
    self._reader = parsers.TextReader(self.handles.handle, **kwds)
  File "pandas/_libs/parsers.pyx", line 521, in pandas._libs.parsers.TextReader.__cinit__
pandas.errors.EmptyDataError: No columns to parse from file

Means our endo_factors.bed Does Not put to the right directory
Endo BED file download
In directory endo_factors, There is a "endo_factors.bed" file (endo_factors/endo_factors.bed).
The "endo_factors.bed" is shared via Pan.baidu: "https://pan.baidu.com/s/1wXlLXSgPXN4tdgTqgq8DeQ" with extraction code: "endo"
One downloaded, extracted from zip file and put "endo_factors.bed" in the director "endo_factors".

3)
No such file or directory:XXXXX/BD_intersect.sh 
Change the interpreter of BD_intersect.sh to your own bash interpreter. Such as in MAC, "#!/bin/env bash" is a bad interpreter.Changed to "#!/bin/bash"

4)
Traceback (most recent call last):
  File "/home/wull01/BE_Endo/Main/Prediction_main.py", line 6, in <module>
    import File_input_main  as PIP_F
  File "/home/wull01/BE_Endo/Main/File_input_main.py", line 7, in <module>
    import BE_endo_Pack.Preprocess.Preprocessing as Prepro
ModuleNotFoundError: No module named 'BE_endo_Pack'

export PYTHONPATH=$PYTHONPATH:/home/wull01/BE_Endo
