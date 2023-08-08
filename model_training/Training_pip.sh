#!/bin/env bash


train_split_path=$1
home_dir=$2
base_editior=$3
bio_number=$4
#source activate py3.8


python /home/wull01/ACBE_endomodel_reanalysis/BE_Endo_Smart/model_training/Eff_training_main.py ${train_split_path} ${home_dir} ${base_editior} ${bio_number}


if [ $? -eq 0 ]; then
echo "S"
else
echo "${train_split_path} F" >&2
fi
