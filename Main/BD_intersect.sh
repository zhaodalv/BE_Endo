#!/bin/env bash

Input=$1
Wdir=$2
TNAME=$3
bedtools intersect -loj -wa -wb -a ${Input} -b ${Wdir}/BE_Endo_Smart/endo_factors/endo_factors.bed >${Wdir}/BE_Endo_Smart/Main/Intersection_temp/${TNAME}


if [ $? -eq 0 ];then
echo "bedtools intersecting Success!"

else
echo "bedtools interscting Failure!" >&2

fi

 
