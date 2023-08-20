#!/bin/env python

import sys

#sys.path.append("/home/wull01/ACBE_endomodel_reanalysis/BE_Endo")
import File_input_main  as PIP_F
import argparse


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    #ap.add_argument("-M","--ModelType",help="Select model type (Sequence/Combination/Expression/RNAPoll/CTCF/H3K4me3/Dnase/ATAC/Methylation/H3K4me1/H3K27ac/H3K36me3)",type=str,required=True)
    ap.add_argument("-T","--Editor",help="BASE editor(ABE/CBE/abe/cbe)",type=str,required=True)
    ap.add_argument("-P","--Package",help="Package path of BE_Endo",type=str,required=True)
    ap.add_argument("-B","--BEDOUT",help="Output BED file for intersection",type=str,required=True)
    ap.add_argument("-O","--Outprefix",help="Outfile prefix to store effiency and propotion result",type=str,required=True)
    ap.add_argument("-E","--Endogenous",help="Choose one of endogenous factor name:expression methylation Dnase H3K27ac H3K4me3 PII H3K4me1 CTCF H3K36me3'",type=str,default="None")
    ap.add_argument("Input_file",help="40bp(10bp upstream + 20bp sgRNA+3bp PAM + 7bp downstream)(Tab seperate)chr:start-end",type=str)
    args=ap.parse_args()
    if args.Endogenous != 'None' and args.Endogenous not in ['expression','methylation','Dnase','H3K27ac','H3K4me3','PII','H3K4me1','CTCF','H3K36me3']:
        print("Endogenous factor must be selected from: ", *['expression','methylation','Dnase','H3K27ac','H3K4me3','PII','H3K4me1','CTCF','H3K36me3'],sep="\t" )
        print("Prediction Failure!")
        exit (1)
    if args.Endogenous in ['expression','methylation','Dnase','H3K27ac','H3K4me3','PII','H3K4me1','CTCF','H3K36me3']:
        print ("Customizing Prediction using {} model".format(args.Endogenous))
        PIP_F.pip_customization(args.Input_file,args.Editor,args.BEDOUT,args.Package,args.Outprefix,args.Endogenous)
    else:
        print ("Default Model Prediction: Seq model and Endo model")
        PIP_F.pip_file(args.Input_file,args.Editor,args.BEDOUT,args.Package,args.Outprefix)
    
    
    

    
