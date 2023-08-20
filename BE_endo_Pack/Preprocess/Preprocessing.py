import numpy as np
import joblib
from subprocess import Popen,PIPE
import re
import pandas as pd
from functools import reduce

class input_process:
    def __init__(self,input_file,editor):
        self.Ipath=input_file
        self.editor=editor
    def process(self):
        seq_40s=[]
        targets=[]
        labels=[]
        positions=[]
        with open(self.Ipath) as f:
            for num,line in enumerate(f):
                seq40,pos=line.strip().split("\t")
                try:
                    self.check_seq(seq40)
                except TypeError:
                    print(seq40,'is error formatted!!')
                    continue
                seq_40s.append(seq40)
                sgRNA=self.sgRNA(seq40)
                targets.append(sgRNA)

                labels.append(self.labels(sgRNA,self.editor))
                positions.append(self.positions(pos,num))

        return seq_40s,targets,labels,positions
    
    def single_line_process(self,seq40,pos):
        try:
            self.check_seq(seq40)
        except TypeError:
            return False
    
    def sgRNA(self,seq):
        return seq[10:30]
    def labels(self,sgRNA,editor):
        label=[]
        intab_ref=""
        outtab_ref=""
        if editor=='CBE':
            intab_ref='ATCGN-'
            outtab_ref='001000'
        if editor=='ABE':
            intab_ref='ATCGN-'
            outtab_ref='100000'
        transtab_ref=str.maketrans(intab_ref,outtab_ref)
        ref_translate=[int(x) for x in str.translate(sgRNA,transtab_ref)]
        target_indexes=[]
        for index,value in enumerate(ref_translate):
            if value==1:
                target_indexes.append(index)
                label.append(True)
            else:
                label.append(False)
        return label
    def positions(self,pos,order):
        sgRNA_pos=re.split(':|-',pos)
        sgRNA_pos.append('L'+str(order))
        return sgRNA_pos
    def check_seq(self,sequence):
        length=len(sequence)
        PAM=sequence[30:33] in ['AGG','CGG','TGG','GGG']
        if not (length==40 and PAM):
            raise TypeError("Sequence Type Error, referring --help for input information")
    def write_out(self,df,out_file):
        df.to_csv(out_file,sep="\t",header=None,index=None)
        return out_file


class Endo_Prepare:
    def __init__(self,Editor,sgRNA_bed):
        self.ET=Editor
        #self.MI=model_info
        self.Iinput=sgRNA_bed
        self.intersection=""
    #def Selecting_model(self):
    #    Model_dict=self.F_to_dict(self.MI)
    #    return Model_dict[self.ET]
    def Prepare_inputs(self):
        input_df=pd.read_csv(self.Iinput,sep="\t",header=None)
        input_df.columns=['sgchrom','sgstart','sgend','line']
        input_df.set_index('line',inplace=True)
        return input_df
    def intersecting(self,main_script,workdir,Tmp_name):
        process1=Popen([main_script,self.Iinput,workdir,Tmp_name],stdout=PIPE,stderr=PIPE)
        process1.communicate()
    def read_ins(self,intersected_path): ##readin intersect results
        result=pd.read_csv(intersected_path,sep="\t",header=None)
        result.columns=['sgchrom','sgstart','sgend','line','chr_T','start_T','end_T','MM_value','factor']
        result.loc[result['chr_T']=='.','MM_value']=0.0
        result['MM_value']=result['MM_value'].astype(float)
        self.intersection=result
    def factor_process(self,index_df): ###index_df is line and 'sgchrom','sgstart','sgend'
        Order=['expression','methylation','Dnase','H3K27ac','H3K4me3','PII','H3K4me1','CTCF','H3K36me3']
        #ABE_order=['H3K27ac']
        factor_list=[index_df]
        input_length=index_df.shape[0]
        #if self.ET=='CBE':
        factor_percentage={}
        for epi in Order:
            out_df=self.epi_process(epi,self.intersection)
            out_df.columns=[epi]
            percentage=len(out_df)/input_length*100
            print("{} is account {}%".format(epi,percentage))
            if epi=='expression' and percentage < 40:
                percentage=0
            factor_percentage[epi]=percentage
            #factor_percentage[epi]=percentage
            
            factor_list.append(out_df)
        ordered_dict={k: v for k, v in sorted(factor_percentage.items(), key=lambda item: item[1],reverse=True)}
        #print("ordered factor is: ", *ordered_dict.keys())
        #else:
        #    for epi in ABE_order:
        #        out_df=self.epi_process(epi,self.intersection)
        #        out_df.columns=[epi]
        #        factor_list.append(out_df)
        #if self.ET=='CBE':
        single_factor=list(ordered_dict.keys())[0]
        merge=reduce(lambda  left,right: pd.merge(left,right,on=['line'],
                                           how='outer'), factor_list)
        #for key,value in self.imputation.items():
        for key in Order:
            merge.loc[pd.isna(merge[key]),key]=0
        #print(merge.columns)
        #print("choose single_factor {} for effiency prediction".format(single_factor))
        #all_factor=merge.copy()
        #else:
        #    merge=reduce(lambda  left,right: pd.merge(left,right,on=['line'],
        #                                    how='outer'), factor_list)
       #     merge.loc[pd.isna(merge['H3K27ac']),'H3K27ac']=self.imputation['H3K27ac']
        return single_factor,merge

    def epi_process(self,factor,df):
        sub_df=df.loc[df['factor']==factor].copy()
        out_df=""
        if factor=='expression':
            out_df=sub_df.groupby('line').max()
        else:
            out_df=sub_df.groupby('line').mean()

        return out_df[['MM_value']]
