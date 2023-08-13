#!/bin/env python
import sys
import pandas as pd
import os
import numpy as np
#sys.path.append("..")
import BE_endo_Pack.Preprocess.Preprocessing as Prepro
import BE_endo_Pack.Effiency.Eff as Effiency
import BE_endo_Pack.Proportion.proportion as Proportion


def out_eff_table(eff_prediction,input_df,seq_40s,target):
    eff_table=pd.DataFrame(eff_prediction)
    eff_table.columns=['base1', 'base2', 'base3', 'base4', 'base5', 'base6', 'base7', 'base8',
   'base9', 'base10', 'base11', 'base12', 'base13', 'base14', 'base15',
   'base16', 'base17', 'base18', 'base19', 'base20']
    eff_table['seq']=seq_40s
    eff_table['target']=target
    eff_table['chr']=input_df['sgchrom'].values
    eff_table['start']=input_df['sgstart'].values
    eff_table['end']=input_df['sgend'].values
    return eff_table

def pip_file(file_in,editor,sgRNA_bed_path,workdir,out_file,prop_cut=0.01):
   model_path=""
   pre_path=""
   Editor=editor.upper()
   #tmpname=basename(out_file)
   #model_expression="../model_data_endo/{}_eff/{}_expression_endo.h5".format(Editor,Editor)
   #model_PII="../model_data_endo/{}_eff/{}_PII_endo.h5".format(Editor,Editor)
   #model_CTCF="../model_data_endo/{}_eff/{}_CTCF_endo.h5".format(Editor,Editor)
   #model_Dnase="../model_data_endo/{}_eff/{}_Dnase_endo.h5".format(Editor,Editor)
   #model_H3K4me1="../model_data_endo/{}_eff/{}_H3K4me1_endo.h5".format(Editor,Editor)
   #model_H3K4me3="../model_data_endo/{}_eff/{}_H3K4me3_endo.h5".format(Editor,Editor)
   #model_H3K27ac="../model_data_endo/{}_eff/{}_H3K27ac_endo.h5".format(Editor,Editor)
   #model_H3K36me3="../model_data_endo/{}_eff/{}_H3K36me3_endo.h5".format(Editor,Editor)
   #model_methylation="../model_data_endo/{}_eff/{}_methylation_endo.h5".format(Editor,Editor)
   #model_all="../model_data_endo/{}_eff/{}_all_endo.h5".format(Editor,Editor)
   pre_path="{}/BE_Endo/model_data_endo/{}_pro/{}_proportion_pre.npy".format(workdir,Editor,Editor)



  # if Editor=='CBE':
      
  #    model_path=workdir+"/BE_endo_smart/model_data/CBE_CBE_all_0_endo.h5"
  #    pre_path=workdir+'/BE_endo_smart/model_data/CBE_proportion_pre.npy'
  # else:
  #    model_path=workdir+"/BE_endo_smart/model_data/ABE_ABE_H3K27ac_3_endo.h5" 
  #    pre_path=workdir+'/BE_endo_smart/model_data/ABE_proportion_pre.npy'
   Input_process_obj=Prepro.input_process(file_in,Editor) ##create preprocess obj 
   seq_40s,targets,labels,positions=Input_process_obj.process() ## obtain seq_40s ,targets,labels,positions
   position_df=pd.DataFrame(positions) ##sgRNA position write to local for intersection
   out_file_path=Input_process_obj.write_out(position_df,sgRNA_bed_path) ###write out sgRNA.bed
   tmpname=os.path.basename(out_file_path)
   Model_prepare_obj=Prepro.Endo_Prepare(Editor,out_file_path) ##create Model_prepare_obj
   Model_prepare_obj.intersecting("{}/BE_Endo/Main/BD_intersect.sh".format(workdir),workdir,tmpname)
   input_df=Model_prepare_obj.Prepare_inputs()
   Model_prepare_obj.read_ins("{}/BE_Endo/Main/Intersection_temp/{}".format(workdir,tmpname))
   single_factor,merge_result=Model_prepare_obj.factor_process(input_df)
   model_path="{}/BE_Endo/model_data_endo/{}_eff/{}_{}_endo.h5".format(workdir,Editor,Editor,single_factor)
   model_all_path="{}/BE_Endo/model_data_endo/{}_eff/{}_all_endo.h5".format(workdir,Editor,Editor)
   ###prediction
   Seq_obj=Effiency.Sequence_prepare(seq_40s)
   X_array=Seq_obj.one_hot_encoding()
   print(model_all_path," is predicting......")
   eff_prediction_all=Effiency.pip_endo(model_all_path,seq_40s,merge_result.iloc[:,3:].values,labels,Editor)
   print(model_path," is predicting......")
   eff_prediction_single=Effiency.pip_endo(model_path,seq_40s,merge_result[[single_factor]].values,labels,Editor)
   eff_table_all=out_eff_table(eff_prediction_all,input_df,seq_40s,targets)
   eff_table_single=out_eff_table(eff_prediction_single,input_df,seq_40s,targets)
   
   #eff_table.columns=['base1', 'base2', 'base3', 'base4', 'base5', 'base6', 'base7', 'base8',
   #'base9', 'base10', 'base11', 'base12', 'base13', 'base14', 'base15',
   #'base16', 'base17', 'base18', 'base19', 'base20']
   #eff_table['seq']=seq_40s
   #eff_table['chr']=input_df['sgchrom'].values
   #eff_table['start']=input_df['sgstart'].values
   #eff_table['end']=input_df['sgend'].values
   eff_table_all[['chr','start','end','seq','target','base1', 'base2', 'base3', 'base4', 'base5', 'base6', 'base7', 'base8',
   'base9', 'base10', 'base11', 'base12', 'base13', 'base14', 'base15',
   'base16', 'base17', 'base18', 'base19', 'base20']].to_csv(out_file+".{}.all.eff".format(Editor,),sep="\t",index=None)
   

   eff_table_single[['chr','start','end','seq','target','base1', 'base2', 'base3', 'base4', 'base5', 'base6', 'base7', 'base8',
   'base9', 'base10', 'base11', 'base12', 'base13', 'base14', 'base15',
   'base16', 'base17', 'base18', 'base19', 'base20']].to_csv(out_file+".{}.{}.eff".format(Editor,single_factor),sep="\t",index=None)
   proportion_pre=np.load(pre_path)
   Proportion_prediction_obj=Proportion.Proportion_prediction(Editor,proportion_pre)
   Proportion_df_all=Proportion_prediction_obj.prediction(seq_40s,targets,labels,eff_prediction_all)
   Proportion_df_all['wildtype_indicator']=Proportion_df_all['outcomes']==Proportion_df_all['target']
   Proportion_df_single=Proportion_prediction_obj.prediction(seq_40s,targets,labels,eff_prediction_single)
   Proportion_df_single['wildtype_indicator']=Proportion_df_single['outcomes']==Proportion_df_single['target']

   Proportion_df_all_clean=Proportion_df_all.loc[Proportion_df_all['Predicted_eff']>prop_cut].copy()
   Proportion_df_single_clean=Proportion_df_single.loc[Proportion_df_single['Predicted_eff']>prop_cut].copy()

   Proportion_df_all_clean.to_csv(out_file+".{}.all.proportion".format(Editor),sep="\t",index=None)
   Proportion_df_single_clean.to_csv(out_file+".{}.{}.proportion".format(Editor,single_factor),sep="\t",index=None)



