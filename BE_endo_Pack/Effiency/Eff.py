import tensorflow as tf
from sklearn.preprocessing import OneHotEncoder
import numpy as np
import pandas as pd

class Seq_encoder:
    def __init__(self,code):
        self.code_str=code
    def get_array(self,seq):
        code_arr = np.array(list(seq))
        return code_arr.reshape(-1,1)
    def enconding(self,string):
        onehot=OneHotEncoder()
        code_arr=self.get_array(self.code_str)
        string_arr=self.get_array(string)
        onehot.fit(code_arr)
        out_onehot=onehot.transform(string_arr)
        return out_onehot.toarray()

class Sequence_prepare:
    def __init__(self,seq_list):
        self.one_hot_obj=Seq_encoder('ATCG')
        self.seqs=seq_list
    def one_hot_encoding(self):
        out_hot_list=[]
        for seq in self.seqs:
            out_hot_list.append(self.one_hot_obj.enconding(seq))
        X=np.stack(out_hot_list)
        X_array=X.reshape(len(self.seqs),1,40,4)
        return X_array


def transform_all(arrary):
    df=pd.DataFrame(arrary)
    df2=df.where(df>0).fillna(0)
    tensors=[]
    for col in df2.columns:
        tensors.append(np.repeat(df2[col].values,40).reshape(len(arrary),1,40,1))
    all_tensors=tf.concat(tensors,-1)
    return all_tensors       

def transform_single(bio):
    repeats=np.repeat(bio,40)
    bio_train_reshape=repeats.reshape(len(bio),1,40,1)
    return bio_train_reshape

class Model_Application: ##input is stacked X,Y
    def __init__(self,mpath,X):
        self.mpath=mpath
        self.X=X
        self.prediction=""
    def load_model(self):
        self.model=tf.keras.models.load_model(self.mpath)
    def predicting(self):
        self.prediction=self.model.predict(self.X)



def pip_endo(mpath,seq_40s,bio,labels,editor):
    Seq_obj=Sequence_prepare(seq_40s)
    X=Seq_obj.one_hot_encoding()
    bio_train_reshape=""
    #if editor=='CBE':
    bio_train_reshape=transform_all(bio)
    #else:
    #     bio_train_reshape=transform_single(bio)
    #print(bio_train_reshape.shape)
    Model=Model_Application(mpath,[X,bio_train_reshape])
    Model.load_model()
    Model.predicting()
    out_eff=[]
    for label,eff in zip(labels,Model.prediction):
        out_eff.append([y if x else 0 for x,y in zip(label,eff) ])
    del(Model) 
    return out_eff


def pip_Inte(mpath,seq_40s,labels,editor):
    Seq_obj=Sequence_prepare(seq_40s)
    X=Seq_obj.one_hot_encoding()
    bio_train_reshape=""
    if editor=='CBE':
         bio_train_reshape=transform_all(bio)
    else:
         bio_train_reshape=transform_single(bio)

    Model=Model_Application(mpath,[X,bio_train_reshape])
    Model.load_model()
    Model.predicting()
    out_eff=[]
    for label,eff in zip(labels,Model.prediction):
        out_eff.append([y if x else 0 for x,y in zip(label,eff) ])
    del(Model)
    return out_eff


