import sys
import scipy.stats as sst
import pandas as pd
import numpy as np
import pickle
import re
import collections
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder
import tensorflow as tf
from tensorflow import keras
import pandas as pd
import os

from sklearn.metrics import r2_score

import statsmodels.api as sm
import seaborn as sns




def Endo_model():
    inputs=keras.Input(shape=(1,40,4),name="main_input")
    input_bio=keras.Input(shape=(1,40,1),name='bio_input')
    seq_bio_concat=tf.concat([inputs,input_bio],-1)
    cov2_1x2=keras.layers.Conv2D(10,(1,2),padding='same',name='cov2_1X2',activation='relu',kernel_regularizer=keras.regularizers.l2(0.0002))(seq_bio_concat)
    cov2_1x3=keras.layers.Conv2D(10,(1,3),padding='same',name='cov2_1X3',activation='relu',kernel_regularizer=keras.regularizers.l2(0.0002))(seq_bio_concat)
    cov2_1x4=keras.layers.Conv2D(10,(1,4),padding='same',name='cov2_1X4',activation='relu',kernel_regularizer=keras.regularizers.l2(0.0002))(seq_bio_concat)
    cov2_1x5=keras.layers.Conv2D(10,(1,5),padding='same',name='cov2_1X5',activation='relu',kernel_regularizer=keras.regularizers.l2(0.0002))(seq_bio_concat)
    #concate_layer=keras.layers.Concatenate(name='concate')([cov2_1x2,cov2_1x3,cov2_1x4,cov2_1x5])
    pool1X2=keras.layers.MaxPool2D((1,2),name="max_pool1")(cov2_1x2)
    pool1X3=keras.layers.MaxPool2D((1,2),name="max_pool2")(cov2_1x3)
    pool1X4=keras.layers.MaxPool2D((1,2),name="max_pool3")(cov2_1x4)
    pool1X5=keras.layers.MaxPool2D((1,2),name="max_pool4")(cov2_1x5)
    concate_layer=keras.layers.Concatenate(name='concate')([pool1X2,pool1X3,pool1X4,pool1X5])
    flatten_layer=tf.keras.layers.Flatten()(concate_layer)
    out_put=tf.keras.layers.Dense(20,activation='sigmoid')(flatten_layer)
    endo_model=keras.models.Model([inputs,input_bio],[out_put])
    return endo_model



def train_home(path,factor):
    directory=path+'/'+factor
    try:
        os.mkdir(directory)
    except:
        print('Directory making fail')
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

def one_hot_array(train,one_hot_obj):
    out_one_hot=[]
    out_index=[]
    for key in train:
        seq=key.split('_')[0]
        out_one_hot.append(one_hot_obj.enconding(seq))
        seq_list=list(seq[10:30])
        out_index.append(seq_list)
    X=np.stack(out_one_hot)
    X_out=X.reshape(len(train),1,40,4)
    df=pd.DataFrame(out_index)
    #row,col=df.shape
    #index=df.values.reshape(row*col,1)
    return X_out,df




class Model_Application: ##input is stacked X,Y
    def __init__(self,mpath,X,Y,editor):
        self.path=mpath
        self.X=X
        self.Y=Y
        self.index=""
        self.indexW=""
        self.editor=editor
        self.model=""
        self.prediction=""
    def get_index(self,index_df): ##index_raw is index_df
        row,col=index_df.shape
        index_raw=index_df.values.reshape(row*col,1)
        if self.editor=='ABE':
            self.index=index_raw=='A'
        else:
            self.index=index_raw=='C'
    def get_window_index(self,index_df):
        row,col=index_df.shape
        if self.editor=='ABE':
            cols=[4,5,6,7]
            index_df.loc[:,~index_df.columns.isin(cols)]='N'
            index_raw=index_df.values.reshape(row*col,1)
            self.indexW=index_raw=='A'
        else:
            cols=[2,3,4,5,6,7,8]
            index_df.loc[:,~index_df.columns.isin(cols)]='N'
            index_raw=index_df.values.reshape(row*col,1)
            self.indexW=index_raw=='C'

    def load_model(self):
        self.model=tf.keras.models.load_model(self.path)
    def predicting(self):
        self.prediction=self.model.predict(self.X)
    def R_valueW(self):
        row,col=self.prediction.shape
        prediction_reshape=self.prediction.reshape(row*col,1)
        real_Y_reshape=self.Y.reshape(row*col,1)

        target_real_Y=real_Y_reshape[self.indexW]
        prediction_Y=prediction_reshape[self.indexW]
        prediction_Y2=sm.add_constant(prediction_Y)

        est = sm.OLS(target_real_Y,prediction_Y2).fit()

        slope, intercept, r_value, p_value, std_err=sst.linregress(prediction_Y,target_real_Y)
        r2=r2_score(target_real_Y,prediction_Y)
        return [slope, intercept,r2,r_value, p_value, std_err,est.rsquared,est.rsquared_adj,est.pvalues]

    def R_value(self):
        row,col=self.prediction.shape
        prediction_reshape=self.prediction.reshape(row*col,1)
        real_Y_reshape=self.Y.reshape(row*col,1)

        target_real_Y=real_Y_reshape[self.index]
        prediction_Y=prediction_reshape[self.index]
        prediction_Y2=sm.add_constant(prediction_Y)

        est = sm.OLS(target_real_Y,prediction_Y2).fit()

        slope, intercept, r_value, p_value, std_err=sst.linregress(prediction_Y,target_real_Y)
        r2=r2_score(target_real_Y,prediction_Y)
        return [slope, intercept,r2,r_value, p_value, std_err,est.rsquared,est.rsquared_adj,est.pvalues]



def pip_endo(mpath,X_raw,Y,bio,one_hot_obj,editor):
    X,raw_index_df=one_hot_array(X_raw,one_hot_obj)
    repeats=np.repeat(bio,40)
    bio_train_reshape=repeats.reshape(len(X),1,40,1)
    Model=Model_Application(mpath,[X,bio_train_reshape],Y,editor)
    Model.load_model()
    Model.get_index(raw_index_df)
    Model.get_window_index(raw_index_df)
    Model.predicting()
    model_result=Model.R_value()
    model_resultW=Model.R_valueW()
    del(Model)
    return model_result,model_resultW


def endo_model_process(X,Y,bio,out_dir,factor,editor):
    home_path=out_dir+'/'+editor+'_'+factor
    if not os.path.isdir(home_path):
        train_home(out_dir,editor+'_'+factor)
    saved_model=out_dir+'/'+editor+'_'+factor+"/"+editor+'_'+factor+'_endo.h5'
    model=Endo_model()
    repeats=np.repeat(bio,40)
    bio_train_reshape=repeats.reshape(len(X),1,40,1)
    model.compile(optimizer=tf.keras.optimizers.Adam(lr=0.001),loss='mse',metrics=['mse'])
    checkpoint_cb=tf.keras.callbacks.ModelCheckpoint(saved_model,save_best_only=True)
    early_stopping_cb=tf.keras.callbacks.EarlyStopping(patience=10,restore_best_weights=True)
    history_convae = model.fit([X,bio_train_reshape],Y,epochs=1000,batch_size=50,validation_split=0.2,callbacks=[checkpoint_cb,early_stopping_cb])
    return saved_model


def Training_pip_endo_all(pkl,home_dir,one_hot_obj,factor,editor,bio_number):
    statistics=open(home_dir+"/"+factor+".statistic",'w+')
    for num,item in enumerate(pkl):
        ID=factor+'_'+str(num)
        X_train,X_train_bio,Y_train,X_test,X_test_bio,Y_test=item
        X_1,train_rawindex=one_hot_array(X_train,one_hot_obj)
        X_2,test_rawindex=one_hot_array(X_test,one_hot_obj)
        X_train_bio[X_train_bio<0]=0
        X_test_bio[X_test_bio<0]=0
        Y_1=np.stack(Y_train)
        endo_path=endo_model_process(X_1,Y_1,X_train_bio,home_dir,ID,editor)
        result,resultw=pip_endo(endo_path,X_train,Y_train,X_train_bio,one_hot_obj,editor)
        print(factor,ID,'endo_train',*result,*resultw,sep="\t",file=statistics)
        result,resultw=pip_endo(endo_path,X_test,Y_test,X_test_bio,one_hot_obj,editor)
        print(factor,ID,'endo_test',*result,*resultw,sep="\t",file=statistics)

    statistics.close()


if __name__=='__main__':
    ##inptut is train_split_data home_dir base_editor_name:CBE/ABE
    train_split_path=sys.argv[1]
    home_dir=sys.argv[2]
    base_editior=sys.argv[3]
    bio_number=sys.argv[4]
    data=""
    one_hot_obj=Seq_encoder('ATCG')
    factor=os.path.basename(train_split_path).split('.')[0]
    with open(train_split_path,'rb') as f:
        data=pickle.load(f)
    Training_pip_endo_all(data,home_dir,one_hot_obj,factor,base_editior,int(bio_number))



