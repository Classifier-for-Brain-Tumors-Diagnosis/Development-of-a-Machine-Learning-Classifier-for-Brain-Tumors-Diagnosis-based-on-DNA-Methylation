# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 09:13:50 2020

@author: Administrator
"""
import numpy as np 
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
label1=pd.read_csv('/home/chenyx/dnameth/diff_feature/gse90496_label.csv')
label2=pd.read_csv('/home/chenyx/dnameth/diff_feature/gse109379_label.csv')
label3=pd.read_csv('/home/chenyx/dnameth/valdataset/valdataset_idat_all/SampleSheet.csv')
beta_dir1='/home/chenyx/dnameth/diff_feature/gse90496_rf_importance_order.csv'
beta_dir2='/home/chenyx/dnameth/gse109379/beta/beta_testdataset_sub_v2.csv'
beta_dir3='/home/chenyx/dnameth/gseval/beta_valdataset_sub_v2.csv'
beta1=pd.read_csv(beta_dir1)
beta2=pd.read_csv(beta_dir2,index_col=0)
beta3=pd.read_csv(beta_dir3,index_col=0)
m=767
beta1_selected=beta1.iloc[:,0:m]
beta2_selected=beta2.iloc[:,0:m]
beta3_selected=beta3.iloc[:,0:m]
mlp_clf= MLPClassifier(activation='relu',solver='sgd',learning_rate='adaptive',hidden_layer_sizes=(1320,),batch_size=2,learning_rate_init=0.001,power_t=0.1,max_iter=1000,random_state=42,early_stopping=True)
mlp_clf.fit(beta1_selected,label1.iloc[:,1])
with open('/home/chenyx/dnameth/train_test_val_v2/predict_splited.txt','a') as f:
    train_pred=mlp_clf.predict(beta1_selected)
    train_prob=mlp_clf.predict_proba(beta1_selected)
    train_prob=pd.DataFrame(train_prob)
    train_prob.columns=mlp_clf.classes_
    train_con=pd.DataFrame(confusion_matrix(label1.iloc[:,1],train_pred,labels=mlp_clf.classes_))
    train_con.columns=mlp_clf.classes_
    train_con.index=mlp_clf.classes_
    train_con.to_csv('/home/chenyx/dnameth/train_test_val_v2/train_confustion.csv')
    train_prob.to_csv("/home/chenyx/dnameth/train_test_val_v2/train_prob_nocorrected.csv",index=False)
    f.write("train_acc={}\n".format(accuracy_score(label1.iloc[:,1],train_pred)))
    pd.concat([pd.DataFrame(label1.iloc[:,1]),pd.DataFrame(train_pred)],axis=1).to_csv('/home/chenyx/dnameth/train_test_val_v2/train_true_pred.csv')
    test_pred=mlp_clf.predict(beta2_selected)
    test_prob=mlp_clf.predict_proba(beta2_selected)
    test_prob=pd.DataFrame(test_prob)
    test_prob.columns=mlp_clf.classes_
    test_con=pd.DataFrame(confusion_matrix(label2.iloc[:,1],test_pred,labels=mlp_clf.classes_))
    test_con.columns=mlp_clf.classes_
    test_con.index=mlp_clf.classes_
    test_con.to_csv('/home/chenyx/dnameth/train_test_val_v2/test_confustion.csv')
    test_prob.to_csv("/home/chenyx/dnameth/train_test_val_v2/test_prob_nocorrected.csv",index=False)
    f.write("test_acc={}\n".format(accuracy_score(label2.iloc[:,1],test_pred)))
    pd.concat([pd.DataFrame(label2.iloc[:,1]),pd.DataFrame(test_pred)],axis=1).to_csv('/home/chenyx/dnameth/train_test_val_v2/test_true_pred.csv')
    val_pred=mlp_clf.predict(beta3_selected)
    val_prob=mlp_clf.predict_proba(beta3_selected)
    val_prob=pd.DataFrame(val_prob)
    val_prob.columns=mlp_clf.classes_
    val_con=pd.DataFrame(confusion_matrix(label3.iloc[:,4],val_pred,labels=mlp_clf.classes_))
    val_con.columns=mlp_clf.classes_
    val_con.index=mlp_clf.classes_
    val_con.to_csv('/home/chenyx/dnameth/train_test_val_v2/val_confustion.csv')
    val_prob.to_csv("/home/chenyx/dnameth/train_test_val_v2/val_prob_nocorrected.csv",index=False)
    f.write("val_acc={}\n".format(accuracy_score(label3.iloc[:,4],val_pred)))
    pd.concat([pd.DataFrame(label3.iloc[:,4]),pd.DataFrame(val_pred)],axis=1).to_csv('/home/chenyx/dnameth/train_test_val_v2/val_true_pred.csv')
    gse_sort=np.unique(label3.iloc[:,3])
    for i in np.arange(len(gse_sort)):
        temp_index=np.where(label3.iloc[:,3]==gse_sort[i])[0]
        temp_val=beta3_selected.iloc[temp_index,:]
        val_true=pd.DataFrame(label3.iloc[temp_index,4])
        val_true.index=np.arange(len(val_true))
        val_pred=mlp_clf.predict(beta3_selected.iloc[temp_index,:])
        t_p=pd.concat([pd.DataFrame(val_true),pd.DataFrame(val_pred)],axis=1)
        t_p.to_csv("/home/chenyx/dnameth/train_test_val_v2/{}_t_p_nocorrected.csv".format(gse_sort[i]))
        f.write("{}_acc={}\n".format(gse_sort[i],accuracy_score(val_true,val_pred)))





