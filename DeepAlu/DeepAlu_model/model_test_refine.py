# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import tensorflow as tf
import pysam
import sys
import getopt
import random
import math
import os.path
import numpy as np
import datetime, os
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import callbacks
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.io import TFRecordWriter
from deepalu_base_parall import *
os.environ["CUDA_VISIBLE_DEVICES"]="-1"


_TRAIN = 'train.tfrecord'
_EVAL = 'eval.tfrecord'
_TEST = 'test.tfrecord'        
hparams=PileupImageOptions(3)
try:
    opts, args = getopt.getopt(sys.argv[1:],"i:s:t:o:q:r:",["input_gt=","split_softclipped_sort=","threshold=","output=","quick_model=","reference="])
    print(opts)
except getopt.GetoptError:
  print('python model_test.py  -i <input_gt> -o <output_predict_file> -s <split_softclipped_sort> -ref <reference> -q <quick model>')
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
    print('python model_test.py  -i <input_gt> -o <output_predict_file> -s <split_softclipped_sort> -ref <reference>')
    sys.exit()
  elif opt in ("-i", "--input_gt"):
    genotype_file= arg
  elif opt in ("-o", "--output"):
    predict_file= arg
  elif opt in ("-q", "--quick_model"):
    quick_model= arg
  elif opt in ("-s", "--split_softclipped_sort"):
    hparams.split_softclipped= arg
    print(arg)
  elif opt in ("-t", "--threshold"):
    threshold_score= arg
  elif opt in ("-r", "--reference"):
    hparams.reference= arg
    print(arg)
def run_test(hparams,load_model,load_model_ckpt,genotype_file,predict_file, quick_model,use_existing_data=True,seed=1):
  random.seed(seed)
  tf.random.set_seed(seed)
#  genotype_file=str(sys.argv[1])
#  predict_file=str(sys.argv[2])
 # os.system('mv HG002_combine.tfrecord test.tfrecord')
  filepath=hparams.model_base_dir+"/weights/val_best_model"
  filepath=hparams.model_base_dir+"/weights/final_model_new"
  filepath=hparams.model_base_dir+"/weights/hist"
  filepath=hparams.model_base_dir+"/weights/val_best_model"
  if load_model=='last':
    if os.path.isfile(filepath+'/variables/variables.data-00000-of-00001'):
      print("Load existed best model ......")   
      model = keras.models.load_model(filepath)
      print("Load existed best model success !")
    else:
      print("model not existed ......")
      sys.exit('load_model is null')
  elif load_model=='best':
    if os.path.isfile(filepath+'/variables/variables.data-00000-of-00001'):
      print("Load existed best model ......")   
      model = keras.models.load_model(filepath)
    else:
      print("model not existed ......")   
      sys.exit('load_model is null')
  elif load_model=='load_weight':
    print("load weight ......")
    optimizer = tf.keras.optimizers.Adam(lr=hparams.learning_rate)
    model = build_model(hparams)
    model.compile(
        optimizer=optimizer,
        loss=tf.keras.losses.sparse_categorical_crossentropy,
        metrics=['accuracy'])
    if load_model=='load_weight':
      model.load_weights(filepath+'/'+load_model_ckpt)
  else:
    sys.exit('load_model is null')
  genotype_total=[]
  with open(genotype_file, 'r', encoding='utf-8') as f:
    for genotype_each in f:
      genotype_total.append(genotype_each)
  genotype_total=np.array(genotype_total)
  group_len=6000
  genotype_notFill=list(np.array(genotype_total[int(len(genotype_total)/group_len)*group_len:]))
  genotype_total=genotype_total[0:int(len(genotype_total)/group_len)*group_len]
  genotype_total=list(np.reshape(genotype_total,(-1,group_len)))
  if len(genotype_notFill)>0:
    genotype_total.append(genotype_notFill)
#  print('###')
#  print(genotype_notFill)
#  print(genotype_total)
#  print('###')
  predict=[]
  record=[]
  record_use_total=[]
  for genotype_batch in genotype_total:
    features=[]
    t,record_use,feature=make_alu_examples(hparams,genotype_batch,0,95)
    features=np.array(feature).reshape([-1,hparams.height,391,6])
#    continue
#    print(features.shape)
    features_pad=features
    test_dataset = tf.data.Dataset.from_tensor_slices((features[:,:,20:371,:], np.array([0]*t)))
    test_dataset = test_dataset.batch(batch_size=128)
    predict_ori=model.predict(test_dataset)
    features_select=[]
    record_use_select=[]
    for predict_ori_each,record_use_ori_each,feature_ori_each in zip(predict_ori,record_use,features):
        if predict_ori_each[0]<0.6:
          features_select.append(feature_ori_each)
          record_use_select.append(record_use_ori_each)
    print('original site check finished')
#    predict.append(predict_ori)
#    print(genotype_batch.shape)
    if len(features_select) ==0 :
      with open(predict_file, 'w', encoding='utf-8') as res_f:
        with open(genotype_file, 'r', encoding='utf-8') as f:
          for genotype_each in f:
            chr,pos,gt,data_type,me=genotype_each.strip().split('\t')
            res_f.write(genotype_each.strip()+'\t0\t0\t0\t0\t0\t0\n')
      sys.exit('not find any insertion')

    features_pad=np.array(features_select)
#    print(features_pad.shape)
    pad_range=[20,21,19,22,18,23,17,24,16,25,15,26,14,27,13,28,12,29,11,30,10,31,9,32,8,33,7,34,6,35,5,36,4,37,3,38,2,39,1]
    if quick_model=='1':
      pad_range=[20,21,19,22,18,23]
      print('quick mode starting ...')
    for padding in pad_range:
      print('padd:'+str(padding))
      t=len(features_pad)
      test_dataset = tf.data.Dataset.from_tensor_slices((features_pad[:,:,padding:(351+padding),:], np.array([0]*t))).batch(batch_size=128)
#      print('new image generate finished')
      predict_pad=model.predict(test_dataset)
      pad_each_i=0
      pad_i=[]
      for predict_pad_each,record_pad_each in zip(predict_pad,record_use_select):
        if (predict_pad_each[1]+predict_pad_each[2]) >=1 or -10*math.log10(1-predict_pad_each[1]-predict_pad_each[2]) > float(threshold_score) or 0.6<predict_pad_each[0]:
          predict.append(predict_pad_each)
          record.append(record_pad_each)
          pad_i.append(pad_each_i)
        pad_each_i=pad_each_i+1
      features_pad=np.delete(features_pad,pad_i,0)
      record_use_select=np.delete(record_use_select,pad_i,0)
      if len(features_pad) ==0:
        break

  with open(predict_file, 'w', encoding='utf-8') as res_f:
    for o_predict_each,o_record_each in zip(predict,record):
        chr,pos,gt,data_type,me=o_record_each.strip().split('\t')
        j=np.argmax(o_predict_each)
        max_t=max(o_predict_each)
        if (o_predict_each[1]+o_predict_each[2])>=1 or max_t>=1:
          phred_score_find=87
        else:
          phred_score_find=-10*math.log10(1-o_predict_each[1]-o_predict_each[2])
        if max_t>=1:
          phred_score_max=87
        else:
          phred_score_max=-10*math.log10(1-max_t)
        res_f.write(str(o_record_each.strip())+'\t'+str(j)+'\t'+str(phred_score_max)+'\t'+str(phred_score_find)+'\t'+str(me)+'\t'+str(o_predict_each[0])+'\t'+str(o_predict_each[1])+'\t'+str(o_predict_each[2])+'\n')


'''
  with open(genotype_file, 'r', encoding='utf-8') as f:
    for record_each in f:
      chr,pos,gt,data_type=record_each.strip().split('\t')
#      if data_type == 'HG002':
      detail_info_test.append(record_each.strip())
  with open(predict_file, 'a+', encoding='utf-8') as res_f:
    for data_each in test_dataset_0:
  #    break
      j=np.argmax(predict[i,])
      max_t=max(predict[i,])
      info_t=detail_info_test[i]
      if (predict[i,1]+predict[i,2])>=1 or max_t>=1:
         phred_score_find=87
         phred_score_max=87
      else:
         phred_score_find=-10*math.log10(1-predict[i,1]-predict[i,2])
         phred_score_max=-10*math.log10(1-max_t)
      res_f.write(info_t+'\t'+str(j)+'\t'+str(phred_score_max)+'\t'+str(phred_score_find)+'\t'+str(predict[i,0])+'\t'+str(predict[i,1])+'\t'+str(predict[i,2])+'\n')
      feature,label=data_each
      #print(feature)
      np.savetxt('tensor'+str(i)+'0.tsv',feature[0:100,0:351,0].numpy())
      np.savetxt('tensor'+str(i)+'1.tsv',feature[0:100,0:351,1].numpy())
      np.savetxt('tensor'+str(i)+'2.tsv',feature[0:100,0:351,2].numpy())
      np.savetxt('tensor'+str(i)+'3.tsv',feature[0:100,0:351,3].numpy())
      np.savetxt('tensor'+str(i)+'4.tsv',feature[0:100,0:351,4].numpy())
      np.savetxt('tensor'+str(i)+'5.tsv',feature[0:100,0:351,5].numpy())
      chr,pos,gt,data_type=detail_info_test[i].split('\t')
      plt.figure(figsize=(10,10))
      plt.subplot(k,1,1)
      plt.imshow(feature[0:100,0:351,0], cmap='gray')
      plt.subplot(k,1,2)
      plt.imshow(feature[0:100,0:351,1], cmap='gray')
      plt.subplot(k,1,3)
      plt.imshow(feature[0:100,0:351,2], cmap='gray')
      plt.subplot(k,1,4)
      plt.imshow(feature[0:100,0:351,3], cmap='gray')
      plt.subplot(k,1,5)
      plt.imshow(feature[0:100,0:351,4], cmap='gray')
      plt.subplot(k,1,6)
      plt.imshow(feature[0:100,0:351,5], cmap='gray')
      plt.xlabel(str(total_label[i])+'-'+str(j))
      plt.savefig('not_match_image/'+str(chr)+':'+str(pos)+'_'+str(total_label[i])+'-'+str(j)+'-'+str(i)+'.png')
'''
def build_model(hparams):
  """Convolutional neural network architecture."""

  model=tf.keras.applications.InceptionV3(
      include_top=True, weights=None, input_tensor=None,
      input_shape=(100,351,6), pooling=None, classes=3,
      classifier_activation='softmax'
  )



  return model


#run_test(hparams,'load_weight','cp-0460.ckpt',True)
run_test(hparams,'last','cp-0088.ckpt',genotype_file,predict_file,quick_model,True,1)
