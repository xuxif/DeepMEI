# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import tensorflow as tf
import pysam
import sys
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
import matplotlib.pyplot as plt
os.environ["CUDA_VISIBLE_DEVICES"]="-1"


_TRAIN = 'train.tfrecord'
_EVAL = 'eval.tfrecord'
_TEST = 'test.tfrecord'        
hparams=PileupImageOptions(3)
#hparams.reference='/DeepAlu/DeepAlu_model/reference/Homo_sapiens_assembly38.fasta'
def run_test(hparams):
  hparams.split_softclipped='/DeepAlu/data_cluster/split_softclipped_sort_test'
  output_dir=str(sys.argv[1])
  features=[]
  record_total=[]
  with open(sys.argv[2], 'r', encoding='utf-8') as geno_f:
    for record_each in geno_f:
      feature,record_each_t =make_alu_examples_each(hparams,record_each)
      features.append(np.array(feature).reshape([hparams.height,351,6]))
      record_total.append(record_each_t)
  
  t=len(features)
  test_dataset = tf.data.Dataset.from_tensor_slices((features, np.array([0]*t)))
  filepath=hparams.model_base_dir+"/weights/val_best_model"
  model = keras.models.load_model(filepath)
  test_dataset = test_dataset.batch(batch_size=128)
  predict_ori=model.predict(test_dataset)
  for predict_info,record_each,feature in zip(predict_ori,record_total,features):
    j=np.argmax(predict_info)
#    print(predict_info)
    max_t=max(predict_info)
    if (predict_info[1]+predict_info[2])>=1:
      phred_score_find=87
    else:
      phred_score_find=-10*math.log10(1-predict_info[1]-predict_info[2])
    if max_t==1:
      phred_score_max=87
    else:
      phred_score_max=-10*math.log10(1-max_t)
#    print(record_each)
    record_each=record_each.strip().split('\t')
    print(str(record_each[0])+'\t'+str(record_each[1])+'\t'+str(j)+'\t'+str(int(phred_score_find))+'\t'+str(int(phred_score_max))+'\t'+str(predict_info[0])+'\t'+str(predict_info[1])+'\t'+str(predict_info[2]).strip())
    k=6
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
    plt.savefig(output_dir+'/'+str(record_each[0])+':'+str(record_each[1])+'_'+str(j)+'_'+str(int(phred_score_find))+'.png')
run_test(hparams)
