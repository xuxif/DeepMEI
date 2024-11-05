# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import tensorflow as tf
import tensorflow 
import pysam
import sys
import getopt
import random
import math
import os.path
import numpy as np
import datetime, os
from tensorflow import keras
from keras import layers
from keras import callbacks
from keras.callbacks import ModelCheckpoint
from tensorflow.io import TFRecordWriter
from deepalu_base_parall import *
os.environ["CUDA_VISIBLE_DEVICES"]="-1"


_TRAIN = 'train.tfrecord'
_EVAL = 'eval.tfrecord'
_TEST = 'test.tfrecord'        
hparams=PileupImageOptions(3)
def run_test(hparams,load_model,load_model_ckpt,seed):
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
      print("Load existing best model ......")   
      model = keras.models.load_model(filepath)
      print("Successfully loaded the existing best model!")
    else:
      print("model do not existed ......")
      sys.exit('load_model is null')
  elif load_model=='best':
    if os.path.isfile(filepath+'/variables/variables.data-00000-of-00001'):
      print("Load existing best model ......")   
      model = keras.models.load_model(filepath)
    else:
      print("model do not existed ......")   
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
    feature=[]


#run_test(hparams,'load_weight','cp-0460.ckpt',True)
run_test(hparams,'last','cp-0088.ckpt',7)
