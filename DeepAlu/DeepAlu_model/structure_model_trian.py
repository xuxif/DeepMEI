from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import tensorflow as tf
import pysam
import sys
import random
import os.path
import numpy as np
import datetime, os
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import callbacks
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.io import TFRecordWriter
from deepalu_base import *
#from deepalu_base import PileupImageOptions
#from deepalu_base import ImageRow


_TRAIN = 'train.tfrecord'
_EVAL = 'eval.tfrecord'
_TEST = 'test.tfrecord'        

def run(hparams, use_existing_data,load_model,load_model_ckpt, seed=1):
  """Creates a model, runs training and evaluation."""

  # Set seed for reproducibility.  
  random.seed(seed)
  tf.random.set_seed(seed)

  os.environ["CUDA_VISIBLE_DEVICEs"]="4,5,6,7"

  strategy=tf.distribute.MirroredStrategy()
  with strategy.scope():
      print("Build model at first time ......")
      optimizer = tf.keras.optimizers.Adam(learning_rate=hparams.learning_rate)
      model = build_model(hparams)
      model.compile(
          optimizer=optimizer,
          loss=tf.keras.losses.sparse_categorical_crossentropy,
          metrics=['accuracy'])
  model.summary()

def build_model(hparams):
  """Convolutional neural network architecture."""
  
  model=tf.keras.applications.InceptionV3(
      include_top=True, weights=None, input_tensor=None,
      input_shape=(100,351,6), pooling=None, classes=3,
      classifier_activation='softmax'
  )

  return model

run(options_,use_existing_data=True,load_model='new',load_model_ckpt='')
#run(options_,use_existing_data=False,load_model='best',load_model_ckpt='')

