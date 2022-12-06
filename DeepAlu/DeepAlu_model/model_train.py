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

  if not use_existing_data:
    print('Generating data ...')
    hparams.reference='reference/Homo_sapiens_assembly38.fasta'
    generate_tfrecord_datasets(hparams)
#    print('saving data ...')
#    os.system('cp *.tfrecord /content/drive/Shareddrives/ibicqupt/deepalu_data/input_tfrecord')
    #sys.exit()


#  sys.exit()
  if os.path.isfile(_TRAIN):
    print("Train dataset exist!") 
  else:
    sys.exit()
#    print("Transfer dataset from google drive ......")
#    os.system('cp /content/drive/Shareddrives/ibicqupt/deepalu_data/input_tfrecord/*.tfrecord .')
#    print("Transfer finished ......")

  train_dataset = get_dataset(
      hparams, filename=_TRAIN, num_epochs=10)
  eval_dataset = get_dataset(
      hparams, filename=_EVAL, num_epochs=10)
  test_dataset=eval_dataset
#  test_dataset = get_dataset(
#      hparams, filename=_TEST, num_epochs=10)
  train_dataset = train_dataset.shuffle(buffer_size=10000, reshuffle_each_iteration=True).batch(batch_size=256)
  eval_dataset = eval_dataset.shuffle(buffer_size=10000, reshuffle_each_iteration=True).batch(batch_size=108)
  test_dataset = test_dataset.shuffle(buffer_size=10000, reshuffle_each_iteration=True).batch(batch_size=128)
  #NA12878_dataset = get_dataset(
  #    hparams, filename=os.path.join(hparams.out_dir, 'NA12878.tfrecord'), num_epochs=10)
  #NA12878_dataset = NA12878_dataset.repeat(hparams.total_epochs).shuffle(buffer_size=10000, reshuffle_each_iteration=True).batch(batch_size=128)
   #print("1")
  #for image_features in test_dataset:
  #  image_raw = image_features['base'].numpy()
 #   print(image_raw)
  os.environ["CUDA_VISIBLE_DEVICEs"]="4,5,6,7"
#  filepath_best=hparams.model_base_dir+"/weights/val_best_model"
#  filepath_last=hparams.model_base_dir+"/weights/final_model"

  filepath_best="weights/val_best_model"
  filepath_last="weights/final_model"
  strategy=tf.distribute.MirroredStrategy()
  with strategy.scope():
    if load_model=='best':
      if os.path.isfile(filepath_best+'/variables/variables.data-00000-of-00001'):
        print("Load existed best model ......")   
        model = keras.models.load_model(filepath_best)
      else:
        sys.exit('best model not existed!')
    elif load_model=='last' :
      if os.path.isfile(filepath_last+'/variables/variables.data-00000-of-00001'):
        print("Load existed last model ......")   
        model = keras.models.load_model(filepath_last)
      else:
        sys.exit('final model not existed!')
    elif load_model=='load_weight' or load_model=='new' :
      print("Build model at first time ......")
      optimizer = tf.keras.optimizers.Adam(learning_rate=hparams.learning_rate)
      #tensorboard_callback = tf.keras.callbacks.TensorBoard(
      #    hparams.log_dir, histogram_freq=1)
      model = build_model(hparams)
    #  model = build_model()
      model.compile(
          optimizer=optimizer,
          loss=tf.keras.losses.sparse_categorical_crossentropy,
          metrics=['accuracy'])
      if load_model=='load_weight':
        model.load_weights("weights/hist/"+load_model_ckpt)
    else:
      sys.exit('load_model is null')

  checkpoint= ModelCheckpoint(filepath_best, monitor='val_accuracy',save_weights_only=False, verbose=1, save_best_only=True, mode='max')
  filepath_hist="weights/hist/cp-{epoch:04d}.ckpt"
  checkpoint_hist= ModelCheckpoint(filepath_hist,save_weights_only=True,save_freq='epoch', verbose=1)
  callbacks_list= [checkpoint]
 # model.summary()
  logdir = os.path.join("logs", datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
  print('Training the model ......')
  auto_decay =tf.keras.callbacks.ReduceLROnPlateau(
    monitor='val_loss', factor=0.1, patience=10, verbose=1, mode='auto',
    jmin_delta=0.0001, cooldown=0, min_lr=0 )
  model.fit(
      train_dataset,
      epochs=hparams.total_epochs,
##      steps_per_epoch=13,
      validation_data=eval_dataset,
      callbacks=[callbacks_list,checkpoint_hist,auto_decay],
      verbose=2)
  model.save(filepath_last)
  print('Training complete. Obtaining final metrics...')
  eval_metrics = model.evaluate(eval_dataset, verbose=2)
  test_metrics = model.evaluate(test_dataset, verbose=2)
 # model.predict()
  print('Final eval metrics - loss: %f - accuracy: %f' %
        (eval_metrics[0], eval_metrics[1]))
  print('Final test metrics - loss: %f - accuracy: %f' %
        (test_metrics[0], test_metrics[1]))
  

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

