from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import tensorflow as tf
import pysam
import sys
import random
import os.path
import numpy as np
#!pip uninstall protobuf -y
#!pip install tensorflow
#!pip install tensorflow --upgrade --force-reinstall
#!pip install protobuf==3.3.0
#!pip3 uninstall protobuf -y
#!pip3 install -U protobuf 
#%load_ext tensorboard
import datetime, os
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import callbacks
from tensorflow.keras.callbacks import ModelCheckpoint
#from nucleus.io.genomics_writer import TFRecordWriter
from tensorflow.io import TFRecordWriter


