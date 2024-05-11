from __future__ import absolute_import
from __future__ import absolute_import, division, print_function
from __future__ import division
from __future__ import print_function
from collections import namedtuple
from multiprocessing import Pool
from tensorflow import keras
from tensorflow.io import TFRecordWriter
from tensorflow.keras import callbacks
from tensorflow.keras import layers
from tensorflow.keras.callbacks import ModelCheckpoint
import datetime
import getopt
import math
import numpy as np;
import os
import pysam
import random
import sys
import tensorflow as tf
import time

