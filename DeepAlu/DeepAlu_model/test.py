#import tensorflow as tf 
import os
import sys
#os.environ["CUDA_VISIBLE_DEVICES"]="-1"
#gpus=tf.config.list_physical_devices('CPU')
#gpus=tf.config.list_physical_devices('GPU')
#print(gpus)
#from tensorflow.python.client import device_lib
#device_lib.list_local_devices()
try:
    arg=sys.argv[1]
    print(sys.argv[1])
except:
    print(sys.argv[0])
