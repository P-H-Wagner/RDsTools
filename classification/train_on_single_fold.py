import os
import sys
from os import path
import sys
import psutil
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf') #used s.t. it can be assigned to the batch
import matplotlib.pyplot as plt
from itertools import product
from time import time
from datetime import datetime
from sklearn.preprocessing import LabelBinarizer
from tensorflow.keras import regularizers
#import shap
from contextlib import redirect_stdout
import pdb

#import seaborn as sns
from matplotlib.ticker import LogLocator, LogFormatterExponent
from matplotlib.ticker import LogFormatter


# do this before importing tf :)

# ----- For CPU usage ----
os.environ["OMP_NUM_THREADS"] = "8"          # OpenMP (used by NumPy, Eigen)
os.environ["TF_NUM_INTRAOP_THREADS"] = "8"   # TensorFlow internal threading
os.environ["TF_NUM_INTEROP_THREADS"] = "2"   # Parallel ops

# ----- For GPU usage ----
#use one gpu :) (make only one visible)
#os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"


import tensorflow as tf
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
from root_numpy import tree2array

import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Input, Dropout, BatchNormalization
from keras.utils import plot_model
from keras.callbacks import EarlyStopping, Callback, ReduceLROnPlateau, ModelCheckpoint
from keras.constraints import unit_norm
from keras.optimizers import SGD, Adam

import sklearn as sk
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, confusion_matrix, auc
from sklearn.preprocessing import RobustScaler, StandardScaler, MinMaxScaler,MaxAbsScaler

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from helper import *
import argparse
import seaborn
import yaml

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('-dt', '--datetime'  , required=True)
parser.add_argument('-f' , '--fold'      , required=True)
parser.add_argument('-d' , '--debug'     , action='store_true' )

args = parser.parse_args()
dt   = args.datetime
fold = args.fold

#toSave
outdir = f"/work/pahwagne/RDsTools/classification/nn_training/{dt}/"
os.system(f"mkdir -p {outdir}")

#class naming logic
mc_ids      = [0       ,1        ,10         ,11          ,-1]
mc_classes  = [2       ,0        ,3          ,1           ,4 ] #changed ! 
label       = {0: r" \mu", 1:r" \tau", 10:r" \mu^{*}", 11: r" \tau^{*}", -1: r" H_{b}" }
#class_label = {class_id: [] for class_id in mc_classes}

class_label = {}
for mc_id, class_id in zip (mc_ids,mc_classes):
  class_label[class_id] = label[mc_id]
class_label[5] = "Comb. Bkg."


# ----- For CPU usage ----
# Set thread counts before running anything
tf.config.threading.set_intra_op_parallelism_threads(8)
tf.config.threading.set_inter_op_parallelism_threads(2)

# Optional: check what they were set to
print("Intra-op threads:", tf.config.threading.get_intra_op_parallelism_threads())
print("Inter-op threads:", tf.config.threading.get_inter_op_parallelism_threads())

# ----- For GPU usage ----
# For all visible gpu (0), the memory should be allocated on demand,
# not per se from the beginning

gpus = tf.config.list_physical_devices('GPU') #should be the one set we os.environ
#tf.config.run_functions_eagerly(True)

if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print("Error setting memory growth:", e)

print("Using GPUs:", gpus)

#print some info
from tensorflow.python.platform import build_info as tf_build_info
print(tf_build_info.build_info)


def defineModel(xIn):
  '''
    Define the NN
  '''
  #NOTE for the moment, everything is hardcoded

  rate_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
  initial_learning_rate=0.00001,
  decay_steps=10000,
  decay_rate=0.9
  )

  l2_rate = regularizers.l2(0.005)
  learning_rate = 0.00005
  model = tf.keras.Sequential()
  model.add(tf.keras.layers.Input((xIn,)))
  model.add(tf.keras.layers.Dense(64,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  #model.add(tf.keras.layers.Dropout(0.2))
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  #model.add(tf.keras.layers.Dropout(0.2))
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  model.add(tf.keras.layers.Dense(64,  activation ='swish'  ,  kernel_regularizer=l2_rate))

  #model.add(tf.keras.layers.Dense(64, activation='swish'))
  #model.add(tf.keras.layers.Dropout(0.2))
  #model.add(tf.keras.layers.Dense(64, activation='swish'))
  #model.add(tf.keras.layers.Dropout(0.2))
  #model.add(tf.keras.layers.Dense(64,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  #model.add(tf.keras.layers.Dropout(0.2))
  #model.add(tf.keras.layers.Dense(64,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  #model.add(tf.keras.layers.Dropout(0.2))
  #model.add(tf.keras.layers.Dense(64,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  #model.add(tf.keras.layers.Dropout(0.2))
  #model.add(tf.keras.layers.Dense(24, activation='relu'))


  model.add(tf.keras.layers.Dense(6  ,  activation= 'softmax'                                            ))
  opt = keras.optimizers.Adam(learning_rate=learning_rate)
  model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=['acc'])

  print(model.summary()) #not enough

  def print_model_summary(model):

    print("============= DETAILED MODEL SUMMARY ===============")
    print("learning rate = ", learning_rate)
    for i, layer in enumerate(model.layers):
      config = layer.get_config()
      print(f"Layer {i+1} ({layer.name})")
      #print(f"  Type: {config['class_name']}")
      print(f"  Output Shape: {layer.output_shape}")
      print(f"  Activation: {config.get('activation', 'N/A')}")

      if 'kernel_regularizer' in config:
          reg = config['kernel_regularizer']
          if reg:
              print(f"  Regularizer: {reg['class_name']}, (rate: {l2_rate})")
      print(f"  Parameters: {layer.count_params()}")
      print("-" * 60)

  print_model_summary(model)

  #write model into file
  f = open(outdir + "/settings.txt", "a")
  with redirect_stdout(f):
    model.summary()
    print_model_summary(model)
  f.close()

  return model


def loadData(dt,fold):
  
  path = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/{dt}" 

  with open(f"{path}/xx_train_{fold}.pck", 'rb') as f:
    xx_train = pickle.load(f)
  with open(f"{path}/y_train_{fold}.pck", 'rb') as f:
    y_train  = pickle.load(f)
  with open(f"{path}/w_train_{fold}.pck", 'rb') as f:
    w_train  = pickle.load(f)

  with open(f"{path}/xx_val_{fold}.pck", 'rb') as f:
    xx_val   = pickle.load(f)
  with open(f"{path}/y_val_{fold}.pck", 'rb') as f:
    y_val    = pickle.load(f)
  with open(f"{path}/w_val_{fold}.pck", 'rb') as f:
    w_val    = pickle.load(f)

  with open(f"{path}/class_w.pck", 'rb') as f:
    class_w  = pickle.load(f)

  return xx_train, y_train, w_train, xx_val, y_val, w_val, class_w 


def defineCallbacks(fold):
  '''
    Define the callbacks
  '''
  # early stopping
  monitor = 'val_loss'
  es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=3000)
  
  # reduce learning rate when at plateau, fine search the minimum
  reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=5, min_lr=0.000001, cooldown=10, verbose=True)
  
  # save the model every now and then
  filepath = '/'.join([outdir, f'fold_{fold}' + '_saved-model-{epoch:04d}_val_loss_{val_loss:.4f}_val_acc_{val_acc:.4f}.h5'])
  save_model = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=False, mode='auto', period=1)
  
  callbacks = [save_model]
  callbacks.append(es)
  #callbacks.append(reduce_lr)
  #callbacks.append(CustomCallback())

  return callbacks



def getParameters():

  epochs = 1000 
  batch_size = 2000
  callbacks = defineCallbacks(fold)


  path = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/{dt}" 
  #get the dimenision of the input layer
  with open(f"{path}/len_features.pck", 'rb') as f:
    len_features = pickle.load(f)


  return epochs, batch_size, callbacks, len_features

def saveFig(plt, name):
  '''
    Save python figure
  '''
  plt.savefig(f'{outdir}/{name}.pdf')
  plt.savefig(f'{outdir}/{name}.png')
  print(f' ========> {outdir}/{name}.png created')


##################################################################

#define parameters
epochs, batch_size, callbacks, len_features = getParameters() 
#define model
model = defineModel(len_features)

#load data to train on for fold {fold}
xx_train, y_train, w_train, xx_val, y_val, w_val, class_w = loadData(dt,fold)

print("====> Flatten weights ...")
w_train_flat = w_train.flatten().astype(np.float32)
w_val_flat   = w_val.flatten().astype(np.float32)
print("====> Flatten done ...")

print("====> Convert to tensor...")
w_train_tf = tf.convert_to_tensor(w_train_flat, dtype=tf.float32)
w_val_tf   = tf.convert_to_tensor(w_val_flat, dtype=tf.float32)
print("====> Conversion done ...")


history = model.fit(xx_train, y_train, validation_data=(xx_val, y_val, w_val_tf), epochs=epochs, callbacks=callbacks, batch_size=batch_size, verbose=True,class_weight = class_w , sample_weight = w_train_tf )

#save history and model into pck
with open(f"nn_training/{dt}/history_{fold}.pck", "wb") as f:
  pickle.dump(history, f)






