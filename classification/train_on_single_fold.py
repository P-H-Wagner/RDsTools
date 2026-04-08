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

#plotting colors
colors = [
"red",
"green",
"blue",
"yellow",
"orange",
"purple",
"pink",
"cyan",
"magenta",
"brown"
]

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
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
  model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
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

  epochs = 100
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


def plotMetric(history, history_key, title, ylabel, fold = None):
  '''
    Plot the loss/acc metric for training and validation sets
  '''
  # folds can have different length due to early stopping!!!
  # TOTAL LOSS
  toPlot = history.history[history_key]
  epochs = range(1, len(toPlot)+1)
  h = plt.plot(epochs, toPlot, colors[int(fold)], label=f'Fold {fold}')
  plt.subplots_adjust(left=0.2, right=0.95, top=0.85, bottom=0.20)
  plt.title(title)
  plt.xlabel('Epochs')
  plt.ylabel(ylabel)
  plt.legend()
  saveFig(plt, history_key.replace(" ", "_") + f"_fold{fold}")
  saveFig(plt, history_key.replace(" ", "_") + f"_fold{fold}")
  #log
  ymin = 0
  ymax = np.max(toPlot)
  plt.yscale('log')  
  ax = plt.gca()
  ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
  ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=20))
  ax.yaxis.set_major_formatter(LogFormatter(base=10.0, labelOnlyBase=False))
  plt.ylim(bottom=ymin * 0.9, top=ymax * 1.1)
  #plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
  saveFig(plt, history_key.replace(" ", "_") + "_log" + f"_fold{fold}")
  saveFig(plt, history_key.replace(" ", "_") + "_log" + f"_fold{fold}")
  #plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))
  plt.clf()
  plt.close()
  return epochs, toPlot


def plotKSTest(model, xx_train, y_train, xx_val, y_val,sig, fold):
  '''
    Plot the outcome of the Kolmogorov test
    Used to test the overfitting
  '''
  #class predictions 1D
  #is of shape [1,2,5,4,3,2,4,2,1, ....] 
  true_1d_val = np.argmax(y_val, axis=1)
  #select only events where the true class is chan...
  x_chan_val    = xx_val[ true_1d_val == sig ]
  #...and predict their score! (data is already scaled!)
  y_chan_val    = model.predict(x_chan_val)[:,sig]
  ## TRAINING
  #class predictions 1D
  #is of shape [1,2,5,4,3,2,4,2,1, ....] 
  true_1d_train = np.argmax(y_train, axis=1)
  #select only events where the true class is chan...
  x_chan_train    = xx_train[ true_1d_train == sig ]
  #...and predict their score! (data is already scaled!)
  y_chan_train    = model.predict(x_chan_train)[:,sig]
  #import pdb
  #pdb.set_trace()  
  np.savetxt(f"{outdir}/scoretrain_{sig}.csv",y_chan_train,delimiter = ",")
  np.savetxt(f"{outdir}/scoretest_{sig}.csv",y_chan_val,delimiter = ",")
  h1 = ROOT.TH1F(f'train_{fold}_{sig}', f'train_{fold}_{sig}', 30, 0, 1)
  h2 = ROOT.TH1F(f'val_{fold}_{sig}', f'val_{fold}_{sig}', 30, 0, 1)
  for st,sv in zip(y_chan_train, y_chan_val):
    #remark st and sv are lists of length 6! -> only keep score of the signal we're interested in
    h1.Fill(st)
    h2.Fill(sv)
  c1=ROOT.TCanvas()
  if h1.Integral()!=0: h1.Scale(1./h1.Integral())
  if h2.Integral()!=0: h2.Scale(1./h2.Integral())
  c1.Draw()
  h1.GetXaxis().SetTitle("score")
  h1.GetYaxis().SetRangeUser(0, max([h2.GetMaximum(),h1.GetMaximum()])*1.6)
  h1.SetFillColor(ROOT.kGreen+2)
  h1.SetLineColor(ROOT.kGreen+2)
  h1.SetFillStyle(3345)
  h1.Draw('HIST')
  h2.SetLineColor(ROOT.kBlack)
  h2.SetMarkerStyle(8)
  h2.SetMarkerSize(0.5)
  h2.SetMarkerColor(ROOT.kBlack)
  h2.Draw('EP SAME')
  ks_score = h1.KolmogorovTest(h2)
  ks_value = ROOT.TPaveText(0.5, 0.76, 0.88, 0.80, 'nbNDC')
  ks_value.AddText(f'KS score of class {sig} = {round(ks_score,3)}')
  ks_value.SetFillColor(0)
  ks_value.Draw('EP SAME')
  leg = ROOT.TLegend(.55,.82,.83,.88)
  leg.AddEntry(h1 ,'Training' ,'F' )
  leg.AddEntry(h2 ,'Validation' ,'EP' )
  leg.Draw("SAME")
  c1.SaveAs(outdir + f'/KS_test_{sig}_fold_{fold}.pdf')
  c1.SaveAs(outdir + f'/KS_test_{sig}_fold_{fold}.png')
  #print('KS score: ',ks_score, len(train_pred),len(test_pred))
  c1.Close()
  del c1
  return h1, h2

def plotScore(model, xx_train, y_train, xx_val, y_val, sig, class_label, fold = None):
  '''
    Plot the score of all channels class sig
  '''
  channels = [0,1,2,3,4,5]
  col ={0:'c',1:'g',2:'m',3:'r',4:'b',5:'k'}
  linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}
  hist_content = {}
  fig = plt.figure()
  hists = {}
  for chan in channels:
    dummy = []
    print("chan is ", chan, ", fold is ", fold, ", sig is ", sig)
    #class predictions 1D
    #is of shape [1,2,5,4,3,2,4,2,1, ....] 
    true_1d = np.argmax(y_val, axis=1)
    #select only events where the true class is chan...
    x_chan    = xx_val[ true_1d == chan ]
    #...and predict their score! (data is already scaled!)
    y_chan    = model.predict(x_chan)[:,sig]
    # Plot data into histo
    hists[chan] = plt.hist(y_chan, bins=np.arange(0,1.025,0.025), color=col[chan], alpha=0.5, label=class_label[chan], histtype='stepfilled',density = True,linestyle = linestyles[chan], linewidth = 1.5)
  plt.legend(loc='upper right')
  plt.title(f'Score for class ' + class_label[sig] + " (Average over folds)")
  plt.xlabel('score')
  plt.ylabel('events')
  #fig.savefig('outputs/score_' + str(sig) + '.pdf')
  #fig.savefig('outputs/score_' + str(sig) + '.png')
  saveFig(fig, f'score_' + str(sig) + f"_fold_{fold}" )
  plt.clf()
  plt.close()
  return hists


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

#save history and model into pck
with open(f"nn_training/{dt}/history_{fold}.pck", "wb") as f:
  pickle.dump(history, f)

#epochs, toPlot_loss     = plotMetric(history, "loss"      ,"Training Loss"       , "Loss"   , fold = fold)
#epochs, toPlot_acc      = plotMetric(history, "acc"       ,"Training Accuracy"   , "Acc."   , fold = fold)
#epochs, toPlot_val_loss = plotMetric(history, "val_loss"  ,"Validation Loss"     , "Loss"   , fold = fold)
#epochs, toPlot_val_acc  = plotMetric(history, "val_acc"   ,"Validation Accuracy" , "Acc."   , fold = fold)
#
#h1_0,h2_0 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 0, fold = fold)
#h1_1,h2_1 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 1, fold = fold)
#h1_2,h2_2 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 2, fold = fold)
#h1_3,h2_3 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 3, fold = fold)
#h1_4,h2_4 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 4, fold = fold)
#h1_5,h2_5 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 5, fold = fold)
#
#score_0 = plotScore (model, xx_train, y_train, xx_val, y_val, 0, class_label, fold = fold)
#score_1 = plotScore (model, xx_train, y_train, xx_val, y_val, 1, class_label, fold = fold)
#score_2 = plotScore (model, xx_train, y_train, xx_val, y_val, 2, class_label, fold = fold)
#score_3 = plotScore (model, xx_train, y_train, xx_val, y_val, 3, class_label, fold = fold)
#score_4 = plotScore (model, xx_train, y_train, xx_val, y_val, 4, class_label, fold = fold)
#score_5 = plotScore (model, xx_train, y_train, xx_val, y_val, 5, class_label, fold = fold)






