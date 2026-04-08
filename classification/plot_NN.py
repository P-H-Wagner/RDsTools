import glob
import re
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

args = parser.parse_args()
dt   = args.datetime

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


def loadData(dt,fold):
  
  path      = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/{dt}" 
  path_work = f"/work/pahwagne/RDsTools/classification/nn_training/{dt}" 

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

  with open(f"{path_work}/history_{fold}.pck", 'rb') as f:
    history = pickle.load(f)


  return xx_train, y_train, w_train, xx_val, y_val, w_val, history


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



def plotAverageMetric(epochs, toPlot, history_key, title, ylabel, nfolds):
  # TOTAL LOSS
  pdb.set_trace()
  average = []
  #max epochs (get max epochs)
  max_epochs = max([len(val) for val in epochs.values()])
  for n in range(nfolds):
    plt.plot(epochs[n], toPlot[n], colors[n], label=f'Fold {n}')
    #create a full array of nans of length max_epochs       
    padded = np.full(max_epochs, np.nan)
    #fill the first part with (the rest stays nan)
    padded[:len(toPlot[n])] = toPlot[n]
    #append the padded array to the average 
    average.append(padded)
  #take the average and ignore nans
  average_toPlot = np.nanmean(np.vstack(average), axis=0)
  plt.plot(range(1, max_epochs + 1), average_toPlot, 'black' , label=f"Fold Average", linestyle= 'dashed')
  plt.subplots_adjust(left=0.2, right=0.95, top=0.85, bottom=0.20)
  plt.title(title)
  plt.xlabel('Epochs')
  plt.ylabel(ylabel)
  plt.legend()
  saveFig(plt, history_key + "_average")
  saveFig(plt, history_key + "_average")
  plt.close()

def plotAverageScore(scores, sig, class_label, nfolds):
  channels = [0,1,2,3,4,5]
  col ={0:'c',1:'g',2:'m',3:'r',4:'b',5:'k'}
  linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}
  hist_content = {}
  #plot average
  for chan in channels:
    dummy = []
    for n in range(nfolds):
      #append the bin content (at index 0) for every fold! 
      dummy.append(scores[n][chan][0])
    #average over the folds to get the average histo
    hist_content[chan] = sum(dummy)/nfolds
  fig = plt.figure()
  pdb.set_trace()
  for chan in channels:
    #bin edges are at index 1 (the same for all, can take the last one)
    bin_edges = scores[0][chan][1]
    bin_width = np.diff(bin_edges)
    # plot the score distributions
    plt.bar(bin_edges[:-1], hist_content[chan], color=col[chan], width = bin_width, align = 'edge', alpha=0.5, label=class_label[chan], linestyle = linestyles[chan], linewidth = 1.5)
  #plot all channels into same plot 
  plt.legend(loc='upper right')
  plt.title(f'Score for class ' + class_label[sig] + " (Average over folds)")
  plt.xlabel('score')
  plt.ylabel('events')
  #fig.savefig('outputs/score_' + str(sig) + '.pdf')
  #fig.savefig('outputs/score_' + str(sig) + '.png')
  saveFig(fig, "score_" + str(sig) + "_average" )
  plt.clf()
  plt.close()


def plotAverageKSTest(h1, h2, sig, nfolds):
  h1_av = ROOT.TH1F(f'train_{sig}', f'train_{sig}', 30, 0, 1)
  h2_av = ROOT.TH1F(f'val_{sig}', f'val_{sig}', 30, 0, 1)
  for n in range(nfolds):
    h1_av.Add(h1[n])
    h2_av.Add(h2[n])
  h1_av.Scale(1/nfolds)
  h2_av.Scale(1/nfolds)
  c1=ROOT.TCanvas()
  if h1_av.Integral()!=0: h1_av.Scale(1./h1_av.Integral())
  if h2_av.Integral()!=0: h2_av.Scale(1./h2_av.Integral())
  c1.Draw()
  h1_av.GetXaxis().SetTitle("score")
  h1_av.GetYaxis().SetRangeUser(0, max([h2_av.GetMaximum(),h1_av.GetMaximum()])*1.6)
  h1_av.SetFillColor(ROOT.kGreen+2)
  h1_av.SetLineColor(ROOT.kGreen+2)
  h1_av.SetFillStyle(3345)
  h1_av.Draw('HIST')
  h2_av.SetLineColor(ROOT.kBlack)
  h2_av.SetMarkerStyle(8)
  h2_av.SetMarkerSize(0.5)
  h2_av.SetMarkerColor(ROOT.kBlack)
  h2_av.Draw('EP SAME')
  ks_score = h1_av.KolmogorovTest(h2_av)
  ks_value = ROOT.TPaveText(0.5, 0.76, 0.88, 0.80, 'nbNDC')
  ks_value.AddText(f'Average KS score of class {sig} = {round(ks_score,3)}')
  ks_value.SetFillColor(0)
  ks_value.Draw('EP SAME')
  leg = ROOT.TLegend(.55,.82,.83,.88)
  leg.AddEntry(h1_av ,'Training' ,'F' )
  leg.AddEntry(h2_av ,'Validation' ,'EP' )
  leg.Draw("SAME")
  c1.SaveAs(outdir + f'/KS_test_{sig}_average.pdf')
  c1.SaveAs(outdir + f'/KS_test_{sig}_average.png')
  #print('KS score: ',ks_score, len(train_pred),len(test_pred))
  c1.Close()
  del c1



##################################################################


#get the number of folds
files = glob.glob("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/xx_train_*.pck")
folds = [int(re.search(r'_(\d+)\.pck', f).group(1)) for f in files]
max_fold = max(folds) #f.e. 9 when we have 10 folds (0, ..., 9)
nfolds = max_fold + 1


###################
# K-folding loop  #
###################
dict_epochs          = {}
dict_toPlot_loss     = {}
dict_toPlot_acc      = {}
dict_toPlot_val_loss = {}
dict_toPlot_val_acc  = {}
dict_score_0         = {}
dict_score_1         = {}
dict_score_2         = {}
dict_score_3         = {}
dict_score_4         = {}
dict_score_5         = {}
dict_h1_0            = {}
dict_h2_0            = {}
dict_h1_1            = {}
dict_h2_1            = {}
dict_h1_2            = {}
dict_h2_2            = {}
dict_h1_3            = {}
dict_h2_3            = {}
dict_h1_4            = {}
dict_h2_4            = {}
dict_h1_5            = {}
dict_h2_5            = {}

for fold in range(nfolds):


  #for every fold, get the best/newest model
  models = glob.glob(f"/work/pahwagne/RDsTools/classification/nn_training/{dt}/fold_{fold}_saved-model*.h5")

  max_nr = 0
  best_model = ""

  for model in models:

    end = model.split("saved-model-")[1]
    nr  = int(end.split("_val_loss")[0])
    #print("number now is:", nr)
    #print("max number is:", max_nr)
    if nr >= max_nr: 
      max_nr = nr
      best_model = model

  print(f"Best model for fold {fold} is: {best_model}")

  
  #load this model
  model = tf.keras.models.load_model(best_model) 


  #load data to train on for fold {fold}
  xx_train, y_train, w_train, xx_val, y_val, w_val, history = loadData(dt,fold)
  
  epochs, toPlot_loss     = plotMetric(history, "loss"      ,"Training Loss"       , "Loss"   , fold = fold)
  epochs, toPlot_acc      = plotMetric(history, "acc"       ,"Training Accuracy"   , "Acc."   , fold = fold)
  epochs, toPlot_val_loss = plotMetric(history, "val_loss"  ,"Validation Loss"     , "Loss"   , fold = fold)
  epochs, toPlot_val_acc  = plotMetric(history, "val_acc"   ,"Validation Accuracy" , "Acc."   , fold = fold)
  
  h1_0,h2_0 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 0, fold = fold)
  h1_1,h2_1 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 1, fold = fold)
  h1_2,h2_2 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 2, fold = fold)
  h1_3,h2_3 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 3, fold = fold)
  h1_4,h2_4 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 4, fold = fold)
  h1_5,h2_5 = plotKSTest(model, xx_train, y_train, xx_val, y_val, 5, fold = fold)
  
  score_0 = plotScore (model, xx_train, y_train, xx_val, y_val, 0, class_label, fold = fold)
  score_1 = plotScore (model, xx_train, y_train, xx_val, y_val, 1, class_label, fold = fold)
  score_2 = plotScore (model, xx_train, y_train, xx_val, y_val, 2, class_label, fold = fold)
  score_3 = plotScore (model, xx_train, y_train, xx_val, y_val, 3, class_label, fold = fold)
  score_4 = plotScore (model, xx_train, y_train, xx_val, y_val, 4, class_label, fold = fold)
  score_5 = plotScore (model, xx_train, y_train, xx_val, y_val, 5, class_label, fold = fold)


  dict_epochs         [fold] = epochs         
  dict_toPlot_loss    [fold] = toPlot_loss    
  dict_toPlot_acc     [fold] = toPlot_acc     
  dict_toPlot_val_loss[fold] = toPlot_val_loss
  dict_toPlot_val_acc [fold] = toPlot_val_acc 
  dict_score_0        [fold] = score_0         
  dict_score_1        [fold] = score_1         
  dict_score_2        [fold] = score_2         
  dict_score_3        [fold] = score_3         
  dict_score_4        [fold] = score_4         
  dict_score_5        [fold] = score_5         
  dict_h1_0           [fold] = h1_0           
  dict_h2_0           [fold] = h2_0           
  dict_h1_1           [fold] = h1_1           
  dict_h2_1           [fold] = h2_1           
  dict_h1_2           [fold] = h1_2           
  dict_h2_2           [fold] = h2_2           
  dict_h1_3           [fold] = h1_3           
  dict_h2_3           [fold] = h2_3           
  dict_h1_4           [fold] = h1_4           
  dict_h2_4           [fold] = h2_4           
  dict_h1_5           [fold] = h1_5           
  dict_h2_5           [fold] = h2_5           


plotAverageMetric (dict_epochs, dict_toPlot_loss    , "loss"     ,"Training Loss"       , "Loss", nfolds)
plotAverageMetric (dict_epochs, dict_toPlot_acc     , "acc"      ,"Training Accuracy"   , "Acc.", nfolds)
plotAverageMetric (dict_epochs, dict_toPlot_val_loss, "val_loss" ,"Validation Loss"     , "Loss", nfolds)
plotAverageMetric (dict_epochs, dict_toPlot_val_acc , "val_acc"  ,"Validation Accuracy" , "Acc.", nfolds)

plotAverageScore  (dict_score_0, 0, class_label, nfolds)
plotAverageScore  (dict_score_1, 1, class_label, nfolds)
plotAverageScore  (dict_score_2, 2, class_label, nfolds)
plotAverageScore  (dict_score_3, 3, class_label, nfolds)
plotAverageScore  (dict_score_4, 4, class_label, nfolds)
plotAverageScore  (dict_score_5, 5, class_label, nfolds)

plotAverageKSTest (dict_h1_0  , dict_h2_0, 0, nfolds)
plotAverageKSTest (dict_h1_1  , dict_h2_1, 1, nfolds)
plotAverageKSTest (dict_h1_2  , dict_h2_2, 2, nfolds)
plotAverageKSTest (dict_h1_3  , dict_h2_3, 3, nfolds)
plotAverageKSTest (dict_h1_4  , dict_h2_4, 4, nfolds)
plotAverageKSTest (dict_h1_5  , dict_h2_5, 5, nfolds)

