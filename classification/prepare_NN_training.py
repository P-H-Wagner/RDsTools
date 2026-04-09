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


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('-dt', '--datetime'  , required=True)
parser.add_argument('-d' , '--debug'     , action='store_true' )
args = parser.parse_args()


scratch_dir = f"/scratch/pahwagne/{args.datetime}"

#shap.initjs()

#now = datetime.now()
#dt  = now.strftime("%d_%m_%Y_%H_%M_%S")
dt = args.datetime


########################################
#                                      #
# multiclassifier with 6 classes to    #
# optimize signal to background        #
# significance                         #
#                                      #
########################################


#------------------input files-------------------


baseline_selection = minimal 

#trigger_data = "&& ((mu7_ip4 == 1) || (mu9_ip6 == 1))"
#trigger_data  = " &&((mu7_ip4==1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6==1))"
trigger_branches = [
    "mu7_ip4", 
    "mu8_ip3", 
    "mu8_ip5", 
    "mu8_ip6",
    "mu9_ip4", 
    "mu9_ip5", 
    "mu9_ip6", 
    "mu12_ip6"
]
trigger_data = ""
trigger_mc   = ""

#bdt is evaluated on the skimmed datasets :) 
base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data}/"
file_data = [base + f for f in os.listdir(base)]

#enough!
file_data = file_data[:50]  

#hammer is evaluated on the skimmed datasets :) 
base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{sig_cons_hammer_25}/" 
file_sig = [base + f for f in os.listdir(base)]

base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_25[0]}/"
file_hb = [base + f for f in os.listdir(base)]


file_data = file_data[:70]  

if args.debug:
  file_data = file_data[:5]  
  file_sig  = file_sig [:5]  
  file_hb   = file_hb  [:5]  

#--------------------define sideband region (mlow, mhigh) and dump it into pickle ----------------

sigma = 0.009

#signal region
mlow   = dsMass_ - nSignalRegion*sigma
mhigh  = dsMass_ + nSignalRegion*sigma

#sideband start
mlow2  = dsMass_ - nSidebands*sigma
mhigh2 = dsMass_ + nSidebands*sigma

#sideband stops
mlow3  = mlow2  - sbWidth*sigma
mhigh3 = mhigh2 + sbWidth*sigma

signalRegion = f"&& ({mlow} < phiPi_m) & (phiPi_m < {mhigh})"
leftSB       = f"&& ({mlow3} < phiPi_m) & (phiPi_m < {mlow2})"
rightSB      = f"&& ({mhigh2} < phiPi_m) & (phiPi_m < {mhigh3})"
lowMass      = f"&& (dsMu_m < {bsMass_})"


#--------------------define training features ----------------

features = [
    "phiPi_deltaR",
    "kk_deltaR",
    "dsMu_deltaR",
    "bs_pt_coll",
    "bs_pt_lhcb_alt",
    "bs_pt_reco_1",
    "bs_pt_reco_2",
    "dsMu_m",
    "phiPi_m",
    "q2_coll",
    "q2_lhcb_alt",
    "q2_reco_1",
    "q2_reco_2",
    "pi_pt",
    "pi_eta",
    "mu_pt",
    "mu_eta",
    "k1_pt",
    "k2_pt",
    "cosPiK1",
    "cosMuW_coll",
    "cosMuW_lhcb_alt",
    "cosMuW_reco_1",
    "cosMuW_reco_2",
    "e_star_coll",
    "e_star_lhcb_alt",
    "e_star_reco_1",
    "e_star_reco_2",
    "fv_chi2",
    "tv_chi2",
    "sv_chi2",
    "lxy_ds",
    "mu_id_medium",
    "rel_iso_03_pv",
    "mu_is_global",
    "ds_vtx_cosine_xyz_pv"

]


kin_var = [

# helicity
#'cosPhiDs_lhcb_alt',
#'cosPhiDs_reco_1',
#'cosPhiDs_reco_2',
#
'abs(cosPiK1)',

'cosMuW_lhcb_alt', 
'cosMuW_reco_1', 
'cosMuW_reco_2', 

# vertex info
'fv_chi2', # dont take prob, they can hit the prec. limit when small!
'tv_chi2', # "
'sv_chi2', # "

# displacement
'lxy_bs_sig',
'lxy_ds_sig',

'dxy_mu_sig_pv',
'dxy_mu_sig_sv',

# kinematics
'mu_pt',
'mu_eta',
'pi_pt',
'pi_eta',
'k1_pt',
'k1_eta',
'k2_pt',
'k2_eta',

# masses
'kk_m',
#'phiPi_m',
'dsMu_m',

# delta R
'kk_deltaR',
'phiPi_deltaR',
'dsMu_deltaR',

# q2 and co
'q2_coll',
'e_star_coll',

#'bs_boost_lhcb_alt',
'bs_pt_lhcb_alt',
'e_star_lhcb_alt',
#'m2_miss_lhcb_alt',
'pt_miss_lhcb_alt',
'm2_miss_lhcb_alt',

#'bs_boost_reco_1',
'bs_pt_reco_1',
'e_star_reco_1',
'q2_reco_1',
'pt_miss_reco_1',
'm2_miss_reco_1',

#'bs_boost_reco_2',
'bs_pt_reco_2',
'e_star_reco_2',
'q2_reco_2',
'pt_miss_reco_2',
'm2_miss_reco_2',

'disc_negativity',
'ds_vtx_cosine_xy_pv',
'kappa',

'signed_decay_ip3d_mu_ds_sv',
'signed_decay_ip3d_mu_bs_sv',

'rel_iso_03',
'rel_iso_03_pv',
'rel_iso_03_sv',
'rel_iso_03_ds_sv',
'rel_iso_03_tv',

]

hammer = [
"central_w",
]

trigger_sf = [
"trigger_sf",
]

class Sample(object):
  '''
    Class to convert the sample into dataframe while applying some selection
  '''

  #init gets called, everytime a class object gets defined
  
  def __init__(self, filename,tree, selection,signal):
    self.filename = filename
    self.selection = selection
    self.tree = tree
    self.signal = signal


    select_trigger = False

    if self.signal in [0,1,2,3]:
      #signals
      #print("hammer branches")
      branches = kin_var + trigger_branches + ['event'] + hammer + trigger_sf #event will be removed later!
      #branches = kin_var + ['event'] #event will be removed later!

    elif self.signal == 4: #hb
      branches = kin_var + trigger_branches + ['event'] + trigger_sf #event will be removed later!

    else: #data
      branches = kin_var + trigger_branches + ['event','sf_weights'] #event will be removed later!
      select_trigger = True

    print("====> analyzing:",self.signal)

    if isinstance(filename, list):

      pd_list = []
      for name in filename:
        f = ROOT.TFile(name)
        tree_obj = f.Get(self.tree)
        arr = tree2array(tree_obj,selection = self.selection, branches = branches) #event will be removed later!

        # NEW! select trigger here
        intermediate = pd.DataFrame(arr)

        # a pandas mask to select the trigger
        trigger_mask = np.logical_or.reduce([intermediate[trigger] == 1 for trigger in trigger_branches])  
        intermediate = intermediate[trigger_mask]

        #append
        pd_list.append(intermediate)

      self.df = pd.concat(pd_list)
      #pdb.set_trace()


    else: 
      f = ROOT.TFile(self.filename)
      tree_obj = f.Get(self.tree)
      arr = tree2array(tree_obj,selection = self.selection, branches = branches ) #event will be removed later!

      # NEW! select trigger here
      intermediate = pd.DataFrame(arr)

      # a pandas mask to select the trigger
      trigger_mask = np.logical_or.reduce([intermediate[trigger] == 1 for trigger in trigger_branches])  
      intermediate = intermediate[trigger_mask]

      self.df = intermediate 


    pd.options.mode.chained_assignment = None
    self.df['is_signal'] = signal


class Trainer(object):

  'train the data'

  def __init__(self, features, epochs, batch_size, scaler_type, do_early_stopping, do_reduce_lr, dirname, baseline_selection, nfolds, frac_sb, frac_sf):
    self.features = features #variables to train on
    self.epochs = epochs #samples / batch_size =  number of iterations to 1 epoch
    self.batch_size = batch_size #number of samples which we feed to our model
    self.scaler_type = scaler_type
    self.do_early_stopping = do_early_stopping 
    self.do_reduce_lr = do_reduce_lr
    self.dirname = dirname + '_' + datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss')
    self.baseline_selection = baseline_selection
    self.nfolds      = nfolds
    self.frac_sb = frac_sb
    self.frac_sf = frac_sf
    self.colors = [
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


  def saveSettings(self):

    #already created in submitter file
    outdir = f'./nn_training/{dt}'
    self.outdir = outdir

    f = open(self.outdir + "/settings.txt", "x")
    f.write("Features: "                  + str(self.features)            + "\n")
    f.write("Epochs: "                    + str(self.epochs)              + "\n")
    f.write("Nfolds: "                    + str(self.nfolds)              + "\n")
    f.write("Batch size: "                + str(self.batch_size)          + "\n")
    f.write("Scaler type: "               + str(self.scaler_type)         + "\n")
    f.write("Early stopping: "            + str(self.do_early_stopping)   + "\n")
    f.write("Reduce lr: "                 + str(self.do_reduce_lr)        + "\n")
    f.write("Baseline selection: "        + str(self.baseline_selection)  + "\n")
    f.write("Signal region up to:"        + str(nSignalRegion) + " sigma" + "\n")
    f.write("Sideband region starts at:"  + str(nSidebands)    + " sigma" + "\n")
    f.write("Sideband region width:"      + str(sbWidth)       + " sigma" + "\n")
    f.close()

    return 


  def saveFig(self, plt, name):
    '''
      Save python figure
    '''
    plt.savefig(f'{self.outdir}/{name}.pdf')    
    plt.savefig(f'{self.outdir}/{name}.png')    
    print(f' ========> {self.outdir}/{name}.png created')


  def getSamples(self):
    '''
      Function that fetches the samples into lists
    '''

    train = {}
    test  = {}


    print('========> start reading the trees')
    now = time()

    #-------------------------- SIGNAL and HB ------------------------------------ 
 
    tree_name = 'tree'
  
    # Lets collect everything which is not signal and not combinatorial into hb = -1 (for now)
    hb_selec = " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) && (gen_match_success) && (gen_same_mother == 1)"
 
    #signals      #ds mu   #ds tau   #dstar mu   #dstar tau   #hb
    mc_ids      = [0       ,1        ,10         ,11          ,-1]
    mc_classes  = [2       ,0        ,3          ,1           ,4 ] #changed ! 
    label       = {0: r" \mu", 1:r" \tau", 10:r" \mu^{*}", 11: r" \tau^{*}", -1: r" H_{b}" }
    class_label = {class_id: [] for class_id in mc_classes}

    for mc_id, class_id in zip (mc_ids,mc_classes):

      class_label[class_id].append(label[mc_id]) 

      if mc_id >= 0: 
        mc_sample            = Sample(filename = file_sig,    selection=self.baseline_selection  + trigger_mc + f'&& (gen_sig == {mc_id})' + lowMass ,                tree = tree_name,signal = class_id).df

      else:
        mc_sample0          = Sample(filename = file_hb,     selection=self.baseline_selection  + trigger_mc + hb_selec + lowMass ,                tree = tree_name,signal = class_id).df
        mc_sample           = pd.concat([mc_sample0], sort = False)
      
      mc_train, mc_test = train_test_split(mc_sample,test_size = 0.2,random_state = 1000)

      f = open(self.outdir + "/settings.txt", "a")
      f.write(f"train sample for signal {mc_id} has class name {class_id} and contains {len(mc_train)} events \n" )
      f.write(f"test sample for signal {mc_id} has class name {class_id} and contains {len(mc_test)} events \n"   )
      f.close()

      print(f"train sample for signal {mc_id} has class name {class_id} and contains {len(mc_train)} events")
      print(f"test sample for signal {mc_id} has class name {class_id} and contains {len(mc_test)} events")
  
      train[class_id] = mc_train
      test [class_id] = mc_test

    #---------------------------  COMB BACKGROUND ------------------------------------ 

    ## Denote comb with id 2
    data_id = 5
    class_label[5] = ["Comb. Bkg."]

    #join the labels with ','
    for key in class_label.keys(): 
      class_label[key] =  ','.join(class_label[key]) 
      class_label[key] =   "$" + class_label[key] + "$"

    #siebands used for the comb. bkg

    sign_flip   = " && ((k1_charge*k2_charge < 0) && (pi_charge*mu_charge>0))"
    #sign_flip   = " && ((k1_charge*k2_charge > 0) || (pi_charge*mu_charge>0))"

    #data_selec  = self.baseline_selection + sign_flip + signalRegion + lowMass
    data_selec  = self.baseline_selection + sign_flip + lowMass 

    data_sample_sf = Sample(filename=file_data,   selection=data_selec + trigger_data , tree = 'tree',signal = data_id).df
    data_sample  = pd.concat([data_sample_sf], sort = False)

    data_train, data_test = train_test_split(data_sample, test_size=0.2, random_state = 1000)

    train[data_id] = data_train
    test [data_id] = data_test

    f = open(self.outdir + "/settings.txt", "a")
    f.write(f"train sample for data has class name {data_id} and  has {len(data_train)} events \n")
    f.write(f"test sample for data has class name {data_id} and  has {len(data_test)} events \n"  )
    f.close()

    print(f"train sample for data has class name {data_id} and  has {len(data_train)} events")
    print(f"test sample for data has class name {data_id} and  has {len(data_test)} events")
    print(f'========> it took {round(time() - now,2)} seconds' )

    return train, test, class_label 

  def createDataframe(self, samples):
    '''
      calculates the weight, drops nans
    '''

    samples_notnan  = {}
    samples_len     = {}
    total_len       = 0
    weight          = {}

    for key in samples.keys():

      #drop nans
      dummy = samples[key].dropna()
      samples_notnan[key] = dummy

      #for weight
      samples_len[key] = len(dummy)
      total_len        += samples_len[key]
      
    
    for key in samples.keys():
      weight[key] = total_len / samples_len[key] 
      print(f"sample for class {key} contains {samples_len[key]} events (without nans)")
    print(weight) 
    return samples_notnan, weight


  def doScaling(self, X):
      '''
        Normalise the input features with a keras scaler 
      '''
      if self.scaler_type == 'robust':
        qt = RobustScaler()
      elif self.scaler_type == 'standard':
        qt = StandardScaler()
      elif self.scaler_type == 'minmax':
        qt = MinMaxScaler(feature_range=(-1,1))
      elif self.scaler_type == 'maxabs':
        qt = MaxAbsScaler()

      else:
        raise RuntimeError(f'Unknown scaler {self.scaler_type} - Aborting...')

      qt.fit(X[self.features])
      xx = qt.transform(X[self.features])

      #xx = np.clip(xx, np.percentile(xx, 0.01, axis=0),np.percentile(xx, 99.9, axis=0) )

      return xx, qt


  def preprocessing(self, train_samples, test_samples):
    '''
      Preprocessing of data before training/testing the NN
      This includes:
        - building the main_df
        - building the scaler
        - get the scaled features xx
        - get the target Y
    '''

    

    # concatenate all the datasets and shuffle events
    # pd which don't have a "central_w" column will just carry a NAN in this column
    train = pd.concat([train_samples[key] for key in train_samples.keys()], sort = False)
    test  = pd.concat([test_samples[key] for key  in test_samples.keys()],  sort = False)

    # since our dataset is cleaned from nans (except the one we get from concatenating),
    # we can replace them with 1.0. Like this, all non-hammered events have simply weight 1
    # same for sf_weight, it will be 1.0 for MC
    # same for trigger_sf, it will be 1.0 for data
    #pdb.set_trace()
    train = train.fillna(1.0)
    test  = test .fillna(1.0)

    # create new column which holds the multiplication of all weights
    #train["total_w"] = train["sf_weights"] * train["central_w"] * train["trigger_sf"]
    #test ["total_w"] = test ["sf_weights"] * test ["central_w"] * test ["trigger_sf"]

    train["total_w"] = train["central_w"] * train["trigger_sf"]
    test ["total_w"] = test ["central_w"] * test ["trigger_sf"]


    # even undo the splitting which is not used for k-folding
    main_df = pd.concat([train,test],sort= False)

    # re-index (this does not shuffle events, only reindex!)
    main_df.index = np.array(range(len(main_df)))

    # shuffle
    main_df = main_df.sample(frac=1, replace=False, random_state=1986) # of course, keep R's seed ;)

    # X and Y, keep event number for kfold splitting! will be removed later!

    #w_cols = ['central_w', 'event']
    #w_cols = ['sf_weights', 'central_w', 'event']
    w_cols = ['total_w', 'event']

    X       = pd.DataFrame(main_df, columns=list(set(self.features)) + ['event'] )
    Y       = pd.DataFrame(main_df, columns=['is_signal', 'event'])
    weights = pd.DataFrame(main_df, columns=w_cols)
    # save the exact list of features
    features_filename = '/'.join([self.outdir, 'input_features.pck'])
    pickle.dump(self.features, open(features_filename, 'wb' ))
    print('========> {} created'.format(features_filename))

    return X, Y, weights

  def prepareFold(self, X, Y, weights, fold):

    # create cyclic list, f.e. for nfold = 5, this is:
    # [0,1,2,3,4,0,1,2,3,4]
    cyclic = list(range(self.nfolds)) * 2  

    # define the test and train indices for this fold <fold>
    train_indices = [cyclic[fold + (i+1)] for i in range(self.nfolds - 2)  ]
    test_index    = cyclic [fold + self.nfolds - 1]

    print(f"For fold/net {fold} we train on folds {train_indices} and test on fold {test_index}")

    #########################################################
    # Prepare concat fold list for training ((n-1)/n) folds #
    #########################################################

    xx_folds_train = []
    y_folds_train  = []
    w_folds_train  = []

    for m in train_indices:

      #select training folds by % evt number 
      x_fold_train = X      [ X["event"]       % self.nfolds == m ]
      y_fold_train = Y      [ Y["event"]       % self.nfolds == m ]
      w_fold_train = weights[ weights["event"] % self.nfolds == m ]

      #delete event columns afterwards!
      x_fold_train = x_fold_train.drop("event", axis = 1)
      y_fold_train = y_fold_train.drop("event", axis = 1)
      w_fold_train = w_fold_train.drop("event", axis = 1)

      # scale the features
      # this is an important step!
      xx_fold_train, qt =  self.doScaling(x_fold_train)

      # and make the xx pandas again after scaling(the rest still ist)
      xx_fold_train     = pd.DataFrame(xx_fold_train,  columns=list(set(self.features))) # alternative to X_train[self.features[:]]

      # reset pandas index to start from 0
      # this doees not shuffle rows but jsut renumbers the index (first columns)
      xx_fold_train = xx_fold_train.reset_index(drop=True)
      y_fold_train  = y_fold_train .reset_index(drop=True)
      w_fold_train  = w_fold_train .reset_index(drop=True)

      #append to big list
      xx_folds_train.append(xx_fold_train)
      y_folds_train .append(y_fold_train)
      w_folds_train .append(w_fold_train)

      # and save_train the scaler, which will have to be used throughout the full process, even at evaluation time
      scaler_filename = '/'.join([self.outdir, f'input_tranformation_weighted_fold_{m}.pck'])
      pickle.dump(qt,open(scaler_filename, 'wb'))
      print( ' ========> {} created'.format(scaler_filename))

    #after looping over all train indices, we pd concat the list      
    xx_train = pd.concat(xx_folds_train)
    y_train  = pd.concat(y_folds_train)
    w_train  = pd.concat(w_folds_train)
 
 
    #########################################################
    # Prepare concat fold list for validation 1/n folds     #
    #########################################################

    x_val = X      [ X["event"]       % self.nfolds == test_index ]
    y_val = Y      [ Y["event"]       % self.nfolds == test_index ]
    w_val = weights[ weights["event"] % self.nfolds == test_index ]
 
    #delete event columns afterwards!
    x_val = x_val.drop("event", axis = 1)
    y_val = y_val.drop("event", axis = 1)
    w_val = w_val.drop("event", axis = 1)

    xx_val, qt =  self.doScaling(x_val)

    # and make the xx pandas again after scaling(the rest still ist)
    xx_val     = pd.DataFrame(xx_val,  columns=list(set(self.features))) # alternative to X_val[self.features[:]]

    # reset pandas index to start from 0
    # this doees not shuffle rows but jsut renumbers the index (first columns)
    xx_val = xx_val.reset_index(drop=True)
    y_val  = y_val .reset_index(drop=True)
    w_val  = w_val .reset_index(drop=True)
   
    #save fold to disk and open later again to save RAM
    #with open(f'{self.outdir}/x_fold_{n}.pck'  , 'wb') as f:
    #  pickle.dump(x_fold  ,f)
    #with open(f'{self.outdir}/y_fold_{n}.pck'  , 'wb') as f:
    #  pickle.dump(y_fold  ,f)
    #with open(f'{self.outdir}/w_fold_{n}.pck'  , 'wb') as f:
    #  pickle.dump(w_fold  ,f)
    #with open(f'{self.outdir}/xx_fold_{n}.pck'  , 'wb') as f:
    #  pickle.dump(xx_fold  ,f)

 
    return xx_train, y_train, w_train, xx_val, y_val, w_val 


  def defineModel(self):
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
    model.add(tf.keras.layers.Input((len(features),)))
    model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    model.add(tf.keras.layers.Dense(128,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    model.add(tf.keras.layers.Dense(6  ,  activation= 'softmax'                                            ))
    opt = keras.optimizers.Adam(learning_rate=learning_rate) 
    model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=['acc'])
   
    # OPTUNA SUGGESTION 
    # Trial 57 finished with value: 0.8491290502301749 and parameters: 
    #{'n_layers': 7, 'regu_rate': 0.0010307639101928131, 
    #'num_nodes_0': 246, 
    #'num_nodes_1': 300, 
    #'num_nodes_2': 44, 
    #'num_nodes_3': 275, 
    #'num_nodes_4': 265, 
    #'num_nodes_5': 299, 
    #'num_nodes_6': 55, 
    #'optimizer': 'Adam', 'adam_learning_rate': 0.00021996309401293347, 'batch_size': 40}. 

    #l2_rate = regularizers.l2(0.001)
    #learning_rate = 0.00022
    #model = tf.keras.Sequential()
    #model.add(tf.keras.layers.Input((len(features),)))
    #model.add(tf.keras.layers.Dense(246,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    #model.add(tf.keras.layers.Dense(300,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    #model.add(tf.keras.layers.Dense(44,   activation ='swish'  ,  kernel_regularizer=l2_rate))
    #model.add(tf.keras.layers.Dense(275,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    #model.add(tf.keras.layers.Dense(265,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    #model.add(tf.keras.layers.Dense(299,  activation ='swish'  ,  kernel_regularizer=l2_rate))
    #model.add(tf.keras.layers.Dense(6  ,  activation= 'softmax'                                            ))
    #opt = keras.optimizers.Adam(learning_rate=learning_rate) 
    #model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=['acc'])


 
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
    f = open(self.outdir + "/settings.txt", "a")
    with redirect_stdout(f): 
      model.summary()
      print_model_summary(model)
    f.close()
  
    return model


  def defineCallbacks(self, fold):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=3000)
    
    # reduce learning rate when at plateau, fine search the minimum
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=5, min_lr=0.000001, cooldown=10, verbose=True)
    
    # save the model every now and then
    filepath = '/'.join([self.outdir, f'fold_{fold}' + '_saved-model-{epoch:04d}_val_loss_{val_loss:.4f}_val_acc_{val_acc:.4f}.h5'])
    save_model = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=False, mode='auto', period=1)
    
    callbacks = [save_model]

    if self.do_early_stopping:
      callbacks.append(es)

    if self.do_reduce_lr:
      callbacks.append(reduce_lr)

    #callbacks.append(CustomCallback())

    return callbacks


  def prepareInputs(self, xx_fold, y_fold, w_fold):

    '''
      Note: the input xx should arlready be scaled
    '''

    # the scaled arrays within xx_folds are not pandas anymore, lets make them pandas again!
    # y_folds still contains pandas, as it has not been scaled!

    for n in range(self.nfolds): 

      xx_folds[n]  = pd.DataFrame(xx_folds[n],  columns=list(set(self.features))) # alternative to X_train[self.features[:]]

      # reset pandas index to start from 0
      # this doees not shuffle rows but jsut renumbers the index (first columns)
      xx_folds[n] = xx_folds[n].reset_index(drop=True)
      y_folds[n]  = y_folds[n].reset_index(drop=True)
      w_folds[n]  = w_folds[n].reset_index(drop=True)

    return xx_folds, y_folds, w_folds 


  def train(self, xx_train, y_train, w_train, xx_val, y_val, w_val, class_label, fold):
    '''
      Perform the training
    '''

    model    = {}
    history  = {}

    # calculate the class weight (not the same as sample weight!
    # this corrects for class imbalance including the effect of the sample
    # weights
    eff_pop = np.zeros(6)

    for i in range(0,6):

      # effective population of classes 0-5 for fold n:
      eff_pop[i] = w_train[y_train['is_signal'] == i].sum()

    total_pop = np.sum(eff_pop)
    class_w   = {i: total_pop / eff_pop[i] for i in range (0,6)} 

    #convert to numpy for training (and use one-hot for y)
    xx_train  = xx_train.to_numpy()
    xx_val    = xx_val.to_numpy()

    y_train   = tf.one_hot(y_train['is_signal'].to_numpy(), 6)
    y_train   = y_train.numpy()

    y_val     = tf.one_hot(y_val['is_signal'].to_numpy(), 6)
    y_val     = y_val.numpy()

    w_train   = w_train.to_numpy()
    w_val     = w_val.to_numpy()

    callbacks = self.defineCallbacks(fold)

    #save the output dictionaries
    #with open(f'{self.outdir}/xx_val_{n}.pck'  , 'wb') as f:
    #  pickle.dump(xx_val[n]  ,f)
    
    #with open(f'{self.outdir}/xx_train_{n}.pck', 'wb') as f:
    #  pickle.dump(xx_train[n],f)
    
    #with open(f'{self.outdir}/y_val_{n}.pck'   , 'wb') as f:
    #  pickle.dump(y_val[n]   ,f)
    
    #with open(f'{self.outdir}/y_train_{n}.pck' , 'wb') as f:
    #  pickle.dump(y_train[n] ,f)

    #now train!

    #define the model every time!! Otherwise you train the same model over and over again lol
    model = self.defineModel()

    import timeit
    start = timeit.timeit()

    print("====> Flatten weights ...")
    w_train_flat = w_train.flatten().astype(np.float32)
    w_val_flat   = w_val.flatten().astype(np.float32)
    print("====> Flatten done ...")
    end = timeit.timeit()
    print(end - start)

    print("====> Convert to tensor...")
    w_train_tf = tf.convert_to_tensor(w_train_flat, dtype=tf.float32)
    w_val_tf   = tf.convert_to_tensor(w_val_flat, dtype=tf.float32)
    print("====> Conversion done ...")


    with open(f"{scratch_dir}/xx_train_{fold}.pck", "wb") as f:
      pickle.dump(xx_train, f)
    with open(f"{scratch_dir}/y_train_{fold}.pck", "wb") as f:
      pickle.dump(y_train, f)
    with open(f"{scratch_dir}/w_train_{fold}.pck", "wb") as f:
      pickle.dump(w_train, f)

    with open(f"{scratch_dir}/xx_val_{fold}.pck", "wb") as f:
      pickle.dump(xx_val, f)
    with open(f"{scratch_dir}/y_val_{fold}.pck", "wb") as f:
      pickle.dump(y_val, f)
    with open(f"{scratch_dir}/w_val_{fold}.pck", "wb") as f:
      pickle.dump(w_val, f)

    with open(f"{scratch_dir}/class_w.pck", "wb") as f:
      pickle.dump(class_w, f)
 
    with open(f"{scratch_dir}/len_features.pck", "wb") as f:
      pickle.dump(len(features), f)

    #sys.exit()

    # history = model.fit(xx_train, y_train, validation_data=(xx_val, y_val, w_val_tf), epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True,class_weight = class_w , sample_weight = w_train_tf )

    # remove sample weight for fast debugging!
    #history = model.fit(xx_train, y_train, validation_data=(xx_val, y_val), epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True,class_weight = class_w  )

    # remove class weights
    #history = model.fit(xx_train, y_train, validation_data=(xx_val, y_val), epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True, sample_weight = pd.Series(w_train_flat))
 
    #epochs, toPlot_loss     = self.plotMetric(history, "loss"      ,"Training Loss"       , "Loss"   , fold = fold)
    #epochs, toPlot_acc      = self.plotMetric(history, "acc"       ,"Training Accuracy"   , "Acc."   , fold = fold)
    #epochs, toPlot_val_loss = self.plotMetric(history, "val_loss"  ,"Validation Loss"     , "Loss"   , fold = fold)
    #epochs, toPlot_val_acc  = self.plotMetric(history, "val_acc"   ,"Validation Accuracy" , "Acc."   , fold = fold)

    #h1_0,h2_0 = self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 0, fold = fold)
    #h1_1,h2_1 = self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 1, fold = fold)
    #h1_2,h2_2 = self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 2, fold = fold)
    #h1_3,h2_3 = self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 3, fold = fold)
    #h1_4,h2_4 = self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 4, fold = fold)
    #h1_5,h2_5 = self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 5, fold = fold)

    #score_0 = self.plotScore (model, xx_train, y_train, xx_val, y_val, 0, class_label, fold = fold)
    #score_1 = self.plotScore (model, xx_train, y_train, xx_val, y_val, 1, class_label, fold = fold)
    #score_2 = self.plotScore (model, xx_train, y_train, xx_val, y_val, 2, class_label, fold = fold)
    #score_3 = self.plotScore (model, xx_train, y_train, xx_val, y_val, 3, class_label, fold = fold)
    #score_4 = self.plotScore (model, xx_train, y_train, xx_val, y_val, 4, class_label, fold = fold)
    #score_5 = self.plotScore (model, xx_train, y_train, xx_val, y_val, 5, class_label, fold = fold)
 
    #return epochs,          \
    #       toPlot_loss,     \
    #       toPlot_acc,      \
    #       toPlot_val_loss, \
    #       toPlot_val_acc,  \
    #       score_0,         \
    #       score_1,         \
    #       score_2,         \
    #       score_3,         \
    #       score_4,         \
    #       score_5,         \
    #       h1_0, h2_0,      \
    #       h1_1, h2_1,      \
    #       h1_2, h2_2,      \
    #       h1_3, h2_3,      \
    #       h1_4, h2_4,      \
    #       h1_5, h2_5,      \


  def plotMetric(self, history, history_key, title, ylabel, fold = None):
    '''
      Plot the loss/acc metric for training and validation sets
    '''

    # folds can have different length due to early stopping!!!

    # TOTAL LOSS

    toPlot = history.history[history_key]
    epochs = range(1, len(toPlot)+1)
    h = plt.plot(epochs, toPlot, self.colors[fold], label=f'Fold {fold}')

    plt.subplots_adjust(left=0.2, right=0.95, top=0.85, bottom=0.20)

    plt.title(title)
    plt.xlabel('Epochs')
    plt.ylabel(ylabel)
    plt.legend()
    self.saveFig(plt, history_key.replace(" ", "_") + f"_fold{fold}")
    self.saveFig(plt, history_key.replace(" ", "_") + f"_fold{fold}")

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
    self.saveFig(plt, history_key.replace(" ", "_") + "_log" + f"_fold{fold}")
    self.saveFig(plt, history_key.replace(" ", "_") + "_log" + f"_fold{fold}")
    #plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))
    plt.clf()
    plt.close()

    return epochs, toPlot

  def plotAverageMetric(self, epochs, toPlot, history_key, title, ylabel):

    # TOTAL LOSS
    #pdb.set_trace() 
    average = []
    #max epochs (get max epochs)
    max_epochs = max([len(val) for val in epochs.values()])

    for n in range(self.nfolds):

      plt.plot(epochs[n], toPlot[n], self.colors[n], label=f'Fold {n}')

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
    self.saveFig(plt, history_key + "_average")
    self.saveFig(plt, history_key + "_average")
  
    plt.close() 


  #  #log
  #  ymin = 0
  #  ymax = np.max(average_toPlot) 
  #  plt.yscale('log')  
  #  ax = plt.gca()
  #  ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
  #  ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=20))
  #  ax.yaxis.set_major_formatter(LogFormatter(base=10.0, labelOnlyBase=False))
  #  plt.ylim(bottom=ymin * 0.9, top=ymax * 1.1)    

  #  self.saveFig(plt, history_key.replace(" ", "_") + "_log")
  #  self.saveFig(plt, history_key.replace(" ", "_") + "_log")
  #  plt.clf()
  #  plt.close()


  def plotLoss(self, history):
    '''
      Plot the loss for training and validation sets
    '''

    # folds can have different length due to early stopping!!!

    #training
    average_loss = []
    #max epochs (get max epochs)
    max_epochs = max(len(history[n].history['loss']) for n in range(self.nfolds))

    for n in range(self.nfolds):

      loss_train = history[n].history['loss']
      #average_loss.append(np.array(loss_train))
      epochs = range(1, len(loss_train)+1)
      plt.plot(epochs, loss_train, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the loss (the rest stays nan)
      padded[:len(loss_train)] = loss_train
      #append the padded array to the average loss
      average_loss.append(padded)

    #take the average loss and ignore nans
    average_loss = np.nanmean(np.vstack(average_loss), axis=0)
    plt.plot(range(1, max_epochs + 1), average_loss, 'black' , label='Average loss', linestyle= 'dashed')

    plt.title('Training Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'training_loss')
    plt.clf()
    plt.close()

    #log
 
    average_loss = []
    for n in range(self.nfolds):

      loss_train = history[n].history['loss']
      #average_loss.append(np.array(loss_train))
      epochs = range(1, len(loss_train)+1)
      plt.plot(epochs, loss_train, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the loss (the rest stays nan)
      padded[:len(loss_train)] = loss_train
      #append the padded array to the average loss
      average_loss.append(padded)

    #import pdb; pdb.set_trace();
    average_loss = np.nanmean(np.vstack(average_loss), axis=0)
    plt.plot(range(1, max_epochs + 1), average_loss, 'black' , label='Average loss', linestyle= 'dashed')
    plt.title('Training Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.yscale('log')  
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
    plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))
    plt.legend()
    self.saveFig(plt, 'training_loss_log')
    plt.clf()
    plt.close()

    #validation
    average_loss = []
    for n in range(self.nfolds):

      loss_val = history[n].history['val_loss']
      epochs = range(1, len(loss_val)+1)
      plt.plot(epochs, loss_val, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the loss (the rest stays nan)
      padded[:len(loss_val)] = loss_val
      #append the padded array to the average loss
      average_loss.append(padded)


    average_loss = np.nanmean(np.vstack(average_loss), axis=0)
    plt.plot(range(1, max_epochs + 1), average_loss, 'black' , label='Average loss', linestyle= 'dashed')
    plt.title('Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'validation_loss')
    plt.clf()
    plt.close()

    #log
    average_loss = []
    for n in range(self.nfolds):

      loss_val = history[n].history['val_loss']
      epochs = range(1, len(loss_val)+1)
      plt.plot(epochs, loss_val, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the loss (the rest stays nan)
      padded[:len(loss_val)] = loss_val
      #append the padded array to the average loss
      average_loss.append(padded)


    average_loss = np.nanmean(np.vstack(average_loss), axis=0)
    plt.plot(range(1, max_epochs + 1), average_loss, 'black' , label='Average loss', linestyle= 'dashed')
    plt.title('Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.yscale('log')  
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
    plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))
    plt.legend()
    self.saveFig(plt, 'validation_loss_log')
    plt.clf()

    plt.close()



  def plotAccuracy(self, history):
    '''
      Plot the accuracy for training and validation sets
    '''

    #training
    average_acc = []

    #max epochs (get max epochs)
    max_epochs = max(len(history[n].history['acc']) for n in range(self.nfolds))

    for n in range(self.nfolds):

      acc_train = history[n].history['acc']
      #average_acc.append(np.array(acc_train))
      epochs = range(1, len(acc_train)+1)
      plt.plot(epochs, acc_train, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the acc (the rest stays nan)
      padded[:len(acc_train)] = acc_train
      #append the padded array to the average acc
      average_acc.append(padded)


    average_acc = np.nanmean(np.vstack(average_acc), axis=0)
    plt.plot(range(1, max_epochs + 1), average_acc, 'black' , label='Average acc', linestyle= 'dashed')

    plt.title('Training Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend()
    self.saveFig(plt, 'training_acc')
    plt.clf()
    plt.close()

    #log
    average_acc = []
    for n in range(self.nfolds):

      acc_train = history[n].history['acc']
      #average_acc.append(np.array(acc_train))
      epochs = range(1, len(acc_train)+1)
      plt.plot(epochs, acc_train, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the acc (the rest stays nan)
      padded[:len(acc_train)] = acc_train
      #append the padded array to the average acc
      average_acc.append(padded)

    average_acc = np.nanmean(np.vstack(average_acc), axis=0)
    plt.plot(range(1, max_epochs + 1), average_acc, 'black' , label='Average acc', linestyle= 'dashed')
    plt.title('Training Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.yscale('log')  
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
    plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))

    plt.legend()
    self.saveFig(plt, 'training_acc_log')
    plt.clf()
    plt.close()

    #validation
    average_acc = []
    for n in range(self.nfolds):

      acc_val = history[n].history['val_acc']
      #average_acc.append(np.array(acc_val))
      epochs = range(1, len(acc_val)+1)
      plt.plot(epochs, acc_val, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the acc (the rest stays nan)
      padded[:len(acc_val)] = acc_val
      #append the padded array to the average acc
      average_acc.append(padded)


    average_acc = np.nanmean(np.vstack(average_acc), axis=0)
    plt.plot(range(1, max_epochs + 1), average_acc, 'black' , label='Average acc', linestyle= 'dashed')
    plt.title('Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend()
    self.saveFig(plt, 'validation_acc')
    plt.clf()
    plt.close()

    #log
    average_acc = []
    for n in range(self.nfolds):

      acc_val = history[n].history['val_acc']
      #average_acc.append(np.array(acc_val))
      epochs = range(1, len(acc_val)+1)
      plt.plot(epochs, acc_val, self.colors[n], label=f'Fold {n}')

      #create a full array of nans of length max_epochs       
      padded = np.full(max_epochs, np.nan)
      #fill the first part with the acc (the rest stays nan)
      padded[:len(acc_val)] = acc_val
      #append the padded array to the average acc
      average_acc.append(padded)

    average_acc = np.nanmean(np.vstack(average_acc), axis=0)
    plt.plot(range(1, max_epochs + 1), average_acc, 'black' , label='Average acc', linestyle= 'dashed')
    plt.title('Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.yscale('log')  
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
    plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))

    plt.legend()
    self.saveFig(plt, 'validation_acc_log')
    plt.clf()
    plt.close()




  def predictScore(self, model, df):
    '''
      Return score with scaled input features
    '''
    x = pd.DataFrame(df, columns=self.features)

    # apply the scaler
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted.pck'])
    qt = pickle.load(open(scaler_filename, 'rb'))
    xx = qt.transform(x[self.features])
    #apply same scaling in this function
    # predict
    score = model.predict(xx)

    return score

  def plotScoreTauOnly(self, model, xx_train, y_train, xx_val, y_val, class_label):
    '''
      Plot the score of all channels class sig
    '''

    sig = 1 #only interested in tau star

    channels = [0,1]
    col ={0:'c',1:'g',2:'m',3:'r',4:'b',5:'k'}
    linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}

    hist_content = {}

    for chan in channels:

      dummy = []

      for n in range(self.nfolds):
  
        #class predictions 1D
        #is of shape [1,2,5,4,3,2,4,2,1, ....] 
        true_1d = np.argmax(y_val[n], axis=1)

        #select only events where the true class is chan...
        x_chan    = xx_val[n][ true_1d == chan ]
        #...and predict their score! (data is already scaled!)
        y_chan    = model[n].predict(x_chan)[:,sig] 

        # Plot data into histo
        hist = plt.hist(y_chan, bins=np.arange(0,1.025,0.025), color=col[chan], alpha=0.5, label=class_label[chan], histtype='stepfilled',density = True,linestyle = linestyles[chan], linewidth = 1.5)
  
        #append the bin content (at index 0) for every fold! 
        dummy.append(hist[0])

      #average over the folds to get the average histo
      hist_content[chan] = sum(dummy)/self.nfolds


    fig = plt.figure()

    for chan in channels:

      #bin edges are at index 1 (the same for all, can take the last one)
      bin_edges = hist[1]
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
    self.saveFig(fig, f'tau_only' )
    plt.clf()
    plt.close()


  def plotScore(self, model, xx_train, y_train, xx_val, y_val, sig, class_label, fold = None):
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
    self.saveFig(fig, f'score_' + str(sig) + f"_fold_{fold}" )
    plt.clf()
    plt.close()
 
    return hists
 
  def plotAverageScore(self, scores, sig, class_label):

    channels = [0,1,2,3,4,5]
    col ={0:'c',1:'g',2:'m',3:'r',4:'b',5:'k'}
    linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}

    hist_content = {}
    #plot average
    for chan in channels:
      dummy = []
    
      for n in range(self.nfolds):
     
        #append the bin content (at index 0) for every fold! 
        dummy.append(scores[n][chan][0])
    
      #average over the folds to get the average histo
      hist_content[chan] = sum(dummy)/self.nfolds
    
    fig = plt.figure()
    #pdb.set_trace() 
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
    self.saveFig(fig, "score_" + str(sig) + "_average" )
    plt.clf()
    plt.close()
  
  
  # COMMENT THIS OUT   
  #def plotScore(self, model, xx_train, y_train, xx_val, y_val, sig, class_label, fold = None):
  #  '''
  #    Plot the score of all channels class sig
  #  '''

  #  channels = [0,1,2,3,4,5]
  #  col ={0:'c',1:'g',2:'m',3:'r',4:'b',5:'k'}
  #  linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}

  #  hist_content = {}

  #  for chan in channels:


  #    if fold == None:
  #      #no explicit fold, take fold average
 
  #      dummy = []
  #
  #      for n in range(self.nfolds):
  #     
  #        print("chan is ", chan, ", n is ", n, ", sig is ", sig) 
  #        #class predictions 1D
  #        #is of shape [1,2,5,4,3,2,4,2,1, ....] 
  #        true_1d = np.argmax(y_val[n], axis=1)
  #
  #        #select only events where the true class is chan...
  #        x_chan    = xx_val[n][ true_1d == chan ]
  #        #...and predict their score! (data is already scaled!)
  #        y_chan    = model[n].predict(x_chan)[:,sig] 
  #
  #        np.savetxt(f"{self.outdir}/fold_{n}_score_of_class_{chan}_for_class_{sig}.csv",y_chan,delimiter = ",")
  #
  #        # Plot data into histo
  #        hist = plt.hist(y_chan, bins=np.arange(0,1.025,0.025), color=col[chan], alpha=0.5, label=class_label[chan], histtype='stepfilled',density = True,linestyle = linestyles[chan], linewidth = 1.5)
  #  
  #        #append the bin content (at index 0) for every fold! 
  #        dummy.append(hist[0])
  #
  #      #average over the folds to get the average histo
  #      hist_content[chan] = sum(dummy)/self.nfolds

  #    else:
  #    
  #      #take fold average
  #      n = fold
  #      #is of shape [1,2,5,4,3,2,4,2,1, ....] 
  #      true_1d = np.argmax(y_val[n], axis=1)
  #      
  #      #select only events where the true class is chan...
  #      x_chan    = xx_val[n][ true_1d == chan ]
  #      #...and predict their score! (data is already scaled!)
  #      y_chan    = model[n].predict(x_chan)[:,sig] 
  #      
  #      np.savetxt(f"{self.outdir}/fold_{n}_score_of_class_{chan}_for_class_{sig}.csv",y_chan,delimiter = ",")

  #      # Plot data into histo
  #      hist = plt.hist(y_chan, bins=np.arange(0,0.425,0.025), color=col[chan], alpha=0.5, label=class_label[chan], histtype='stepfilled',density = True,linestyle = linestyles[chan], linewidth = 1.5)
  #
  #        hist_content[chan] = hist[0] 
  #
  #
  #    #fig = plt.figure()
  #
  #    #for chan in channels:
  #
  #    #  #bin edges are at index 1 (the same for all, can take the last one)
  #    #  bin_edges = hist[1]
  #    #  bin_width = np.diff(bin_edges) 
  # 
  #    #  # plot the score distributions
  #    #  plt.bar(bin_edges[:-1], hist_content[chan], color=col[chan], width = bin_width, align = 'edge', alpha=0.5, label=class_label[chan], linestyle = linestyles[chan], linewidth = 1.5)
  #  
  #    ##plot all channels into same plot 
  #    #plt.legend(loc='upper right')
  #    #plt.title(f'Score for class ' + class_label[sig] + " (Average over folds)")
  #    #plt.xlabel('score')
  #    #plt.ylabel('events')
  #    ##fig.savefig('outputs/score_' + str(sig) + '.pdf')
  #    ##fig.savefig('outputs/score_' + str(sig) + '.png')
  #    #self.saveFig(fig, f'score_' + str(sig) )
  #    #plt.clf()
  #    #plt.close()
  #
  #
  #    #fig = plt.figure()
  #
  #    #for chan in channels:
  #
  #    #  #bin edges are at index 1 (the same for all, can take the last one)
  #    #  bin_edges = hist[1]
  #    #  bin_width = np.diff(bin_edges) 
  # 
  #    #  # plot the score distributions
  #    #  plt.bar(bin_edges[:-1], hist_content[chan], color=col[chan], width = bin_width, align = 'edge', alpha=0.5, label=class_label[chan], linestyle = linestyles[chan], linewidth = 1.5)
  #  
  #    ##plot all channels into same plot 
  #    #plt.legend(loc='upper right')
  #    #plt.title(f'Score for class ' + class_label[sig] + " (Average over folds)")
  #    #plt.xlabel('score')
  #    #plt.ylabel('events')
  #
  #    ## log scale
  #    #plt.yscale('log')
  #    #self.saveFig(fig, f'score_' + str(sig) + '_log')
  #    #plt.clf()
  #    #plt.close()



  def plotCorr(self, model, xx_train, y_train, xx_val, y_val, class_label, fold):

    #get the score
    score  =  model[fold].predict(xx_val[fold]) #this is a list of length #classes!

    test_x = xx_val[fold].copy()
    test_x = pd.DataFrame(test_x, columns=list(set(self.features)))

    #append it to pd such that we can plot it in the corr matrix
    test_x["score0"] = score[:,0]
    test_x["score1"] = score[:,1]
    test_x["score2"] = score[:,2]
    test_x["score3"] = score[:,3]
    test_x["score4"] = score[:,4]
    test_x["score5"] = score[:,5]

    corr = test_x.corr()


    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))
    
    # Generate a custom diverging colormap
    cmap = seaborn.diverging_palette(220, 10, as_cmap=True)
    
    # Draw the heatmap with the mask and correct aspect ratio
    g = seaborn.heatmap(corr, cmap=cmap, vmax=1., vmin=-1, center=0, annot=True, fmt='.1f', square=True, linewidths=.8, cbar_kws={"shrink": .8},  annot_kws={"size": 2})
    
    # force all ticklabels to show (otherwise it gives an error with set_xticklabels(..)
    g.set_xticks(range(len(test_x.keys())))
    g.set_yticks(range(len(test_x.keys())))

    # rotate axis labels
    g.set_xticklabels(test_x.keys().tolist(), rotation='vertical')
    g.set_yticklabels(test_x.keys().tolist(), rotation='horizontal')
    
    # plt.show()
    #plt.figure(figsize=(12, 12))
    plt.title(f'Linear correlation matrix - Fold {fold}')
    plt.tight_layout()
    self.saveFig(plt, f'corr_sig_fold_{fold}')
    plt.clf()
    plt.close()

  def plotScoreOneVsAll(self, model, xx_train, y_train, xx_val, y_val, sig, class_label):
    '''
      Plot the score of all channels class sig
    '''

    channels = [0,1,2,3,4,5]
    col ={0:ROOT.kCyan,1:ROOT.kGreen,2:ROOT.kMagenta,3:ROOT.kRed,4:ROOT.kBlue,5:ROOT.kGray}
    linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}

    hist_content_val = {}
    hist_content_train = {}

    for chan in channels:

      dummy_val = []
      dummy_train = []

      for n in range(self.nfolds):
  
        #class predictions 1D
        #is of shape [1,2,5,4,3,2,4,2,1, ....] 
        true_1d_val   = np.argmax(y_val[n], axis=1)
        true_1d_train = np.argmax(y_train[n], axis=1)

        #select only events where the true class is chan...
        x_chan_val    = xx_val[n][ true_1d_val == chan ]
        x_chan_train  = xx_train[n][ true_1d_train == chan ]
        #...and predict their score! (data is already scaled!)
        y_chan_val    = model[n].predict(x_chan_val)[:,sig] 
        y_chan_train  = model[n].predict(x_chan_train)[:,sig] 


        # Plot data into histo
        hist_val   = plt.hist(y_chan_val, bins=np.arange(0,1.025,0.025))
        hist_train = plt.hist(y_chan_train, bins=np.arange(0,1.025,0.025))
  
        #append the bin content (at index 0) for every fold! 
        dummy_val.append(hist_val[0])
        dummy_train.append(hist_train[0])

      #average over the folds to get the average histo
      hist_content_val[chan]   = sum(dummy_val)/self.nfolds
      hist_content_train[chan] = sum(dummy_train)/self.nfolds


    #holds all histos excect the one which has chan = sig
    the_rest_val   = [0] * len(hist_content_val[0]) 
    the_rest_train = [0] * len(hist_content_train[0])

    fig = plt.figure()

    for chan in channels:

      if chan != sig:
        the_rest_val   = [a+b for a, b in zip(the_rest_val,   hist_content_val[chan])]
        the_rest_train = [a+b for a, b in zip(the_rest_train, hist_content_train[chan])]


    #bin edges are at index 1 (the same for all, can take the last one)
    bin_edges = hist_val[1]
    bin_width = np.diff(bin_edges) 
    n_bins = len(bin_edges) - 1 


    sig_val = ROOT.TH1D("","",n_bins, bin_edges)
    sig_train = ROOT.TH1D("","",n_bins, bin_edges)
    back_val = ROOT.TH1D("","",n_bins, bin_edges)
    back_train = ROOT.TH1D("","",n_bins, bin_edges)

    #fill root histo for plotting and KS
    for i in range(n_bins):
      sig_val.SetBinContent(i,hist_content_val[sig][i])
      sig_train.SetBinContent(i,hist_content_train[sig][i])
      back_val.SetBinContent(i,the_rest_val[i])
      back_train.SetBinContent(i,the_rest_train[i])

    c1=ROOT.TCanvas()
    if sig_val.Integral()!=0: sig_val.Scale(1./sig_val.Integral())
    if sig_train.Integral()!=0: sig_train.Scale(1./sig_train.Integral())
    if back_val.Integral()!=0: back_val.Scale(1./back_val.Integral())
    if back_train.Integral()!=0: back_train.Scale(1./back_train.Integral())

    c1.Draw()

    #label axes
    sig_train.GetXaxis().SetTitle("score")
    sig_train.GetYaxis().SetTitle("a.u.")
    sig_train.GetYaxis().SetRangeUser(0, max([sig_val.GetMaximum(),sig_train.GetMaximum(),back_val.GetMaximum(),back_train.GetMaximum()])*1.6)


    #signal channel is like given
    sig_train.SetFillColor(col[sig])
    sig_train.SetLineColor(col[sig])
    sig_train.SetFillStyle(3345)
    sig_train.Draw('HIST')
 
    sig_val.SetLineColor(col[sig])
    sig_val.SetMarkerColor(col[sig])
    sig_val.SetMarkerStyle(8)
    sig_val.Draw('EP SAME')
 
    back_train.SetFillColor(ROOT.kBlack)
    back_train.SetLineColor(ROOT.kBlack)
    back_train.SetFillStyle(3345)
    back_train.Draw('HIST SAME')
 
    back_val.SetLineColor(ROOT.kBlack)
    back_val.SetMarkerColor(ROOT.kBlack)
    back_val.SetMarkerStyle(8)
    back_val.Draw('EP SAME')
 
    ks_score_sig  = sig_val.KolmogorovTest(sig_train)
    print(ks_score_sig)
    ks_score_back = back_val.KolmogorovTest(back_train)
    print(ks_score_back)
    ks_value = ROOT.TPaveText(0.25, 0.76, 0.88, 0.80, 'nbNDC')
    ks_value.AddText(f'KS score of class {sig} = {round(ks_score_sig,3)} vs KS score of the rest = {round(ks_score_back,3)}')
    ks_value.SetFillColor(0)
    ks_value.Draw('EP SAME')

    leg = ROOT.TLegend(.11,.82,.83,.88)
    leg.AddEntry(sig_train ,f'Train Class {sig}' ,'F' )
    leg.AddEntry(sig_val ,f'Test Class {sig}' ,'EP' )
    leg.AddEntry(back_train ,f'Train - the rest' ,'F' )
    leg.AddEntry(back_val ,f'Test - the rest' ,'EP' )
    leg.SetNColumns(4)
    leg.Draw("SAME")


    c1.SaveAs(self.outdir + f'/one_vs_all_KS_test_{sig}.pdf')
    c1.SaveAs(self.outdir + f'/one_vs_all_KS_test_{sig}.png')
    #print('KS score: ',ks_score, len(train_pred),len(test_pred))

  def plotCM(self, model, test_df, class_label):
    '''
      Get and save everything rto later produce confusion matrix
    '''

    test_x = pd.DataFrame(test_df, columns=list(set(self.features)))
    test_y = pd.DataFrame(test_df, columns=['is_signal'])
    test_y = test_y.to_numpy().flatten()
    
    #get the score
    score  =  self.predictScore(model,test_x) #this is a list of length 6!
    #save it
    np.savetxt(f"{self.outdir}/scoreforCM.csv",score,delimiter = ",")

    score_list = []
    for i,event in enumerate(score):
      #get class with highest score
      score_list.append(np.argmax(event))
    score_list = np.array(score_list) 

    #save it
    np.savetxt(f"{self.outdir}/scorelistforCM.csv",score_list,delimiter = ",")
    np.savetxt(f"{self.outdir}/testyforCM.csv",test_y,delimiter = ",")

    cm = confusion_matrix(test_y, score_list)   
    print("Confusion matrix: \n", cm)


    cm_recall = []
    cm_presi= []

    for i,row in enumerate(cm):
    	cm_recall.append(row / np.sum(row ))
    	cm_presi.append(cm[:,i]/np.sum(cm[:,i]))
    print("Confusion matrix recall: \n " , cm_recall)	
    print("Confusion matrix presi: \n " , cm_presi)

    #get labels
    fig = plt.figure()
    sb = seaborn.heatmap(cm_recall, annot = True)

    sb.set_xticklabels(list(class_label.values()))
    sb.set_yticklabels(list(class_label.values()))

    #plt.title(r'Class 0 = $D^{(*)}_\mathrm{s}\mu$, Class 1 = $D^{(*)}_\mathrm{s}\tau$, Class 2 = $H_\mathrm{b}$, Class 3 = $B_\mathrm{c}$')
    plt.ylabel('Actual')
    plt.xlabel('Predicted') 
    fig.savefig(f'{self.outdir}/cm_recall.pdf')
    fig.savefig(f'{self.outdir}/cm_recall.png')
	
    plt.clf()
    plt.close()

    fig = plt.figure()
    sb = seaborn.heatmap(cm_presi, annot = True)

    sb.set_xticklabels(list(class_label.values()))
    sb.set_yticklabels(list(class_label.values()))

    #plt.title(r'Class 0 = $D^{(*)}_\mathrm{s}\mu$, Class 1 = $D^{(*)}_\mathrm{s}\tau$, Class 2 = $H_\mathrm{b}$, Class 3 = $B_\mathrm{c}$')
    plt.ylabel('Actual')
    plt.xlabel('Predicted') 
    fig.savefig(f'{self.outdir}/cm_presi.pdf')
    fig.savefig(f'{self.outdir}/cm_presi.png')
	
    plt.clf()
    plt.close()


  def plotROCbinary(self,model,xx_train,y_train,xx_val,y_val,key):

    'plot one vs all roc curve'

    #plot train or test roc curve
    if key == 'Train':
      x = xx_train
      y = y_train
    else:
      x = xx_val
      y = y_val

    label_binarizer = {}
    y_score         = {}
    y_onehot_test   = {}

    #for one-vs-all roc curve
    for n in range(self.nfolds):

      label_binarizer[n] = LabelBinarizer().fit(y[n])
      y_score[n] = model[n].predict(x[n])
      y_onehot_test[n] = label_binarizer[n].transform(y[n])
     
    col =['c','g','m','r','b','k']
    linestyles = ['solid','dotted','dashed', 'dashdot',(0, (1, 10)),(0,(3, 5, 1, 5, 1, 5))]

    #plot roc for every class 
    for sig in range(6):
      #and plot all folds in one plot
      plt.figure()
      for n in range(self.nfolds):


        #plot roc for every class 
        class_id = np.flatnonzero(label_binarizer[n].classes_ == sig)[0]
        fpr,tpr,_ = roc_curve(y_onehot_test[n][:, class_id],y_score[n][:, class_id])    
        #save it
        np.savetxt(f"{self.outdir}/tprfpr_{sig}_{key}_fold_{n}.csv",np.array([fpr,tpr]),delimiter = ",")
        roc_auc = auc(fpr,tpr)
        plt.plot(fpr, tpr, label=f'Fold {n}, AUC = {round(roc_auc,2)}',color = self.colors[n] ,linestyle = linestyles[sig])

      plt.title(f'{key} ROC of class {sig}')
      plt.xlabel('False positive rate')
      plt.ylabel('True positive rate')
      xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
      plt.plot(xy, xy, color='grey', linestyle='--',label = 'No skill')
      plt.yscale('linear')
      plt.legend()
      self.saveFig(plt, f'ROC_{key}_class_{sig}')
      plt.clf()
      plt.close()

      #log scale
      plt.figure()
      for n in range(self.nfolds):


        #plot roc for every class 
        class_id = np.flatnonzero(label_binarizer[n].classes_ == sig)[0]
        fpr,tpr,_ = roc_curve(y_onehot_test[n][:, class_id],y_score[n][:, class_id])    
        #save it
        np.savetxt(f"{self.outdir}/tprfpr_{sig}_{key}_fold_{n}.csv",np.array([fpr,tpr]),delimiter = ",")
        roc_auc = auc(fpr,tpr)
        plt.plot(fpr, tpr, label=f'Fold {n}, AUC = {round(roc_auc,2)}',color = self.colors[n] ,linestyle = linestyles[sig])


      plt.title(f'{key} ROC of class {sig}')
      plt.xlabel('False positive rate')
      plt.ylabel('True positive rate')
      xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
      plt.plot(xy, xy, color='grey', linestyle='--', label='No skill')
      plt.yscale('log')  
      plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0]))
      plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))

      plt.legend()
      self.saveFig(plt, f'ROC_{key}_class_{sig}_log')
      plt.clf()
      plt.close()
  
  def plotKSTest(self, model, xx_train, y_train, xx_val, y_val,sig, fold):
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

 
    np.savetxt(f"{self.outdir}/scoretrain_{sig}.csv",y_chan_train,delimiter = ",")
    np.savetxt(f"{self.outdir}/scoretest_{sig}.csv",y_chan_val,delimiter = ",")
  
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


    c1.SaveAs(self.outdir + f'/KS_test_{sig}_fold_{fold}.pdf')
    c1.SaveAs(self.outdir + f'/KS_test_{sig}_fold_{fold}.png')
    #print('KS score: ',ks_score, len(train_pred),len(test_pred))
    c1.Close()
    del c1

    return h1, h2

  def plotAverageKSTest(self, h1, h2, sig):

    h1_av = ROOT.TH1F(f'train_{sig}', f'train_{sig}', 30, 0, 1)
    h2_av = ROOT.TH1F(f'val_{sig}', f'val_{sig}', 30, 0, 1)


    for n in range(self.nfolds):
  
      h1_av.Add(h1[n])
      h2_av.Add(h2[n])

    h1_av.Scale(1/self.nfolds) 
    h2_av.Scale(1/self.nfolds) 
  
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
  
  
    c1.SaveAs(self.outdir + f'/KS_test_{sig}_average.pdf')
    c1.SaveAs(self.outdir + f'/KS_test_{sig}_average.png')
    #print('KS score: ',ks_score, len(train_pred),len(test_pred))
    c1.Close()
    del c1
  

  def process(self):
    print( '------------------- MVA Trainer --------------------')
    
    print('\n========> save settings')
    self.saveSettings()

    # get the samples (dictionaries with signal key)
    print('\n========> getting the samples')
    train_samples, test_samples, class_label = self.getSamples()


    # create dataframes
    print('\n========> creating the train dataframes')
    train_notnan, weight      = self.createDataframe(train_samples)
    print('\n========> creating the test dataframes')
    test_notnan , _           = self.createDataframe(test_samples)


    # preprocessing the dataframes
    print('\n========> preprocessing and defining folds' )
    X, Y, weights = self.preprocessing(train_notnan, test_notnan)

    ###################
    # K-folding loop  #
    ###################

    epochs          = {}
    toPlot_loss     = {} 
    toPlot_acc      = {}
    toPlot_val_loss = {}
    toPlot_val_acc  = {}
    score0          = {}
    score1          = {}
    score2          = {}
    score3          = {}
    score4          = {}
    score5          = {}
    h1_0            = {}
    h2_0            = {}
    h1_1            = {}
    h2_1            = {}
    h1_2            = {}
    h2_2            = {}
    h1_3            = {}
    h2_3            = {}
    h1_4            = {}
    h2_4            = {}
    h1_5            = {}
    h2_5            = {}

    for n in range(self.nfolds):

      # prepare training and testing set for fold n
      xx_train, y_train, w_train, xx_val, y_val, w_val = self.prepareFold(X,Y,weights, n)

      # do the training for fold n
      #print('\n========> training...') 
      #epochs         [n], \
      #toPlot_loss    [n], \
      #toPlot_acc     [n], \
      #toPlot_val_loss[n], \
      #toPlot_val_acc [n], \
      #score0         [n], \
      #score1         [n], \
      #score2         [n], \
      #score3         [n], \
      #score4         [n], \
      #score5         [n], \
      #h1_0           [n], \
      #h2_0           [n], \
      #h1_1           [n], \
      #h2_1           [n], \
      #h1_2           [n], \
      #h2_2           [n], \
      #h1_3           [n], \
      #h2_3           [n], \
      #h1_4           [n], \
      #h2_4           [n], \
      #h1_5           [n], \
      #h2_5           [n]  \
      #= self.train(xx_train, y_train, w_train, xx_val, y_val, w_val, class_label, n)

      self.train(xx_train, y_train, w_train, xx_val, y_val, w_val, class_label, n)


    #pdb.set_trace()
   
    # plot fold averages
    #self.plotAverageMetric(epochs, toPlot_loss    , "loss"     ,"Training Loss"       , "Loss")
    #self.plotAverageMetric(epochs, toPlot_acc     , "acc"      ,"Training Accuracy"   , "Acc.")
    #self.plotAverageMetric(epochs, toPlot_val_loss, "val_loss" ,"Validation Loss"     , "Loss")
    #self.plotAverageMetric(epochs, toPlot_val_acc , "val_acc"  ,"Validation Accuracy" , "Acc.")

    #self.plotAverageScore (score0, 0, class_label)
    #self.plotAverageScore (score1, 1, class_label)
    #self.plotAverageScore (score2, 2, class_label)
    #self.plotAverageScore (score3, 3, class_label)
    #self.plotAverageScore (score4, 4, class_label)
    #self.plotAverageScore (score5, 5, class_label)

    #self.plotAverageKSTest (h1_0, h2_0, 0)
    #self.plotAverageKSTest (h1_1, h2_1, 1)
    #self.plotAverageKSTest (h1_2, h2_2, 2)
    #self.plotAverageKSTest (h1_3, h2_3, 3)
    #self.plotAverageKSTest (h1_4, h2_4, 4)
    #self.plotAverageKSTest (h1_5, h2_5, 5)

    # plotting
    #print('\n========> plotting...' )
    #self.plotLoss(history)
    #self.plotAccuracy(history)
    #self.plotScoreTauOnly(model, xx_train, y_train, xx_val, y_val,class_label )
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 0,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 1,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 2,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 3,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 4,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 5,class_label)
    ##self.plotScoreOneVsAll(model, xx_train, y_train, xx_val, y_val, 1,class_label)

    #####self.plotCM(model, main_test_df, class_label)
    #self.plotROCbinary(model,xx_train,y_train,xx_val,y_val,'Train')
    #self.plotROCbinary(model,xx_train,y_train,xx_val,y_val,'Test')
    #self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 0)
    #self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 1)
    #self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 2)
    #self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 3)
    #self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 4)
    #self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 5)
 
    ##one corr for every fold (not class!)
    #self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 0)
    #self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 1)
    #self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 2)
    #self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 3)
    #self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 4)




    return train_notnan, test_notnan


if __name__ == '__main__':


  ROOT.gROOT.SetBatch(True)
  '''
  #limiting cores/CPU/GPU
  num_cores = 1 #one operation at a time and only one thread per operation  
  num_CPU = 1 
  num_GPU = 0 
  config = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=num_cores, inter_op_parallelism_threads=num_cores, allow_soft_placement = False, device_count = {'CPU' : num_CPU, 'GPU' : num_GPU} )
  session = tf.compat.v1.Session(config=config) 
  tf.compat.v1.keras.backend.set_session(session)         
  '''

  #tf.random.set_seed(1000)
  #np.random.seed(1000)
  
  features = kin_var 
  epochs = 3  #30 here
  #batch_size = 128 #128 here
  batch_size = 2000 #128 here
  scaler_type = 'standard'
  do_early_stopping = True  
  do_reduce_lr = False
  dirname = 'test'
  baseline_selection = baseline_selection 
  nfolds = 10
  frac_sb = 0.0
  frac_sf = 1.0

  trainer = Trainer(
      features = features, 
      epochs = epochs,
      batch_size = batch_size,
      scaler_type = scaler_type,
      do_early_stopping = do_early_stopping,
      do_reduce_lr = do_reduce_lr,
      dirname = dirname,
      baseline_selection = baseline_selection,
      nfolds = nfolds,
      frac_sb = frac_sb,
      frac_sf = frac_sf
      )

  train_samples, test_samples = trainer.process()




