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
import shap
from contextlib import redirect_stdout

#import seaborn as sns

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
def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--constrained', type=boolean_string, required = True)
args = parser.parse_args()

print(f"====> Running constrained fit? {args.constrained}")
shap.initjs()
########################################
#                                      #
# multiclassifier with 4 classes to    #
# optimize signal to background        #
# significance                         #
#                                      #
########################################


#------------------input files-------------------

if args.constrained:
  #this is only constrained data! (cut on fv only)
  file_data  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_cons ]
  # we need to train the NN on the pre-hammered signal!
  # not-hammered
  #file_sig   = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_cons  ]
  # hammered
  direc      = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/signals_default/"
  file_sig      = [os.path.join(direc, f) for f in os.listdir(direc)]

  file_hb    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_cons   ]
  file_b0    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_cons   ]
  file_bs    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_cons   ]
  file_bplus = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_cons]

else:
  #this is only unconstrained data! (cut on fv only)
  file_data  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_unc ]
  file_sig   = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_unc  ]
  file_hb    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_unc   ]
  file_b0    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_unc   ]
  file_bs    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_unc   ]
  file_bplus = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_unc]

#def getRdf(inputfiles, debug = False, skimmed = ""):
#  print(inputfiles)
#  chain = ROOT.TChain("tree")
#
#  for f in inputfiles:
#    chain.AddFile(f)
#
#  rdf = ROOT.RDataFrame(chain)
#  import pdb
#  pdb.set_trace()  
#
#  return rdf

def getRdf(dateTimes, debug = None, skimmed = None):


  chain = ROOT.TChain("tree")

  if not isinstance(dateTimes, list):
    print("dateTimes must be a list of strings")

  if debug:
    print(f"Picking {debug} file(s) for debugging ...")
    fileList = os.listdir(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTimes[0]}/")
    #files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/" +  fileList[0] # pick just the first one 
    for i in range(debug):
      try: chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTimes[0]}/" +  fileList[i])
      except: chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTimes[0]}/*")

    rdf = ROOT.RDataFrame(chain)
    return (chain,rdf)


  for dateTime in dateTimes:

    if skimmed:
      files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/skimmed_{skimmed}_{dateTime}.root" #data skimmed with selection 's     kimmed'
      #files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/skimmed_bkg_{dateTime}.root"  # data skimmed with kkpimu > Bs for      closure
      print(f"Appending {files}")

    else:
      #access the flat ntuples
      files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*" #test

    #chain them all
    chain.Add(files)

  #create rdf from tree

  rdf = ROOT.RDataFrame(chain)

  return (chain,rdf)

#--------------------define sideband region (mlow, mhigh) and dump it into pickle ----------------

if args.constrained:
  chainSig,rdfSig     = getRdf(sig_cons, skimmed = "base_wout_tv") 
  chainData,rdfData   = getRdf(data_cons, skimmed = "base_wout_tv")

else:
  chainSig,rdfSig     = getRdf(sig_unc, skimmed = "base_wout_tv") 
  chainData,rdfData   = getRdf(data_unc, skimmed = "base_wout_tv")


sigma, h          = getSigma(rdfSig, "phiPi_m", base + "& (gen_sig == 0)")

massfit = {}
massfit["sigma"]  = sigma
#massfit["h"]      = h 


if not args.constrained:
  
  A, B, C, S        = getABCS( rdfData, base , "phiPi_m", sigma, h, binsFake = 21,      nSig = nSignalRegion, nSb = nSidebands, width = sbWidth)

  massfit["A"] = A 
  massfit["B"] = B
  massfit["C"] = C
  massfit["S"] = S

#write them into file to use the same constants for the evulation (for later: also use the same for pre-fit plots!)


with open("massfit","wb") as f:
  pickle.dump(massfit,f)

#signal region
mlow   = dsMass_ - nSignalRegion*sigma
mhigh  = dsMass_ + nSignalRegion*sigma

#sideband start
mlow2  = dsMass_ - nSidebands*sigma
mhigh2 = dsMass_ + nSidebands*sigma

#sideband stops
mlow3  = mlow2  - sbWidth*sigma
mhigh3 = mhigh2 + sbWidth*sigma

signalRegion = f"& ({mlow} < phiPi_m) & (phiPi_m < {mhigh})"
leftSB       = f"& ({mlow3} < phiPi_m) & (phiPi_m < {mlow2})"
rightSB      = f"& ({mhigh2} < phiPi_m) & (phiPi_m < {mhigh3})"
lowMass      = f"& (dsMu_m < {bsMass_})"

#------------------limit CPU--------------------------------------

def limitCPU(n):
  print("========> limiting CPU")
  num_threads = n
  os.environ['OMP_NUM_THREADS'] = f'{n}'
  os.environ['TF_NUM_INTRAOP_THREADS'] = f'{n}'
  os.environ['TF_NUM_INTEROP_THREADS'] = f'{n}'

  tf.config.threading.set_inter_op_parallelism_threads(
    num_threads
  )
  tf.config.threading.set_intra_op_parallelism_threads(
    num_threads
  )
  tf.config.set_soft_device_placement(True)

#limitCPU(1)

#--------------------define features ----------------------
 
##the feature vector
kin_var = [
'bs_boost_reco_weighted',
#'bs_boost_reco_1',
#'bs_boost_reco_2',
#'bs_boost_lhcb_alt',
#'bs_boost_coll',

'bs_pt_reco_weighted',
#'bs_pt_reco_1',
#'bs_pt_reco_2',
#'bs_pt_lhcb_alt',
#'bs_pt_coll',

#'cosMuW_reco_weighted', #better separates all signals
#'cosMuW_reco_1', #better separates all signals
#'cosMuW_reco_2', #better separates all signals
#'cosMuW_lhcb_alt', #better separates all signals
#'cosMuW_coll', #better separates all signals

'cosPhiDs_lhcb',
'abs(cosPiK1)',
#'dsMu_deltaR',
'kk_deltaR',

'e_gamma',

#'e_star_reco_weighted',
#'e_star_reco_1',
#'e_star_reco_2',
'e_star_lhcb_alt',
#'e_star_coll',

'm2_miss_coll',
'm2_miss_lhcb_alt',

'mu_rel_iso_03',
'phiPi_deltaR',
'dsMu_m',
#'pt_miss_....',        #too similar to m2 miss?

'q2_reco_weighted',
#'q2_reco_1',
#'q2_reco_2',
#'q2_coll',
'q2_lhcb_alt',
'mu_pt',
'mu_eta',
'mu_phi',
'pi_pt',
#'pi_eta',
#'pi_phi',
#'kk_eta',
#'kk_phi',
#'phi_fitted_pt',
#'ds_fitted_pt',
#'bs_fitted_pt',

'fv_prob',
'tv_prob',
'sv_prob',
#'mu_eta',
#'mu_phi',
#'pi_eta',
#'pi_phi',
#'phiPi_pt',
#'phiPi_phi',
#'phiPi_eta',
#'dsMu_pt',
#'dsMu_eta',
#'dsMu_phi',
'disc_negativity',
'lxy_ds_sig',
'ds_vtx_cosine'
]

if args.constrained: kin_var.append("phiPi_m") 
#if not args.constrained: kin_var.append("phiPi_m") 


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

    if isinstance(filename, list):

      pd_list = []
      for name in filename:
        f = ROOT.TFile(name)
        tree_obj = f.Get(self.tree)

        arr = tree2array(tree_obj,selection = self.selection, branches = kin_var + ['event']) #event will be removed later!
        pd_list.append(pd.DataFrame(arr))

      self.df = pd.concat(pd_list)


    else: 
      f = ROOT.TFile(self.filename)
      tree_obj = f.Get(self.tree)
      arr = tree2array(tree_obj,selection = self.selection, branches = kin_var + ['event']) #event will be removed later!
      self.df = pd.DataFrame(arr)


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


  def createOutDir(self):
    '''
      This function creates the output directory
    '''
    outdir = f'./outputs/{self.dirname}'
    self.outdir = outdir

    if not path.exists(outdir):
      os.system(f'mkdir -p ./outputs/{self.dirname}')

    
    f = open(self.outdir + "/settings.txt", "x")
    f.write("Features: "                  + str(self.features)            + "\n")
    f.write("Epochs: "                    + str(self.epochs)              + "\n")
    f.write("Nfolds: "                    + str(self.nfolds)              + "\n")
    f.write("Batch size: "                + str(self.batch_size)          + "\n")
    f.write("Scaler type: "               + str(self.scaler_type)         + "\n")
    f.write("Early stopping: "            + str(self.do_early_stopping)   + "\n")
    f.write("Reduce lr: "                 + str(self.do_reduce_lr)        + "\n")
    f.write("Baseline selection: "        + str(self.baseline_selection)  + "\n")
    f.write("Constrained fit?"            + str(args.constrained)         + "\n")
    f.write("Signal region up to:"        + str(nSignalRegion) + " sigma" + "\n")
    f.write("Sideband region starts at:"  + str(nSidebands)    + " sigma" + "\n")
    f.write("Sideband region width:"      + str(sbWidth)       + " sigma" + "\n")
    f.write("frac of sideband events:"    + str(self.frac_sb)             + "\n")
    f.write("frac of signflip events:"    + str(self.frac_sf)             + "\n")
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
    hb_selec = " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) && (gen_match_success)"
 
    #signals      #ds mu   #ds tau   #dstar mu   #dstar tau   #hb
    mc_ids      = [0       ,1        ,10         ,11          ,-1]
    mc_classes  = [2       ,0        ,3          ,1           ,4 ] #changed ! 
    label       = {0: r" \mu", 1:r" \tau", 10:r" \mu^{*}", 11: r" \tau^{*}", -1: r" H_{b}" }
    class_label = {class_id: [] for class_id in mc_classes}

    for mc_id, class_id in zip (mc_ids,mc_classes):

      class_label[class_id].append(label[mc_id]) 

      if mc_id >= 0: 
        mc_sample           = Sample(filename = file_sig,    selection=self.baseline_selection  + f'& gen_sig == {mc_id}' + signalRegion + lowMass,                tree = tree_name,signal = class_id).df

      else:
        #mc_sample          = Sample(filename = file_hb,     selection=self.baseline_selection  + hb_selec,               + signalRegion + lowMass,                tree = tree_name,signal = class_id)
        mc_sample1          = Sample(filename = file_b0,     selection=self.baseline_selection  + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df
        mc_sample2          = Sample(filename = file_bs,     selection=self.baseline_selection  + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df
        mc_sample3          = Sample(filename = file_bplus,  selection=self.baseline_selection  + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df

        mc_sample           = pd.concat([mc_sample1,mc_sample2,mc_sample3], sort = False)

      mc_train, mc_test = train_test_split(mc_sample,test_size = 0.2,random_state = 1000)
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

    if args.constrained:
      sign_flip   = " && ((k1_charge*k2_charge > 0) || (pi_charge*mu_charge>0))"
      data_selec  = self.baseline_selection + sign_flip + signalRegion + lowMass

      data_sample = Sample(filename=file_data,   selection=data_selec, tree = 'tree',signal = data_id).df

    else:
      sign_flip    = " && ((k1_charge*k2_charge > 0) || (pi_charge*mu_charge>0))"
      SB           = f'& (( {mlow3} < phiPi_m & phiPi_m < {mlow2}) || ({mhigh2} < phiPi_m & phiPi_m < {mhigh3}))'

      data_selec_sb = self.baseline_selection + SB + lowMass
      data_selec_sf = self.baseline_selection + sign_flip + signalRegion + lowMass

      data_sample_sb = Sample(filename=file_data,   selection=data_selec_sb, tree = 'tree',signal = data_id).df
      data_sample_sf = Sample(filename=file_data,   selection=data_selec_sf, tree = 'tree',signal = data_id).df

      #randomly select frac% of the signflip events
      data_sample_sb = data_sample_sb.sample(frac = self.frac_sb, random_state = 1968)
      data_sample_sf = data_sample_sf.sample(frac = self.frac_sf, random_state = 1968)

      f = open(self.outdir + "/settings.txt", "a")
      f.write("Number of sideband events is: "  + str(len(data_sample_sb)) + "\n")
      f.write("Number of signflip events is: "  + str(len(data_sample_sf)) + "\n")
      f.close()
    
      print("Number of sideband events is: "  + str(len(data_sample_sb)))
      print("Number of signflip events is: "  + str(len(data_sample_sf)))

      data_sample  = pd.concat([data_sample_sb, data_sample_sf], sort = False)

    data_train, data_test = train_test_split(data_sample, test_size=0.2, random_state = 1000)

    train[data_id] = data_train
    test [data_id] = data_test

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
        qt = MinMaxScaler()
      elif self.scaler_type == 'maxabs':
        qt = MaxAbsScaler()

      else:
        raise RuntimeError(f'Unknown scaler {self.scaler_type} - Aborting...')

      qt.fit(X[self.features])
      xx = qt.transform(X[self.features])

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
    train = pd.concat([train_samples[key] for key in train_samples.keys()], sort = False)
    test  = pd.concat([test_samples[key] for key  in test_samples.keys()],  sort = False)

    # even undo the splitting which is not used for k-folding
    main_df = pd.concat([train,test],sort= False)

    # re-index
    main_df.index = np.array(range(len(main_df)))

    # shuffle
    main_df = main_df.sample(frac=1, replace=False, random_state=1986) # of course, keep R's seed ;)

    # X and Y, keep event number for kfold splitting! will be removed before later!
    X = pd.DataFrame(main_df, columns=list(set(self.features)) + ['event'] )
    Y = pd.DataFrame(main_df, columns=['is_signal', 'event'])

    # save the exact list of features
    features_filename = '/'.join([self.outdir, 'input_features.pck'])
    pickle.dump(self.features, open(features_filename, 'wb' ))
    print('========> {} created'.format(features_filename))

    x_folds = {}
    y_folds = {}

    xx_folds = {}

    for n in range(self.nfolds):
      
      #trick to split into folds using event number! :D
      x_folds[n] = X[ X["event"] % self.nfolds == n ]
      y_folds[n] = Y[ Y["event"] % self.nfolds == n ]

      #delete event columns afterwards!
      x_folds[n].drop("event", axis = 1)
      y_folds[n].drop("event", axis = 1)

      # scale the features
      # this is an important step!
  
      xx_folds[n], qt =  self.doScaling(x_folds[n])

      # and save the scaler, which will have to be used throughout the full process, even at evaluation time
      scaler_filename = '/'.join([self.outdir, f'input_tranformation_weighted_fold_{n}.pck'])
      pickle.dump(qt,open(scaler_filename, 'wb'))
      print( ' ========> {} created'.format(scaler_filename))
  
    return main_df, qt, xx_folds, y_folds


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

    l2_rate = 0.01

    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Input((len(features),)))
    #model.add(tf.keras.layers.Dense(40, activation= 'swish'))
    #model.add(tf.keras.layers.Dense(32, activation= 'relu'))
    #model.add(tf.keras.layers.Dropout(.2))
    #model.add(tf.keras.layers.Dense(20, activation ='relu'))
    #model.add(tf.keras.layers.Dense(128, activation ='relu', kernel_regularizer=regularizers.l2(0.01)))
    #model.add(tf.keras.layers.Dropout(.2))
    #model.add(tf.keras.layers.Dense(64, activation ='relu', kernel_regularizer=regularizers.l2(0.01)))
    #model.add(tf.keras.layers.Dense(159, activation ='relu', kernel_regularizer=regularizers.l2(0.01)))
    #model.add(tf.keras.layers.Dense(128, activation ='relu', kernel_regularizer=regularizers.l2(0.01)))
    #model.add(tf.keras.layers.Dense(64, activation ='relu', kernel_regularizer=regularizers.l2(0.01)))
    #model.add(tf.keras.layers.Dense(256,  activation ='swish' )) #, kernel_regularizer=regularizers.l2(0.01)))
    model.add(tf.keras.layers.Dense(64,   activation ='swish', kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(128,  activation ='swish', kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(256,  activation ='swish', kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(128,   activation ='swish', kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(64,   activation ='swish', kernel_regularizer=regularizers.l2(l2_rate)))
    #model.add(tf.keras.layers.Dropout(.2))
    #model.add(tf.keras.layers.Dense(20, activation ='relu')) #, kernel_regularizer=regularizers.l2(0.01)))#, kernel_regularizer=regularizers.l2(0.01))) #also 64 works
    #model.add(tf.keras.layers.Dropout(.2))
    model.add(tf.keras.layers.Dense(6, activation= 'softmax'))


    learning_rate = 0.00005
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
    f = open(self.outdir + "/settings.txt", "a")
    with redirect_stdout(f): 
      model.summary()
      print_model_summary(model)
    f.close()
  
    return model


  def defineCallbacks(self, n):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=8)
    
    # reduce learning rate when at plateau, fine search the minimum
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=5, min_lr=0.000001, cooldown=10, verbose=True)
    
    # save the model every now and then
    filepath = '/'.join([self.outdir, f'fold_{n}' + '_saved-model-{epoch:04d}_val_loss_{val_loss:.4f}_val_acc_{val_acc:.4f}.h5'])
    save_model = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=False, mode='auto', period=1)
    
    callbacks = [save_model]

    if self.do_early_stopping:
      callbacks.append(es)

    if self.do_reduce_lr:
      callbacks.append(reduce_lr)

    #callbacks.append(CustomCallback())

    return callbacks


  def prepareInputs(self, xx_folds, y_folds):
    '''
      Note: the input xx should arlready be scaled
    '''

    # the scaled arrays within xx_folds are not pandas anymore, lets make them pandas again!
    # y_folds still contains pandas, as it has not been scaled!

    for n in range(self.nfolds): 

      xx_folds[n]  = pd.DataFrame(xx_folds[n],  columns=list(set(self.features))) # alternative to X_train[self.features[:]]

      #reset pandas index to start from 0
      xx_folds[n] = xx_folds[n].reset_index(drop=True)
      y_folds[n]  = y_folds[n].reset_index(drop=True)

    return xx_folds, y_folds 


  def train(self, xx_folds, y_folds, weight):
    '''
      Perform the training
    '''

    xx_train = {}
    y_train  = {}
    xx_val   = {}
    y_val    = {}
    model    = {}
    history  = {}

    # create cyclic list, f.e. for nfold = 5, this is:
    # [0,1,2,3,4,0,1,2,3,4]
    cyclic = list(range(self.nfolds)) * 2  

    for n in range(self.nfolds):

      train_indices = [cyclic[n + (i+1)] for i in range(self.nfolds - 2)  ]
      test_index    = cyclic[ n + self.nfolds - 1]

      print(f"For fold/net {n} we train on folds {train_indices} and test on fold {test_index}")

      xx_train[n] = pd.concat([ xx_folds[i] for i in train_indices])
      y_train[n]  = pd.concat([ y_folds[i]  for i in train_indices])
     
      xx_val[n]   = xx_folds[test_index] 
      y_val[n]    = y_folds[test_index] 


      #convert to numpy for training (and use one-hot for y)
      xx_train[n] = xx_train[n].to_numpy()
      xx_val[n]   = xx_val[n].to_numpy()

      y_train[n] = tf.one_hot(y_train[n]['is_signal'].to_numpy(), 6)
      y_train[n] = y_train[n].numpy()

      y_val[n] = tf.one_hot(y_val[n]['is_signal'].to_numpy(), 6)
      y_val[n] = y_val[n].numpy()

      callbacks = self.defineCallbacks(n)

      #define the model every time!! Otherwise you train the same model over and over again lol
      model[n] = self.defineModel()

      history[n] = model[n].fit(xx_train[n], y_train[n], validation_data=(xx_val[n], y_val[n]), epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True,class_weight = weight)
    
    print(f"class weights: {weight}")

    return model, history, xx_train, y_train, xx_val, y_val


  def plotLoss(self, history):
    '''
      Plot the loss for training and validation sets
    '''

    average_loss = []
    

    for n in range(self.nfolds):

      loss_train = history[n].history['loss']
      average_loss.append(np.array(loss_train))
      epochs = range(1, len(loss_train)+1)
      plt.plot(epochs, loss_train, self.colors[n], label=f'Fold {n}')

    #import pdb; pdb.set_trace();
    average_loss = sum(average_loss)/self.nfolds
    plt.plot(epochs, average_loss, 'black' , label='Average loss', linestyle= 'dashed')

    plt.title('Training Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'training_loss')
    plt.clf()



    average_loss = []
    for n in range(self.nfolds):

      loss_val = history[n].history['val_loss']
      average_loss.append(np.array(loss_val))
      epochs = range(1, len(loss_train)+1)
      epochs = range(1, len(loss_train)+1)
      plt.plot(epochs, loss_val, self.colors[n], label=f'Fold {n}')

    average_loss = sum(average_loss)/self.nfolds
    plt.plot(epochs, average_loss, 'black' , label='Average loss', linestyle= 'dashed')
    plt.title('Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'validation_loss')
    plt.clf()


  def plotAccuracy(self, history):
    '''
      Plot the accuracy for training and validation sets
    '''

 
    average_accuracy = []
    for n in range(self.nfolds):

      acc_train = history[n].history['acc']
      average_accuracy.append(np.array(acc_train))
      epochs = range(1, len(acc_train)+1)
      plt.plot(epochs, acc_train, self.colors[n], label='Fold {n}')


    average_accuracy = sum(average_accuracy)/self.nfolds
    plt.plot(epochs, average_accuracy, 'black' , label='Average accuracy', linestyle= 'dashed')

    plt.title('Training Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend()
    self.saveFig(plt, 'training_accuracy')
    plt.clf()

    average_accuracy = []
    for n in range(self.nfolds):

      acc_val = history[n].history['val_acc']
      average_accuracy.append(np.array(acc_val))
      epochs = range(1, len(acc_train)+1)
      plt.plot(epochs, acc_val, self.colors[n], label=f'Fold {n}')

    
    average_accuracy = sum(average_accuracy)/self.nfolds 
    plt.plot(epochs, average_accuracy, 'black' , label='Average accuracy', linestyle= 'dashed')

    plt.title('Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend()
    self.saveFig(plt, 'validation_accuracy')
    plt.clf()




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

   
  def plotScore(self, model, xx_train, y_train, xx_val, y_val, sig, class_label):
    '''
      Plot the score of all channels class sig
    '''

    channels = [0,1,2,3,4,5]
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

        np.savetxt(f"{self.outdir}/fold_{n}_score_of_class_{chan}_for_class_{sig}.csv",y_chan,delimiter = ",")

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
    self.saveFig(fig, f'score_' + str(sig) )
    plt.clf()

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
    
    # rotate axis labels
    g.set_xticklabels(test_x.keys().tolist(), rotation='vertical')
    g.set_yticklabels(test_x.keys().tolist(), rotation='horizontal')
    
    # plt.show()
    #plt.figure(figsize=(12, 12))
    plt.title(f'Linear correlation matrix - Fold {fold}')
    plt.tight_layout()
    self.saveFig(plt, f'corr_sig_fold_{fold}')
    plt.clf()

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
     
    plt.figure()
    col =['c','g','m','r','b','k']
    linestyles = ['solid','dotted','dashed', 'dashdot',(0, (1, 10)),(0,(3, 5, 1, 5, 1, 5))]

    #plot roc for every class 
    for sig in range(5):
      #and plot all folds in one plot
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



  
  def plotKSTest(self, model, xx_train, y_train, xx_val, y_val,sig):
    '''
      Plot the outcome of the Kolmogorov test
      Used to test the overfitting
    '''

    h1 = ROOT.TH1F(f'train_{sig}', f'train_{sig}', 30, 0, 1)
    h2 = ROOT.TH1F(f'val_{sig}', f'val_{sig}', 30, 0, 1)


    for n in range(self.nfolds):

      ## VALIDATION
  
      #class predictions 1D
      #is of shape [1,2,5,4,3,2,4,2,1, ....] 
      true_1d_val = np.argmax(y_val[n], axis=1)
  
      #select only events where the true class is chan...
      x_chan_val    = xx_val[n][ true_1d_val == sig ]
      #...and predict their score! (data is already scaled!)
      y_chan_val    = model[n].predict(x_chan_val)[:,sig] 
  
  
      ## TRAINING
  
      #class predictions 1D
      #is of shape [1,2,5,4,3,2,4,2,1, ....] 
      true_1d_train = np.argmax(y_train[n], axis=1)
  
      #select only events where the true class is chan...
      x_chan_train    = xx_train[n][ true_1d_train == sig ]
      #...and predict their score! (data is already scaled!)
      y_chan_train    = model[n].predict(x_chan_train)[:,sig] 
 
      #import pdb
      #pdb.set_trace()  

 
      np.savetxt(f"{self.outdir}/scoretrain_{sig}.csv",y_chan_train,delimiter = ",")
      np.savetxt(f"{self.outdir}/scoretest_{sig}.csv",y_chan_val,delimiter = ",")
  
      h1_dummy = ROOT.TH1F(f'train_{n}_{sig}', f'train_{n}_{sig}', 30, 0, 1)
      h2_dummy = ROOT.TH1F(f'val_{n}_{sig}', f'val_{n}_{sig}', 30, 0, 1)
      for st,sv in zip(y_chan_train, y_chan_val):
  
        #remark st and sv are lists of length 6! -> only keep score of the signal we're interested in
        h1_dummy.Fill(st) 
        h2_dummy.Fill(sv)
  
      h1.Add(h1_dummy)
      h2.Add(h2_dummy)

    h1.Scale(1/self.nfolds) 
    h2.Scale(1/self.nfolds) 

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
    ks_value.AddText(f'Average KS score of signal {sig} = {round(ks_score,3)}')
    ks_value.SetFillColor(0)
    ks_value.Draw('EP SAME')

    leg = ROOT.TLegend(.55,.82,.83,.88)
    leg.AddEntry(h1 ,'Training' ,'F' )
    leg.AddEntry(h2 ,'Validation' ,'EP' )
    leg.Draw("SAME")


    c1.SaveAs(self.outdir + f'/KS_test_{sig}.pdf')
    c1.SaveAs(self.outdir + f'/KS_test_{sig}.png')
    #print('KS score: ',ks_score, len(train_pred),len(test_pred))



  def process(self):
    print( '------------------- MVA Trainer --------------------')
    
    # create output directory
    print('\n========> creating output directory')
    self.createOutDir()

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
    main_df, qt, xx_folds, y_folds = self.preprocessing(train_notnan, test_notnan)


    print('\n========> preparing the inputs') 
    xx_folds, y_folds = self.prepareInputs(xx_folds, y_folds)

    # do the training
    print('\n========> training...') 
    model, history, xx_train, y_train, xx_val, y_val = self.train(xx_folds, y_folds, weight)

    #save the output dictionaries
    with open(f'{self.outdir}/xx_val.pck', 'wb') as f:
      pickle.dump(xx_val,f)

    with open(f'{self.outdir}/xx_train.pck', 'wb') as f:
      pickle.dump(xx_train,f)

    with open(f'{self.outdir}/y_val.pck', 'wb') as f:
      pickle.dump(y_val,f)

    with open(f'{self.outdir}/y_train.pck', 'wb') as f:
      pickle.dump(y_train,f)

    # plotting
    print('\n========> plotting...' )
    self.plotLoss(history)
    self.plotAccuracy(history)
    self.plotScoreTauOnly(model, xx_train, y_train, xx_val, y_val,class_label )
    self.plotScore(model, xx_train, y_train, xx_val, y_val, 0,class_label)
    self.plotScore(model, xx_train, y_train, xx_val, y_val, 1,class_label)
    self.plotScore(model, xx_train, y_train, xx_val, y_val, 2,class_label)
    self.plotScore(model, xx_train, y_train, xx_val, y_val, 3,class_label)
    self.plotScore(model, xx_train, y_train, xx_val, y_val, 4,class_label)
    self.plotScore(model, xx_train, y_train, xx_val, y_val, 5,class_label)
    #self.plotScoreOneVsAll(model, xx_train, y_train, xx_val, y_val, 1,class_label)
    self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 0)
    self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 1)
    self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 2)
    self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 3)
    self.plotCorr(model,  xx_train, y_train, xx_val, y_val,class_label, 4)
    ####self.plotCM(model, main_test_df, class_label)
    self.plotROCbinary(model,xx_train,y_train,xx_val,y_val,'Train')
    self.plotROCbinary(model,xx_train,y_train,xx_val,y_val,'Test')
    self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 0)
    self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 1)
    self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 2)
    self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 3)
    self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 4)
    self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 5)




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

  tf.random.set_seed(1000)
  np.random.seed(1000)
  
  features = kin_var 
  epochs = 30 #30 here
  batch_size = 128 #128 here
  scaler_type = 'robust'
  do_early_stopping = False
  do_reduce_lr = False
  dirname = 'test'
  baseline_selection = base_wout_tv
  nfolds = 5
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

  trainer.process()




