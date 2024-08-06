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

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--constrained', type=boolean_string, required = True)
args = parser.parse_args()

print(f"====> Running constrained fit? {args.constrained}")

########################################
#                                      #
# multiclassifier with 4 classes to    #
# optimize signal to background        #
# significance                         #
#                                      #
########################################


#------------------input files-------------------

#this is only unconstrained data! (cut on fv and tv)
files_data  = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in data_unc ] 
files_sig   = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in sig_unc  ]
files_hb    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in hb_unc   ]
files_b0    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in b0_unc   ]
files_bs    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in bs_unc   ]
files_bplus = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in bplus_unc]
#this is only unconstrained data! (cut on fv only)
file_data  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_unc ]
file_sig   = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_unc  ]
file_hb    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_unc   ]
file_b0    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_unc   ]
file_bs    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_unc   ]
file_bplus = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_unc]
#this is only unconstrained data! (cut on fv only)

#file_data  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_cons ]
#file_sig   = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_cons  ]
#file_hb    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_cons   ]
#file_b0    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_cons   ]
#file_bs    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_cons   ]
#file_bplus = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_cons]

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
  chainSig,rdfSig     = getRdf(sig_cons, skimmed = "base") 
  chainData,rdfData   = getRdf(data_cons, skimmed = "base")

else:
  chainSig,rdfSig     = getRdf(sig_unc, skimmed = "base") 
  chainData,rdfData   = getRdf(data_unc, skimmed = "base")


sigma, h          = getSigma(rdfSig, "phiPi_m", ma_cut_wout_tv + "& (gen_sig == 0)")

massfit = {}
massfit["sigma"]  = sigma
#massfit["h"]      = h 


if not args.constrained:
  
  A, B, C, S        = getABCS( rdfData, ma_cut_wout_tv , "phiPi_m", sigma, h, binsFake = 21,      nSig = nSignalRegion, nSb = nSidebands, width = sbWidth)

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
'bs_boost_lhcb_alt',
'bs_boost_coll',

'bs_pt_reco_weighted',
#'bs_pt_reco_1',
#'bs_pt_reco_2',
'bs_pt_lhcb_alt',
'bs_pt_coll',

'cosMuW_reco_weighted', #better separates all signals
#'cosMuW_reco_1', #better separates all signals
#'cosMuW_reco_2', #better separates all signals
'cosMuW_lhcb_alt', #better separates all signals
'cosMuW_coll', #better separates all signals

'cosPhiDs_lhcb',
'cosPiK1',
'dsMu_deltaR',
'kk_deltaR',

'e_gamma',

'e_star_reco_weighted',
#'e_star_reco_1',
#'e_star_reco_2',
'e_star_lhcb_alt',
'e_star_coll',

'm2_miss_coll',
'm2_miss_lhcb_alt',

'mu_rel_iso_03',
'phiPi_deltaR',
#'phiPi_m',              #only for constrained fitter!
'dsMu_m',
#'pt_miss_....',        #too similar to m2 miss?

'q2_reco_weighted',
#'q2_reco_1',
#'q2_reco_2',
'q2_coll',
'q2_lhcb_alt',
'mu_pt',
'pi_pt',

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
'ds_vtx_cosine'
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

    if isinstance(filename, list):

      pd_list = []
      for name in filename:
        f = ROOT.TFile(name)
        tree_obj = f.Get(self.tree)

        arr = tree2array(tree_obj,selection = self.selection, branches = kin_var)
        pd_list.append(pd.DataFrame(arr))

      self.df = pd.concat(pd_list)


    else: 
      f = ROOT.TFile(self.filename)
      tree_obj = f.Get(self.tree)
      arr = tree2array(tree_obj,selection = self.selection, branches = kin_var)
      self.df = pd.DataFrame(arr)


    pd.options.mode.chained_assignment = None
    self.df['is_signal'] = signal


class Trainer(object):

  'train the data'

  def __init__(self, features, epochs, batch_size, scaler_type, do_early_stopping, do_reduce_lr, dirname, baseline_selection):
    self.features = features #variables to train on
    self.epochs = epochs #samples / batch_size =  number of iterations to 1 epoch
    self.batch_size = batch_size #number of samples which we feed to our model
    self.scaler_type = scaler_type
    self.do_early_stopping = do_early_stopping 
    self.do_reduce_lr = do_reduce_lr
    self.dirname = dirname + '_' + datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss')
    self.baseline_selection = baseline_selection


  def createOutDir(self):
    '''
      This function creates the output directory
    '''
    outdir = f'./outputs/{self.dirname}'
    if not path.exists(outdir):
      os.system(f'mkdir -p ./outputs/{self.dirname}')

    return outdir


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
 
    #signals     #ds mu   #ds tau   #dstar mu   #dstar tau   #hb
    mc_ids     = [0       ,1        ,10         ,11          ,-1]
    mc_classes = [0       ,0        ,0          ,0           ,1 ] #changed ! 

    for mc_id, class_id in zip (mc_ids,mc_classes):

      if mc_id >= 0: 
        mc_sample          = Sample(filename = file_sig, selection=self.baseline_selection  + f'& gen_sig == {mc_id}', tree = tree_name,signal = class_id).df

      else:
        #mc_sample          = Sample(filename = file_hb,     selection=self.baseline_selection  + hb_selec,                tree = tree_name,signal = class_id)
        mc_sample1          = Sample(filename = file_b0,     selection=self.baseline_selection  + hb_selec,                tree = tree_name,signal = class_id).df
        mc_sample2          = Sample(filename = file_bs,     selection=self.baseline_selection  + hb_selec,                tree = tree_name,signal = class_id).df
        mc_sample3          = Sample(filename = file_bplus,  selection=self.baseline_selection  + hb_selec,                tree = tree_name,signal = class_id).df
        mc_sample           = pd.concat([mc_sample1,mc_sample2,mc_sample3], sort = False)

      mc_train, mc_test = train_test_split(mc_sample,test_size = 0.2,random_state = 1000)
      print(f"train sample for signal {mc_id} has class name {class_id} and contains {len(mc_train)} events")
      print(f"test sample for signal {mc_id} has class name {class_id} and contains {len(mc_test)} events")
  
      train[class_id] = mc_train
      test [class_id] = mc_test

    #---------------------------  COMB BACKGROUND ------------------------------------ 

    ## Denote comb with id 2
    data_id = 2

    #siebands used for the comb. bkg

    if args.constrained:
      sign_flip  = "((k1_charge*k2_charge > 0) || (pi_charge*mu_charge>0))"
      data_selec = sign_flip

    else:
      SB = f'& (( {mlow3} < phiPi_m & phiPi_m < {mlow2}) || ({mhigh2} < phiPi_m & phiPi_m < {mhigh3}))'
      data_selec = self.baseline_selection + SB

    data_sample           = Sample(filename=file_data, selection=data_selec, tree = 'tree',signal = data_id)
    data_train, data_test = train_test_split(data_sample.df, test_size=0.2, random_state = 1000)

    train[data_id] = data_train
    test [data_id] = data_test

    print(f"train sample for data has class name {data_id} and  has {len(data_train)} events")
    print(f"test sample for data has class name {data_id} and  has {len(data_test)} events")


    
    print(f'========> it took {round(time() - now,2)} seconds' )

    return train, test 

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

      qt.fit(X[features])
      xx = qt.transform(X[features])

      return xx, qt


  def preprocessing(self, samples):
    '''
      Preprocessing of data before training/testing the NN
      This includes:
        - building the main_df
        - building the scaler
        - get the scaled features xx
        - get the target Y
    '''

    # concatenate all the datasets and shuffle events
    main_df = pd.concat([samples[key] for key in samples.keys()],sort= False)

    # re-index
    main_df.index = np.array(range(len(main_df)))

    # shuffle
    main_df = main_df.sample(frac=1, replace=False, random_state=1986) # of course, keep R's seed ;)

    # X and Y
    X = pd.DataFrame(main_df, columns=list(set(self.features)))
    Y = pd.DataFrame(main_df, columns=['is_signal'])

    # scale the features
    # this is an important step!
    xx, qt = self.doScaling(X)

    # and save the scaler, which will have to be used throughout the full process, even at evaluation time
    scaler_filename = '/'.join([self.outdir, 'input_tranformation_weighted.pck'])
    pickle.dump(qt,open(scaler_filename, 'wb'))
    print( ' ========> {} created'.format(scaler_filename))

    # save the exact list of features
    features_filename = '/'.join([self.outdir, 'input_features.pck'])
    pickle.dump(self.features, open(features_filename, 'wb' ))
    print('========> {} created'.format(features_filename))

    return main_df, qt, xx, Y


  def defineModel(self):
    '''
      Define the NN
    '''
    #NOTE for the moment, everything is hardcoded

    rate = 0.0001 #learning rate (not applied in our model :))
    rate_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
    initial_learning_rate=0.0001,
    decay_steps=10000,
    decay_rate=0.9
)

    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Input((len(features),)))
    #model.add(tf.keras.layers.Dense(40, activation= 'swish'))
    #model.add(tf.keras.layers.Dense(32, activation= 'relu'))
    #model.add(tf.keras.layers.Dropout(.2))
    #model.add(tf.keras.layers.Dense(20, activation ='relu'))
    model.add(tf.keras.layers.Dense(128, activation ='swish'))#, kernel_regularizer=regularizers.l2(0.01)))
    #model.add(tf.keras.layers.Dropout(.2))
    model.add(tf.keras.layers.Dense(64, activation ='swish', kernel_regularizer=regularizers.l2(0.01)))
    model.add(tf.keras.layers.Dense(64, activation ='swish', kernel_regularizer=regularizers.l2(0.01)))
    #model.add(tf.keras.layers.Dropout(.2))
    model.add(tf.keras.layers.Dense(32, activation ='swish'))#, kernel_regularizer=regularizers.l2(0.01))) #also 64 works
    #model.add(tf.keras.layers.Dropout(.2))
    model.add(tf.keras.layers.Dense(3, activation= 'softmax'))

    opt = keras.optimizers.SGD(learning_rate=rate_schedule) 
    model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=['acc'])
    
    print(model.summary())

    return model


  def defineCallbacks(self):
    '''
      Define the callbacks
    '''
    # early stopping
    monitor = 'val_loss'
    es = EarlyStopping(monitor=monitor, mode='auto', verbose=1, patience=8)
    
    # reduce learning rate when at plateau, fine search the minimum
    reduce_lr = ReduceLROnPlateau(monitor=monitor, mode='auto', factor=0.2, patience=5, min_lr=0.000001, cooldown=10, verbose=True)
    
    # save the model every now and then
    filepath = '/'.join([self.outdir, 'saved-model-{epoch:04d}_val_loss_{val_loss:.4f}_val_acc_{val_acc:.4f}.h5'])
    save_model = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, save_weights_only=False, mode='auto', period=1)
    
    callbacks = [save_model]

    if self.do_early_stopping:
      callbacks.append(es)

    if self.do_reduce_lr:
      callbacks.append(reduce_lr)

    #callbacks.append(CustomCallback())

    return callbacks


  def prepareInputs(self, xx, Y,xx_test,Y_test):
    '''
      Note: the input xx should arlready be scaled
    '''

    x_train, x_val, y_train, y_val = xx, xx_test,Y,Y_test    
    # the x should only contain the features and not all the branches of main_df
    x_train = pd.DataFrame(x_train, columns=list(set(self.features))) # alternative to X_train[self.features[:]]
    x_val   = pd.DataFrame(x_val,   columns=list(set(self.features)))

    x_train = x_train.reset_index(drop=True)
    x_val   = x_val.reset_index(drop=True)
    y_train = y_train.reset_index(drop=True)
    y_val   = y_val.reset_index(drop=True)

    return x_train, x_val, y_train, y_val


  def train(self, model, x_train, y_train, x_val, y_val, callbacks,weight):
    '''
      Perform the training
    '''

    x_train = x_train.to_numpy()
    y_train = tf.one_hot(y_train['is_signal'].to_numpy(), 3)
    y_train = y_train.numpy()

    x_val = x_val.to_numpy()
    y_val= tf.one_hot(y_val['is_signal'].to_numpy(), 3)
    y_val = y_val.numpy()
    
    print(f"class weights: {weight}")

    history = model.fit(x_train, y_train, validation_data=(x_val, y_val), epochs=self.epochs, callbacks=callbacks, batch_size=self.batch_size, verbose=True,class_weight = weight)

    return history


  def plotLoss(self, history):
    '''
      Plot the loss for training and validation sets
    '''
    loss_train = history.history['loss']
    loss_val = history.history['val_loss']
    epochs = range(1, len(loss_train)+1)
    plt.plot(epochs, loss_train, 'g', label='Training loss')
    plt.plot(epochs, loss_val, 'b', label='Validation loss')
    plt.title('Training and Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    self.saveFig(plt, 'loss')
    plt.clf()


  def plotAccuracy(self, history):
    '''
      Plot the accuracy for training and validation sets
    '''
    acc_train = history.history['acc']
    acc_val = history.history['val_acc']
    epochs = range(1, len(acc_train)+1)
   
    plt.plot(epochs, acc_train, 'g', label='Training accuracy')
    plt.plot(epochs, acc_val, 'b', label='Validation accuracy')
    plt.title('Training and Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend()
    self.saveFig(plt, 'accuracy')
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


  def plotScoreTau(self, model, mc_test_samples,data_test_samples):
    '''
      Plot the score distributions for signal and background
    '''
    signals = [0,1,2,3,4,5]

    # get the score for the test dataframe

    mc_test_samples = pd.concat(mc_test_samples,sort= False)
    data_test_samples = pd.concat(data_test_samples, sort = False)

    #mu 
    df0 =  mc_test_samples.loc[mc_test_samples['is_signal'] == 0]
    df1 =  mc_test_samples.loc[mc_test_samples['is_signal'] == 2]
    
    #back
    df4 =  mc_test_samples.loc[mc_test_samples['is_signal'] == 4]
    df5 =  data_test_samples.loc[data_test_samples['is_signal'] == 5]
  
    #tau
    df2 = mc_test_samples.loc[mc_test_samples['is_signal'] == 1]
    df3 = mc_test_samples.loc[mc_test_samples['is_signal'] == 3]

    #for part two :)
    df = pd.concat([df2,df3],sort = False)
    df_rest = pd.concat([df0,df1,df4,df5],sort = False)

    #collect initial length 
    lmu = len(df0) + len(df1)
    ltau2 = len(df2)
    ltau3 = len(df3)
    ltau = ltau2 + ltau3
    lback = len(df4) + len(df5) 

    print(f"we start off with {len(df0)+len(df1)} mu samples, {len(df4)} hb samples, {len(df5)} comb bkg samples and {len(df2) + len(df3)} tau samples")

 
    #calculat mu score
    score_mu0 = self.predictScore(model,df0) 
    score_mu1 = self.predictScore(model,df1) 
    #only keep score columns of tau and tau* and comb! We're interested in these!
    score_mu2 = np.concatenate((score_mu0[:,2],score_mu1[:,2])) #mu score at tau
    score_mu3 =  np.concatenate((score_mu0[:,3],score_mu1[:,3])) #mu_score at tau*
    score_mu5 = np.concatenate((score_mu0[:,5],score_mu1[:,5]))

    #same for tau itself
    score_tau2 =  self.predictScore(model,df2)
    score_tau3 =  self.predictScore(model,df3)
    score_tau2_v2 = np.concatenate((score_tau2[:,2],score_tau3[:,2])) #tau score at tau
    score_tau3_v2 = np.concatenate((score_tau2[:,3],score_tau3[:,3])) #tau score at tau*
    score_tau5 = np.concatenate((score_tau2[:,5],score_tau3[:,5]))
    print(len(score_tau2))
 
    #and same for background
    score_back4 = self.predictScore(model,df4)   
    score_back5 = self.predictScore(model,df5)   
    score_back2 = np.concatenate((score_back4[:,2],score_back5[:,2])) #back score at tau
    score_back3 = np.concatenate((score_back4[:,3],score_back5[:,3])) #back score at tau*
    score_back5_v2 = np.concatenate((score_back4[:,5],score_back4[:,5]))

    #cut array over which we optimize
    cuts = np.arange(0,0.8,0.05) 
     
    #prepare 2D arrays, 2D because we want to cut tau and tau* independently 
    tau2_left = np.zeros((len(cuts),len(cuts),len(cuts),len(cuts),len(cuts)))
    tau3_left = np.zeros_like(tau2_left)
    back5_left = np.zeros_like(tau2_left)
    back4_left = np.zeros_like(tau2_left)
    mu0_left = np.zeros_like(tau2_left)
    mu1_left = np.zeros_like(tau2_left)
    mu_left = np.zeros_like(tau2_left)
    tau_left = np.zeros_like(tau2_left)
    back_left = np.zeros_like(tau2_left)
    sig_to_back = np.zeros_like(tau2_left)
    sig_to_real_back = np.zeros_like(tau2_left)
    all_sig_to_back =  np.zeros_like(tau2_left)

    for i,c2 in enumerate(cuts):
      print(c2)
      for j,c3 in enumerate(cuts):
        for k,c5 in enumerate(cuts):
          for l,c0 in enumerate(cuts):
            for m,c1 in enumerate(cuts):

              #ugly but it works! :) 
              tau2_count = len(np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(score_tau2[:,2] > c2).flatten(),np.argwhere(score_tau2[:,3] > c3).flatten()),np.intersect1d(np.argwhere(score_tau2[:,5] < c5).flatten() ,np.argwhere(score_tau2[:,0] > c0).flatten())) ,np.argwhere(score_tau2[:,1] > c1).flatten()))
              tau3_count = len(np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(score_tau3[:,2] > c2).flatten(),np.argwhere(score_tau3[:,3] > c3).flatten()),np.intersect1d(np.argwhere(score_tau3[:,5] < c5).flatten() ,np.argwhere(score_tau3[:,0] > c0).flatten())) ,np.argwhere(score_tau3[:,1] > c1).flatten()))
              mu0_count = len(np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(score_mu0[:,2] > c2).flatten(),np.argwhere(score_mu0[:,3] > c3).flatten()),np.intersect1d(np.argwhere(score_mu0[:,5] < c5).flatten() ,np.argwhere(score_mu0[:,0] > c0).flatten())) ,np.argwhere(score_mu0[:,1] > c1).flatten()))
              mu1_count = len(np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(score_mu1[:,2] > c2).flatten(),np.argwhere(score_mu1[:,3] > c3).flatten()),np.intersect1d(np.argwhere(score_mu1[:,5] < c5).flatten() ,np.argwhere(score_mu1[:,0] > c0).flatten())) ,np.argwhere(score_mu1[:,1] > c1).flatten()))
              back4_count = len(np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(score_back4[:,2] > c2).flatten(),np.argwhere(score_back4[:,3] > c3).flatten()),np.intersect1d(np.argwhere(score_back4[:,5] < c5).flatten() ,np.argwhere(score_back4[:,0] > c0).flatten())) ,np.argwhere(score_back4[:,1] > c1).flatten()))
              back5_count = len(np.intersect1d(np.intersect1d(np.intersect1d(np.argwhere(score_back5[:,2] > c2).flatten(),np.argwhere(score_back5[:,3] > c3).flatten()),np.intersect1d(np.argwhere(score_back5[:,5] < c5).flatten() ,np.argwhere(score_back5[:,0] > c0).flatten())) ,np.argwhere(score_back5[:,1] > c1).flatten()))

              #see how many tau* are left if we cut tau on c2 and tau* on c3: 
              tau2_left[i][j][k][l][m] = tau2_count
              tau3_left[i][j][k][l][m] = tau3_count
              back4_left[i][j][k][l][m] = back4_count
              back5_left[i][j][k][l][m] = back5_count
              mu0_left[i][j][k][l][m] = mu0_count
              mu1_left[i][j][k][l][m] = mu1_count
              tau_left[i][j][k][l][m] = tau2_count + tau3_count
              mu_left[i][j][k][l][m] = mu0_count + mu1_count
              back_left[i][j][k][l][m] = back4_count + back5_count

              tot_count = tau2_count + tau3_count + back4_count + back5_count + mu0_count + mu1_count
          
              if tot_count == 0:
                sig_to_back[i][j][k][l][m] = np.nan
                all_sig_to_back[i][j][k][l][m] = np.nan

              else:
                sig_to_back[i][j][k][l][m] = (tau2_count + tau3_count) / np.sqrt(tot_count)
                all_sig_to_back[i][j][k][l][m] = (tau2_count + tau3_count + mu0_count + mu1_count) / np.sqrt(tot_count)
 
              if back4_count + back5_count + tau2_count + tau3_count == 0:
                sig_to_real_back[i][j][k][l][m] = np.nan
              else:
                sig_to_real_back[i][j][k][l][m] = (tau2_count + tau3_count) / np.sqrt( tau2_count + tau3_count + back_left[i][j][k][l][m])

    

    c2max = np.unravel_index(np.nanargmax(sig_to_back), sig_to_back.shape)[0]
    c3max = np.unravel_index(np.nanargmax(sig_to_back), sig_to_back.shape)[1]
    c5max = np.unravel_index(np.nanargmax(sig_to_back), sig_to_back.shape)[2]
    c0max = np.unravel_index(np.nanargmax(sig_to_back), sig_to_back.shape)[3]
    c1max = np.unravel_index(np.nanargmax(sig_to_back), sig_to_back.shape)[4]

    c2max_real = np.unravel_index(np.nanargmax(sig_to_real_back), sig_to_real_back.shape)[0]
    c3max_real = np.unravel_index(np.nanargmax(sig_to_real_back), sig_to_real_back.shape)[1]
    c5max_real = np.unravel_index(np.nanargmax(sig_to_real_back), sig_to_real_back.shape)[2]
    c0max_real = np.unravel_index(np.nanargmax(sig_to_real_back), sig_to_real_back.shape)[3]
    c1max_real = np.unravel_index(np.nanargmax(sig_to_real_back), sig_to_real_back.shape)[4]

    c2max_all = np.unravel_index(np.nanargmax(all_sig_to_back), all_sig_to_back.shape)[0]
    c3max_all = np.unravel_index(np.nanargmax(all_sig_to_back), all_sig_to_back.shape)[1]
    c5max_all = np.unravel_index(np.nanargmax(all_sig_to_back), all_sig_to_back.shape)[2]
    c0max_all = np.unravel_index(np.nanargmax(all_sig_to_back), all_sig_to_back.shape)[3]
    c1max_all = np.unravel_index(np.nanargmax(all_sig_to_back), all_sig_to_back.shape)[4]



    print(f"Sig to Back is {sig_to_back[c2max][c3max][c5max][c0max][c1max]} maximal at c2 = {cuts[c2max]}, c3 = {cuts[c3max]}, c5 = {cuts[c5max]}, c0 = {cuts[c0max]}, c1 = {cuts[c1max]}")

    print(f"and there we have {back_left[c2max][c3max][c5max][c0max][c1max]/lback} bkg. left")
    print(f"and there we have {mu_left[c2max][c3max][c5max][c0max][c1max]/lmu} mu left")
    print(f"and there we have {tau_left[c2max][c3max][c5max][c0max][c1max]/ltau} tau left")


    print(f"Sig to real Back is {sig_to_real_back[c2max_real][c3max_real][c5max_real][c0max_real][c1max_real]} maximal at c2 = {cuts[c2max_real]}, c3 = {cuts[c3max_real]}, c5 = {cuts[c5max_real]}, c0 = {cuts[c0max_real]}, c1 = {cuts[c1max_real]}")

    print(f"and there we have {back_left[c2max_real][c3max_real][c5max_real][c0max_real][c1max_real]/lback} bkg. left")
    print(f"and there we have {mu_left[c2max_real][c3max_real][c5max_real][c0max_real][c1max_real]/lmu} mu left")
    print(f"and there we have {tau_left[c2max_real][c3max_real][c5max_real][c0max_real][c1max_real]/ltau} tau left")


    print(f"All signals to Back, {sig_to_real_back[c2max_real][c3max_real][c5max_real][c0max_real][c1max_real]} maximal at c2 = {cuts[c2max_all]}, c3 = {cuts[c3max_all]}, c5 = {cuts[c5max_all]}, c0 = {cuts[c0max_all]}, c1 = {cuts[c1max_all]}")

    print(f"and there we have {back_left[c2max_all][c3max_all][c5max_all][c0max_all][c1max_all]/lback} bkg. left")
    print(f"and there we have {mu_left[c2max_all][c3max_all][c5max_all][c0max_all][c1max_all]/lmu} mu left")
    print(f"and there we have {tau_left[c2max_all][c3max_all][c5max_all][c0max_all][c1max_all]/ltau} tau left")



    ############

    print("=========> collapsing multiclassifier to binary (tau against the rest)")


    max_ind = []    
    max_ind_rest = []


    #score of tau signals, get index of maximal score
    for event in np.concatenate((score_tau2,score_tau3)):
      max_ind.append(np.nanargmax(event))
    #score of all other signals, get index of maximal score
    for event in np.concatenate((score_mu0,score_mu1,score_back4,score_back5)):
      max_ind_rest.append(np.nanargmax(event))

    max_ind = np.array(max_ind)
    max_ind_rest = np.array(max_ind_rest)

    true_pos = (len(max_ind[max_ind == 2]) + len(max_ind[max_ind == 3]))/len(max_ind)
    false_neg =(len(max_ind) - len(max_ind[max_ind == 2]) - len(max_ind[max_ind == 3]))/len(max_ind)

    false_pos = (len(max_ind_rest[max_ind_rest == 2]) + len(max_ind_rest[max_ind_rest == 3])) /len(max_ind_rest)
    true_neg = (len(max_ind_rest) - len(max_ind_rest[max_ind_rest == 2]) - len(max_ind_rest[max_ind_rest == 3]))/len(max_ind_rest)


    print(f"True positive rate: {true_pos}, False positive rate: {false_pos}")
    print(f"True negative rate: {true_neg}, False negative rate: {false_neg}")

  def plotScore(self, model, mc_test_samples,data_test_samples,sig):
    '''
      Plot the score of all channels class sig
    '''

    mc_test_samples = pd.concat(mc_test_samples,sort= False)
    data_test_samples = pd.concat(data_test_samples, sort = False)

    #mu 
    df0 =  mc_test_samples.loc[mc_test_samples['is_signal'] == 0]
    df1 =  mc_test_samples.loc[mc_test_samples['is_signal'] == 10]
    
    #back
    df4 =  mc_test_samples.loc[mc_test_samples['is_signal'] == 1]
    df5 =  data_test_samples.loc[data_test_samples['is_signal'] == 11]
  
    #tau
    df2 = mc_test_samples.loc[mc_test_samples['is_signal'] == -1]
    df3 = mc_test_samples.loc[mc_test_samples['is_signal'] == -2]

    scores = []

    for i,df in enumerate([df0,df1,df2,df3,df4,df5]):
      #only keep score of signal sig
      score = self.predictScore(model,df)[:,sig]
      #save it
      np.savetxt(f"{self.outdir}/score_{i}_for_{sig}.csv",score,delimiter = ",")
      scores.append(self.predictScore(model,df)[:,sig]) #only keep sig class score 
  
    name={0:r'$D_s \mu$',1:r'$D_s \tau$',10:r'$D^*_s \mu$',11:r'$D^*_s \tau$',-1:r'$H_b$',-2:'comb. Bkg.'}
    col ={0:'c',1:'g',10:'m',11:'r',5:'b',6:'k'}
    linestyles = {0:'solid',1:'dotted',10:'dashed',11:'dashdot',-1:(0, (1, 10)),-2:(0, (3, 5, 1, 5, 1, 5))}

    fig = plt.figure()

    for i,score in enumerate(scores):
      # plot the score distributions
      plt.hist(score, bins=np.arange(0,1.025,0.025), color=col[i], alpha=0.5, label=f'actual class {i}', histtype='stepfilled',density = True,linestyle = linestyles[i], linewidth = 1.5)
    plt.legend(loc='upper right')
    plt.title(f'Predicition for class {sig}')
    plt.xlabel('score')
    plt.ylabel('events')
    fig.savefig('outputs/score_' + str(sig) + '.pdf')
    fig.savefig('outputs/score_' + str(sig) + '.png')
    self.saveFig(fig, 'score_' + str(sig) )
    plt.clf()



  def plotCM(self, model, mc_test_samples,data_test_samples):
    '''
      Get and save everything rto later produce confusion matrix
    '''

    mc = pd.concat(mc_test_samples,sort = False)
    data = pd.concat(data_test_samples,sort = False)

    test = pd.concat([mc,data],sort = False)
    test_x = pd.DataFrame(test, columns=list(set(self.features)))
    test_y = pd.DataFrame(test, columns=['is_signal'])
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
    print(cm)


  def plotROCbinary(self,model,x_train,y_train,x_val,y_val,key):

    'plot one vs all roc curve'

    #for one-vs-all roc curve
    label_binarizer = LabelBinarizer().fit(y_train)

    #plot train or test roc curve
    if key == 'Train':
      x = x_train
      y = y_train
    else:
      x = x_val
      y = y_val

    #for one-vs-all roc curve
    y_score = model.predict(x)
    y_onehot_test = label_binarizer.transform(y)
   
    plt.figure()
    col =['c','g','m','r'] #,'b','k']
    linestyles = ['solid','dotted','dashed', 'dashdot' ]#,(0, (1, 10)),(0,(3, 5, 1, 5, 1, 5))]


    for sig in range(3):
      #plot roc for every class
      class_id = np.flatnonzero(label_binarizer.classes_ == sig)[0]
      fpr,tpr,_ = roc_curve(y_onehot_test[:, class_id],y_score[:, class_id])    
      #save it
      np.savetxt(f"{self.outdir}/tprfpr_{sig}_{key}.csv",np.array([fpr,tpr]),delimiter = ",")
      roc_auc = auc(tpr,fpr)
      plt.plot(fpr, tpr, label=f'{key} ROC of class {sig}, AUC = {round(roc_auc,2)}',color = col[sig],linestyle = linestyles[sig])

    #baseline
    #true_pos_rate_baseline = 0.5827968065928406 #from ./NNmodels/log_plotter.txt
    #false_pos_rate_baseline = 0.1972626960662999 # "

    #plt.scatter(false_pos_rate_baseline,true_pos_rate_baseline,marker = 'o', label = 'Baseline',color = 'orange')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    xy = [i*j for i,j in product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
    plt.plot(xy, xy, color='grey', linestyle='--',label = 'No skill')
    plt.yscale('linear')
    plt.legend()
    self.saveFig(plt, f'ROC_{key}')
    plt.clf()



  
  def plotKSTest(self, model, x_train, x_val, y_train, y_val,sig):
    '''
      Plot the outcome of the Kolmogorov test
      Used to test the overfitting
    '''


    indices_train =list(np.where(y_train["is_signal"] == sig)[0])
    indices_val = list(np.where(y_val["is_signal"] == sig)[0])

    #only keep signal sig
    x_train_part = x_train.iloc[indices_train] 
    x_val_part =x_val.iloc[indices_val] 

 
    #input has to be numpy array when directly given to .predict:
    x_train_part = x_train_part.to_numpy()
    x_val_part = x_val_part.to_numpy()

    #training data
    score_train = model.predict(x_train_part)
    #validation data
    score_val = model.predict(x_val_part)

    np.savetxt(f"{self.outdir}/scoretrain_{sig}.csv",score_train,delimiter = ",")
    np.savetxt(f"{self.outdir}/scoretest_{sig}.csv",score_val,delimiter = ",")

    h1 = ROOT.TH1F(f'train_{sig}', f'train_{sig}', 30, 0, 1)
    h2 = ROOT.TH1F(f'val_{sig}', f'val_{sig}', 30, 0, 1)
    for st,sv in zip(score_train, score_val):

      #remark st and sv are lists of length 6! -> only keep score of the signal we're interested in
      h1.Fill(st[sig]) 
      h2.Fill(sv[sig])

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
    ks_value.AddText(f'KS score of signal {sig} = {round(ks_score,3)}')
    ks_value.SetFillColor(0)
    ks_value.Draw('EP SAME')

    leg = ROOT.TLegend(.55,.82,.83,.88)
    leg.AddEntry(h1 ,'Training' ,'F' )
    leg.AddEntry(h2 ,'Validation' ,'EP' )
    leg.Draw("SAME")


    c1.SaveAs(self.outdir + f'/KS_test_{sig}.pdf')
    #print('KS score: ',ks_score, len(train_pred),len(test_pred))



  def process(self):
    print( '------------------- MVA Trainer --------------------')
    
    # create output directory
    print('\n========> creating output directory')
    self.outdir = self.createOutDir()

    # get the samples (dictionaries with signal key)
    print('\n========> getting the samples')
    train_samples, test_samples = self.getSamples()

    # create dataframes
    print('\n========> creating the train dataframes')
    train_notnan, weight      = self.createDataframe(train_samples)
    print('\n========> creating the test dataframes')
    test_notnan , _           = self.createDataframe(test_samples)

    # preprocessing the dataframes
    print('\n========> preprocessing the train dataframes' )
    main_df, qt, xx, Y = self.preprocessing(train_notnan)
    print('\n========> preprocessing the test dataframes' )
    main_test_df, qt_test, xx_test, Y_test = self.preprocessing(test_notnan)

    # define the NN
    print('\n========> defining the model' )
    model = self.defineModel()

    # define the callbacks
    print('\n========> defining the callbacks' )
    callbacks = self.defineCallbacks()

    print('\n========> preparing the inputs') 
    x_train, x_val, y_train, y_val = self.prepareInputs(xx,Y,xx_test,Y_test)
    # do the training
    print('\n========> training...') 
    history = self.train(model, x_train, y_train, x_val, y_val, callbacks,weight)

    # plotting
    print('\n========> plotting...' )
    self.plotLoss(history)
    self.plotAccuracy(history)
    #self.plotScoreTau(model,mc_test_df,data_test_df)
    #self.plotScore(model, mc_test_df,data_test_df,-2)
    #self.plotScore(model, mc_test_df,data_test_df,0)
    #self.plotScore(model, mc_test_df,data_test_df,1)
    #self.plotScore(model, mc_test_df,data_test_df,10)
    #self.plotScore(model, mc_test_df,data_test_df,11)
    #self.plotScore(model, mc_test_df,data_test_df,-1)
    #self.plotCM(model, mc_test_df,data_test_df)
    self.plotROCbinary(model,x_train,y_train,x_val,y_val,'Train')
    self.plotROCbinary(model,x_train,y_train,x_val,y_val,'Test')
    #self.plotKSTest(model, x_train, x_val, y_train, y_val, -2)
    self.plotKSTest(model, x_train, x_val, y_train, y_val, 0)
    self.plotKSTest(model, x_train, x_val, y_train, y_val, 1)
    self.plotKSTest(model, x_train, x_val, y_train, y_val, 2)
    #self.plotKSTest(model, x_train, x_val, y_train, y_val, 3)




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
  epochs = 50
  batch_size = 32
  scaler_type = 'robust'
  do_early_stopping = True
  do_reduce_lr = False
  dirname = 'test'
  baseline_selection = ma_cut_wout_tv

  trainer = Trainer(
      features = features, 
      epochs = epochs,
      batch_size = batch_size,
      scaler_type = scaler_type,
      do_early_stopping = do_early_stopping,
      do_reduce_lr = do_reduce_lr,
      dirname = dirname,
      baseline_selection = baseline_selection,
      )

  trainer.process()




