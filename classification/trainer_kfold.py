import uproot
import os 
import sys
from os import path
import sys
#import psutil
import pickle
import numpy as np
import pandas as pd
import matplotlib
import uproot
import awkward as ak
matplotlib.use('pdf') #used s.t. it can be assigned to the batch


# do this before importing tf :)

# ----- For CPU usage ----
os.environ["OMP_NUM_THREADS"] = "8"          # OpenMP (used by NumPy, Eigen)
os.environ["TF_NUM_INTRAOP_THREADS"] = "8"   # TensorFlow internal threading
os.environ["TF_NUM_INTEROP_THREADS"] = "2"   # Parallel ops

# ----- For GPU usage ----
#use one gpu :) (make only one visible)
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

# more libraries
import matplotlib.pyplot as plt
from itertools import product
from time import time 
from datetime import datetime
from sklearn.preprocessing import LabelBinarizer
from tensorflow.keras import regularizers
#import shap
from contextlib import redirect_stdout
import pdb
from matplotlib.ticker import LogLocator, LogFormatterExponent

import tensorflow as tf
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
#import root_numpy
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

#to debug the train_step etc. Attention! Extremely slows down the training!!
#tf.config.run_functions_eagerly(True)
#tf.debugging.enable_check_numerics()



with open("/work/pahwagne/RDsTools/hammercpp/development_branch/weights/average_weights.yaml","r") as f:
  averages = yaml.safe_load(f)


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

# parsing
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--prod"   , required = True, help = "Specify '24' or '25' to specify the data production")
parser.add_argument("-t", "--trigger", required = True, help = "Specify 'mu7' or 'mu9' ")
parser.add_argument('-d', '--debug'  , action='store_true' )
args = parser.parse_args()

trigger = args.trigger

#shap.initjs()

now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")



########################################
#                                      #
# multiclassifier with 4 classes to    #
# optimize signal to background        #
# significance                         #
#                                      #
########################################


#------------------input files-------------------

if args.prod == "24":

  baseline_selection = base_wout_tv_24
 
  #bdt is evaluated on the skimmed datasets :) 
  base_data = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_24}/"
  file_data = [base_data + f for f in os.listdir(base_data)]
 
  #hammer is evaluated on the skimmed datasets :) 
  base_sig = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/{sig_cons_hammer_24}/" 
  file_sig = [base_sig + f for f in os.listdir(base_sig)]
  
  base_hb = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_24[0]}/"
  file_hb = [base_hb + f for f in os.listdir(base_hb)]
  
  base_bs = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bs_cons_24[0]}/"
  file_bs = [base_bs + f for f in os.listdir(base_bs)]
  
  base_b0 = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{b0_cons_24[0]}/"
  file_b0 = [base_b0 + f for f in os.listdir(base_b0)]
  
  base_bplus = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bplus_cons_24[0]}/"
  file_bplus = [base_bplus + f for f in os.listdir(base_bplus)]

elif args.prod == "25":

  baseline_selection = base_wout_tv_25

  if trigger == "mu7":
    trigger_data = " && (mu7_ip4)" #on data we select on the trigger
    trigger_mc   = "" #" && (event % 2 == 0) " #on mc we use event nr (all triggers are on for mc!)
  else:
    trigger_data = " && ((mu9_ip6) && !(mu7_ip4)) "
    trigger_mc   = " && (event % 2 == 1) " 

  #bdt is evaluated on the skimmed datasets :) 
  base_data = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_25}/"
  file_data = [base_data + f for f in os.listdir(base_data)]
  
  #hammer is evaluated on the skimmed datasets :) 
  base_sig = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{sig_cons_hammer_25}/" 
  file_sig = [base_sig + f for f in os.listdir(base_sig)]
  
  base_hb = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_25[0]}/"
  file_hb = [base_hb + f for f in os.listdir(base_hb)]
  
  base_bs = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bs_cons_25[0]}/"
  file_bs = [base_bs + f for f in os.listdir(base_bs)]
  
  base_b0 = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{b0_cons_25[0]}/"
  file_b0 = [base_b0 + f for f in os.listdir(base_b0)]
  
  base_bplus = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bplus_cons_25[0]}/"
  file_bplus = [base_bplus + f for f in os.listdir(base_bplus)]

else:
  raise ValueError ("Not a valid production year! Choose 24 or 25")

def getRdf(dateTimes, debug = None, skimmed = None, hammer = None):


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
    if skimmed and hammer:
      files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/signal_default_13_03_2025_08_42_43/*" #hammered signals 


    else:
      #access the flat ntuples
      files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*" #test

    #chain them all
    chain.Add(files)

  #create rdf from tree

  rdf = ROOT.RDataFrame(chain)

  return (chain,rdf)


def dist_corr(output, hammer, true_label, weight, scalar = None, vector = None):

  #print("======== Calculate distance correlation penalty ========")

  #output is a list of 6dim arrays, holding the score predictions.
  #hammer is a list of n-dims arrays, holding the n-variables we want to decorrrelate from the scores
  #in our case its the e{i}_up, e{i}_down, i = 1, .., 10 for the truth label.
  #true_label is a list of the true_labels.

  #pdb.set_trace()

  # Step 0. Move from one hot to single label
  true_label        = tf.argmax(true_label, axis=1)

  # STEP 1. Restrict the metric to signal events only (only have hammer variations for those!) 
 
  if scalar is not None:
    #dstau and dsmu (class 0 and class 2)
    mask              = tf.logical_or(tf.equal(true_label, 0), tf.equal(true_label, 2)) 

  elif vector is not None: 
    #dssttau and dsstmu (class 1 and class 3)
    #mask              = tf.logical_or(tf.equal(true_label, 1), tf.equal(true_label, 3)) 
    mask              = tf.equal(true_label, 3) 

  output_filt   = tf.cast(tf.boolean_mask(output, mask), tf.float32) #ensure same data type 
  hammer_filt   = tf.cast(tf.boolean_mask(hammer, mask), tf.float32) 
  weight_filt   = tf.cast(tf.boolean_mask(weight, mask), tf.float32)

  # STEP 2. Restrict the output to the corresponding scores and e{i} directions only 

  if scalar is not None:
    #dstau and dsmu (class 0 and class 2) and e1 - e6
    hammer_filt = tf.gather(hammer_filt, indices=[0], axis=1)
    scores      = tf.gather(output_filt, indices=[0, 2], axis=1)

  elif vector is not None: 
    #dssttau and dsstmu (class 1 and class 3) and e1-e10 (all indices)
    #hammer_filt = tf.gather(hammer_filt, indices=[0], axis=1)
    scores      = tf.gather(output_filt, indices=[0,1,2,3], axis=1)
    #scores       = output_filt
  #pdb.set_trace()

  # STEP 3. Calculate the two distance matrices
  score_i           = tf.expand_dims(scores, axis=1) 
  score_j           = tf.expand_dims(scores, axis=0) 
  score_diff        = score_i - score_j

  # ADD EPISLON FOR STABLE GRADIENT BACKPROPAGATION! CRASH OTHERWISE
  score_dist_mat    = tf.sqrt(tf.reduce_sum(tf.square(score_diff), axis=2) + 1e-8) 

  score_row_mean    = tf.reduce_mean(score_dist_mat, axis=1, keepdims=True) 
  score_col_mean    = tf.reduce_mean(score_dist_mat, axis=0, keepdims=True) 
  score_total_mean  = tf.reduce_mean(score_dist_mat)

  score_dist_mat    = score_dist_mat - score_row_mean - score_col_mean + score_total_mean 

  hammer_i          = tf.expand_dims(hammer_filt, axis=1) 
  hammer_j          = tf.expand_dims(hammer_filt, axis=0) 
  hammer_diff       = hammer_i - hammer_j

  # ADD EPISLON FOR STABLE GRADIENT BACKPROPAGATION! CRASH OTHERWISE
  hammer_dist_mat   = tf.sqrt(tf.reduce_sum(tf.square(hammer_diff), axis=2) + 1e-8)
 
  hammer_row_mean   = tf.reduce_mean(hammer_dist_mat, axis=1, keepdims=True) 
  hammer_col_mean   = tf.reduce_mean(hammer_dist_mat, axis=0, keepdims=True) 
  hammer_total_mean = tf.reduce_mean(hammer_dist_mat)

  hammer_dist_mat   = hammer_dist_mat - hammer_row_mean - hammer_col_mean + hammer_total_mean 

  #tf.print(" ===> score mat:", score_dist_mat, summarize=-1)
  #tf.print(" ===> hammer mat:", hammer_dist_mat, summarize=-1)

  # STEP 4. Calculate distance correlation coefficient
  size         = tf.cast(tf.shape(score_dist_mat)[0], tf.float32) # get size of matrices (quadratic, same size)
  denominator  = tf.maximum(size * size, 1e-8)
  dCov2        = tf.reduce_sum(score_dist_mat  * hammer_dist_mat) / denominator 
  dVar2_score  = tf.reduce_sum(score_dist_mat  * score_dist_mat ) / denominator 
  dVar2_hammer = tf.reduce_sum(hammer_dist_mat * hammer_dist_mat) / denominator 

  #tf.print(" ===> size :", size, summarize=-1)
  #tf.print(" ===> dCov2:", dCov2, summarize=-1)
  #tf.print(" ===> dVar2_score:", dVar2_score, summarize=-1)
  #tf.print(" ===> dVar2_hammer:", dVar2_hammer, summarize=-1)

  eps = 1e-10 #small epsilon to avoid zero div.
  
  dCorr = tf.sqrt(dCov2 + eps) / tf.sqrt(tf.sqrt(dVar2_score + eps) * tf.sqrt(dVar2_hammer + eps))

  #pdb.set_trace()
  #tf.print(" ====> Distance Correlation at this training step:", dCorr, summarize=-1)

  return dCorr 

class Losses(tf.keras.callbacks.Callback):
    def on_train_batch_end(self, batch, logs=None):
        logs = logs or {}
        total_loss = logs.get('loss')
        penalty_loss = logs.get('dist_corr_penalty')
        # main loss = total loss - penalty (approximate)
        main_loss = None
        if penalty_loss is not None and total_loss is not None:
            main_loss = total_loss - penalty_loss
        print(f"Batch {batch}: total_loss={total_loss:.4f}, penalty_loss={penalty_loss:.4f}, main_loss={main_loss:.4f}" if main_loss is not None else f"Batch {batch}: total_loss={total_loss:.4f}, penalty_loss={penalty_loss}")

class NaNOutputChecker(tf.keras.callbacks.Callback):
    def on_train_batch_end(self, batch, logs=None):
        # `self.model` is your model
        # Unfortunately, Keras callbacks don’t receive batch inputs directly,
        # so to check outputs, you need to fetch them manually (or from logs if you compute them)
        
        # One way: Run prediction on the current batch input (you need to supply inputs to the callback somehow)
        # This requires passing batch data to the callback or overriding the training loop.
        
        # Alternatively: You can check for NaNs in the loss/logs dictionary
        if logs is not None:
            loss = logs.get('loss')
            if loss is not None and tf.math.is_nan(loss):
                print(f"NaN loss detected at batch {batch}!")

    def on_test_batch_end(self, batch, logs=None):
        if logs is not None:
            loss = logs.get('loss')
            if loss is not None and tf.math.is_nan(loss):
                print(f"NaN loss detected at validation batch {batch}!")


class WeightsCheckCallback(Callback):
    def on_train_batch_end(self, batch, logs=None):
        for i, layer in enumerate(self.model.layers):
            weights = layer.get_weights()
            if weights:  # only layers with weights
                w = weights[0]  # usually the kernel weights
                # Check for NaNs
                if tf.math.reduce_any(tf.math.is_nan(w)):
                    print(f"NaN detected in weights of layer {i} ({layer.name}) at batch {batch+1}")
                else:
                    # Print summary statistics
                    print(f"Epoch {batch+1} - Layer {i} ({layer.name}): weights stats -> mean: {w.mean():.5f}, std: {w.std():.5f}, min: {w.min():.5f}, max: {w.max():.5f}")


#--------------------define sideband region (mlow, mhigh) and dump it into pickle ----------------

#chainSig,rdfSig     = getRdf(sig_cons_24) 
#chainData,rdfData   = getRdf(data_cons)
#
#
#sigma, h          = getSigma(rdfSig, "phiPi_m", base + "& (gen_sig == 0)", dt)
sigma = 0.009

#massfit = {}
#massfit["sigma"]  = sigma
##massfit["h"]      = h 
#
#
#if not args.constrained:
#  dt = datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss') 
#  A, B, C, S        = getABCS( rdfData, base , "phiPi_m", sigma, h, dt, binsFake = 21,      nSig = nSignalRegion, nSb = nSidebands, width = sbWidth)
#
#  massfit["A"] = A 
#  massfit["B"] = B
#  massfit["C"] = C
#  massfit["S"] = S
#
##write them into file to use the same constants for the evulation (for later: also use the same for pre-fit plots!)


#with open("massfit","wb") as f:
#  pickle.dump(massfit,f)

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

 
##the feature vector
kin_var = [

# helicity
'cosPhiDs_lhcb_alt',
'cosPhiDs_reco_1',
'cosPhiDs_reco_2',
#
'abs(cosPiK1)',
#'cosPiK1',

'cosMuW_lhcb_alt', 
'cosMuW_reco_1', 
'cosMuW_reco_2', 

# vertex info
'fv_chi2', # dont take prob, they can hit the prec. limit when small!
'tv_chi2', # "
'sv_chi2', # "

# displacement
'lxy_ds_sig',

# kinematics
'mu_pt',
'mu_eta',
'pi_pt',
'pi_eta',

# masses
'kk_m',
'phiPi_m',
'dsMu_m',

# delta R
'kk_deltaR',
'phiPi_deltaR',
'dsMu_deltaR',

# q2 and co
'q2_coll',

'bs_boost_lhcb_alt',
'bs_pt_lhcb_alt',
'e_star_lhcb_alt',
'm2_miss_lhcb_alt',
'pt_miss_lhcb_alt',

'bs_boost_reco_1',
'bs_pt_reco_1',
'e_star_reco_1',
'q2_reco_1',

'bs_boost_reco_2',
'bs_pt_reco_2',
'e_star_reco_2',
'q2_reco_2',

'disc_negativity',
]

hammer_var = [
"central_w",
#"e1_up",
"e2_up",
"e3_up",
#"e4_up",
"e5_up",
"e6_up",
"e7_up",
"e8_up",
#"e9_up",
#"e10_up",
#"e1_down",
#"e2_down",
#"e3_down",
#"e4_down",
#"e5_down",
#"e6_down",
#"e7_down",
#"e8_down",
#"e9_down",
#"e10_down"
]

more_vars = ["k1_pt","k2_pt","lxy_ds","mu_id_medium","fv_prob","mu_is_global",
"ds_vtx_cosine_xyz_pv","mu7_ip4","mu9_ip6",
"k1_charge", "k2_charge", "pi_charge", "mu_charge"
]


if args.prod == "24":

  #old variable names

  #displacement
  kin_var.append("ds_vtx_cosine")
  #kinematics
  kin_var.append("e_gamma")
  #isolation
  kin_var.append("mu_rel_iso_03")

if args.prod == "25":

  #old defintions
  #kin_var.append("ds_vtx_cosine_xy_pv")
  #kinematics
  #kin_var.append("e_gamma")
  #isolation
  #kin_var.append("rel_iso_03")




  #displacement
  #kin_var.append("lxy_bs_sig")
  #kin_var.append("ds_vtx_cosine_xy_pv")
  #kin_var.append("ds_vtx_cosine_xy")
  #kin_var.append("signed_decay_ip3d_mu_ds_sv")

  ##kinematics
  ##kin_var.append("e_gamma/photon_pt")
  #kin_var.append("e_gamma")
  ##kin_var.append("photon_pt")
  #kin_var.append("bs_mass_corr")
  ##kin_var.append("bs_mass_corr_photon")
  #kin_var.append("ds_perp")
  ##kin_var.append("ds_perp_photon")
  #kin_var.append("ds_mu_perp")
  ##kin_var.append("ds_mu_perp_photon")
   
  #isolation
  #kin_var.append("rel_iso_03_pv")
  print("blub")

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

    if self.signal in [0,1,2,3]:
      #signals
      #print("hammer branches")
      branches = kin_var + ['event'] + hammer_var + ["cosPhiDs_coll"] #event will be removed later!
      #branches = kin_var + ['event'] #event will be removed later!
      toLoad =  more_vars + ["gen_sig"]

    elif self.signal == 4: #hb
      branches = kin_var + ['event'] + ["cosPhiDs_coll"] #event will be removed later!
      toLoad =  more_vars + ["gen_sig", "gen_match_success"]

    else: #data
      branches = kin_var + ['event','sf_weights'] + ["cosPhiDs_coll"] #event will be removed later!
      toLoad =  more_vars



    pd_list = []
    
    if args.debug: filename = filename[:50]
    #pdb.set_trace()
    uproot_sel = self.selection.replace("&&", "&").replace("||", "|")

    start = time()
    #print (f"Loading files from: {filename}")
    #for i,ak_df in enumerate(uproot.iterate(filename + "/*.root:" + self.tree, expressions = branches , cut=uproot_sel, library="ak", step_size="1000 MB",num_workers=8)):

    #    if (i % 10 == 0): print(f"Loading file {i}, it took {round(time() - start,2)} seconds")
    #    #pdb.set_trace()
    #    df = ak.to_dataframe(ak_df[branches])
    #    # Apply selection as awkward array mask (avoid pandas query for speed)
    #    #mask = eval(selection, {}, arrays)  # this only works if selection is valid python expression
    #    #filtered = arrays[mask]
    #    # Convert to pandas only after filtering
    #    #df = arrays[branches].to_pandas()
    #    pd_list.append(df)

    for name in filename:

      f = ROOT.TFile(name)
      tree_obj = f.Get(self.tree)
      arr = tree2array(tree_obj,selection = self.selection, branches = branches) #event will be removed later!
      pd_list.append(pd.DataFrame(arr))

    #  f    = uproot.open(name)
    #  tree = f[self.tree]
    #  arrays = tree.arrays(branches + toLoad, library="np")
    #  df = pd.DataFrame(arrays)
    #  uproot_sel = self.selection.replace("&&", "and").replace("||", "or")
    #  df = df.query(uproot_sel)
    #  df = df[branches]
    #  pd_list.append(df)

    #  #f = ROOT.TFile(name)
    #  #tree_obj = f.Get(self.tree)
    #  #arr = tree2array(tree_obj,selection = self.selection, branches = branches) #event will be removed later!
    #  #pd_list.append(pd.DataFrame(arr))

    self.df = pd.concat(pd_list)

    pd.options.mode.chained_assignment = None
    self.df['is_signal'] = signal


class Trainer(object):

  'train the data'

  def __init__(self, features, epochs, batch_size, learning_rate, lambda_penalty_scalar, lambda_penalty_vector, scaler_type, es_patience, do_reduce_lr, dirname, baseline_selection, nfolds, frac_sb, frac_sf):
    self.features           = features #variables to train on
    self.epochs             = epochs #samples / batch_size =  number of iterations to 1 epoch
    self.batch_size         = batch_size 
    self.learning_rate      = learning_rate 
    self.lambda_penalty_scalar     = lambda_penalty_scalar 
    self.lambda_penalty_vector     = lambda_penalty_vector 
    self.scaler_type        = scaler_type
    self.es_patience        = es_patience 
    self.do_reduce_lr       = do_reduce_lr
    self.dirname            = dirname + '_' + datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss')
    self.baseline_selection = baseline_selection
    self.nfolds             = nfolds
    self.frac_sb            = frac_sb
    self.frac_sf            = frac_sf
    self.colors             = [
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
    f.write("Early stopping patience: "   + str(self.es_patience)         + "\n")
    f.write("Reduce lr: "                 + str(self.do_reduce_lr)        + "\n")
    f.write("Baseline selection: "        + str(self.baseline_selection)  + "\n")
    f.write("penalty scalar: "            + str(self.lambda_penalty_scalar)      + "\n")
    f.write("penalty vector: "            + str(self.lambda_penalty_vector)      + "\n")
    f.write("year of data: "              + str(args.prod)                + "\n")
    f.write("trigger: "                   + str(args.trigger)             + "\n")
    f.write("Signal region up to:"        + str(nSignalRegion) + " sigma" + "\n")
    f.write("Sideband region starts at:"  + str(nSidebands)    + " sigma" + "\n")
    f.write("Sideband region width:"      + str(sbWidth)       + " sigma" + "\n")
    #f.write("frac of sideband events:"    + str(self.frac_sb)             + "\n")
    #f.write("frac of signflip events:"    + str(self.frac_sf)             + "\n")
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
    hb_selec = " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) && (gen_match_success == 1)"
 
    #signals      #ds mu   #ds tau   #dstar mu   #dstar tau   #hb
    mc_ids      = [0       ,1        ,10         ,11          ,-1]
    mc_classes  = [2       ,0        ,3          ,1           ,4 ] #changed ! 
    label       = {0: r" \mu", 1:r" \tau", 10:r" \mu^{*}", 11: r" \tau^{*}", -1: r" H_{b}" }
    class_label = {class_id: [] for class_id in mc_classes}

    for mc_id, class_id in zip (mc_ids,mc_classes):

      class_label[class_id].append(label[mc_id]) 

      if mc_id >= 0: 
        mc_sample            = Sample(filename = file_sig,    selection=self.baseline_selection  + trigger_mc + f'&& (gen_sig == {mc_id})' + signalRegion + lowMass,                tree = tree_name,signal = class_id).df
        #mc_sample            = Sample(filename = base_sig,     selection=self.baseline_selection  + trigger_mc + f'&& (gen_sig == {mc_id})' + signalRegion + lowMass,                tree = tree_name,signal = class_id).df

      else:
        #mc_sample            = Sample(filename = file_hb,     selection=self.baseline_selection  + trigger_mc + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df

        mc_sample1          = Sample(filename = file_b0,     selection=self.baseline_selection  + trigger_mc + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df
        mc_sample2          = Sample(filename = file_bs,     selection=self.baseline_selection  + trigger_mc + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df
        mc_sample3          = Sample(filename = file_bplus,  selection=self.baseline_selection  + trigger_mc + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df

        #mc_sample1          = Sample(filename = base_b0,     selection=self.baseline_selection  + trigger_mc + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df
        #mc_sample2          = Sample(filename = base_bs,     selection=self.baseline_selection  + trigger_mc + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df
        #mc_sample3          = Sample(filename = base_bplus,  selection=self.baseline_selection  + trigger_mc + hb_selec                + signalRegion + lowMass,                tree = tree_name,signal = class_id).df

        mc_sample           = pd.concat([mc_sample1,mc_sample2,mc_sample3], sort = False)
      
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
    data_selec  = self.baseline_selection + sign_flip + signalRegion + lowMass

    data_sample_sf = Sample(filename=file_data,   selection=data_selec + trigger_data , tree = 'tree',signal = data_id).df
    #data_sample_sf = Sample(filename=base_data,   selection=data_selec + trigger_data , tree = 'tree',signal = data_id).df
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
    # pd which don't have a "central_w" column will just carry a NAN in this column
    train = pd.concat([train_samples[key] for key in train_samples.keys()], sort = False)
    test  = pd.concat([test_samples[key] for key  in test_samples.keys()],  sort = False)

    # since our dataset is cleaned from nans (except the one we get from concatenating),
    # we can replace them with 1.0. Like this, all non-hammered events have simply weight 1
    # same for sf_weight, it will be 1.0 for MC
    #pdb.set_trace()
    train = train.fillna(1.0)
    test  = test .fillna(1.0)

    # create new column which holds the multiplication of all weights

    #trick

    averages["central_w_hb"]     = 1.0
    averages["central_w_comb"]   = 1.0

    #strip only the directions e1up, ... from hammer
    hammer_dir = [h for h in hammer_var if h != "central_w"]

    for e in hammer_dir:
      averages[f"{e}_hb"]    = 1.0
      averages[f"{e}_comb"]  = 1.0


    keys = {0: 'dstau', 1: 'dsstartau', 2:'dsmu' , 3: 'dsstarmu', 4:'hb', 5:'comb'}

    #for every event, pick the correct central_av accoding to is_signal
    central_av_train = [averages["central_w_" + keys[sig]] for sig in train["is_signal"]]
    central_av_test  = [averages["central_w_" + keys[sig]] for sig in test ["is_signal"]]

    #REMARK: the entries e7-e7 are nonsense for dsmu and dstau here!

    #for every event, pick the correct variational av according to signal
    var_av_train     = {}
    var_av_test      = {}

    for e in hammer_dir:

      var_av_train[f"average_{e}"] = [averages[f"{e}_" + keys[sig]] for sig in train["is_signal"]]
      var_av_test [f"average_{e}"] = [averages[f"{e}_" + keys[sig]] for sig in test ["is_signal"]]

      #normalize weights
      train[e] =  train[e] /  var_av_train[f"average_{e}"]
      test [e] =  test [e] /  var_av_test [f"average_{e}"]


    #pdb.set_trace()
    train["total_w"] = train["sf_weights"] * train["central_w"] / central_av_train
    test ["total_w"] = test ["sf_weights"] * test ["central_w"] / central_av_test
    

    # even undo the splitting which is not used for k-folding
    main_df = pd.concat([train,test],sort= False)
    #pdb.set_trace()
    # re-index (this does not shuffle events, only reindex!)
    main_df.index = np.array(range(len(main_df)))

    # shuffle
    main_df = main_df.sample(frac=1, replace=False, random_state=5366) # of course, keep R's seed ;)

    # X and Y, keep event number for kfold splitting! will be removed later!

    #w_cols = ['central_w', 'event']
    #w_cols = ['sf_weights', 'central_w', 'event']
    w_cols   = ['total_w', 'event']
    #hammer_var_cols = hammer_dir + ["event"]
    hammer_var_cols = ["cosPhiDs_coll", "event"]

    X       = pd.DataFrame(main_df, columns=list(set(self.features)) + ['event'] )
    Y       = pd.DataFrame(main_df, columns=['is_signal', 'event'])
    weights = pd.DataFrame(main_df, columns=w_cols)
    hammer  = pd.DataFrame(main_df, columns=hammer_var_cols)

    # save the exact list of features
    features_filename = '/'.join([self.outdir, 'input_features.pck'])
    pickle.dump(self.features, open(features_filename, 'wb' ))
    print('========> {} created'.format(features_filename))

    x_folds = {}
    y_folds = {}
    w_folds = {}
    h_folds = {}

    xx_folds = {}
    for n in range(self.nfolds):
      
      #trick to split into folds using event number! :D
      x_folds[n] = X      [ X["event"]       % self.nfolds == n ]
      y_folds[n] = Y      [ Y["event"]       % self.nfolds == n ]
      w_folds[n] = weights[ weights["event"] % self.nfolds == n ]
      h_folds[n] = hammer [ hammer ["event"] % self.nfolds == n ]

      #delete event columns afterwards!
      x_folds[n] = x_folds[n].drop("event", axis = 1)
      y_folds[n] = y_folds[n].drop("event", axis = 1)
      w_folds[n] = w_folds[n].drop("event", axis = 1)
      h_folds[n] = h_folds[n].drop("event", axis = 1)


      # scale the features
      # this is an important step!
  
      xx_folds[n], qt =  self.doScaling(x_folds[n])

      # and save the scaler, which will have to be used throughout the full process, even at evaluation time
      scaler_filename = '/'.join([self.outdir, f'input_tranformation_weighted_fold_{n}.pck'])
      pickle.dump(qt,open(scaler_filename, 'wb'))
      print( ' ========> {} created'.format(scaler_filename))
 
 
    return main_df, qt, xx_folds, y_folds, w_folds, h_folds

 
    

  def defineModel(self):
    '''
      Define the NN
    '''

    # params
    l2_rate = 0.01

    model = tf.keras.Sequential()
    # features
    model.add(tf.keras.layers.Input((len(features),)))
    # body
    model.add(tf.keras.layers.Dense(64 ,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(128,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(256,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    #odel.add(tf.keras.layers.Dense(256,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    #model.add(tf.keras.layers.Dense(512,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    #model.add(tf.keras.layers.Dense(256,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(256,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(128,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    model.add(tf.keras.layers.Dense(64 ,   activation ='swish',   kernel_regularizer=regularizers.l2(l2_rate)))
    # output
    model.add(tf.keras.layers.Dense(6  ,   activation= 'softmax'))

    # optimizer
    opt = keras.optimizers.Adam(learning_rate=self.learning_rate)

    model.compile(optimizer=opt, loss='categorical_crossentropy', metrics=['acc'])

    # print model

    print(model.summary()) #not enough

    def print_model_summary(model):

      print("============= DETAILED MODEL SUMMARY ===============")
      print("learning rate = ", self.learning_rate)
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

    # write model into file
    f = open(self.outdir + "/settings.txt", "a")
    with redirect_stdout(f): 
      model.summary()
      print_model_summary(model)
    f.close()
  
    return model


  def prepareInputs(self, xx_folds, y_folds, w_folds, h_folds):
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
      y_folds[n]  = y_folds[n] .reset_index(drop=True)
      w_folds[n]  = w_folds[n] .reset_index(drop=True)
      h_folds[n]  = h_folds[n] .reset_index(drop=True)

    return xx_folds, y_folds, w_folds, h_folds


  def train(self, xx_folds, y_folds, w_folds, h_folds, weight, class_label):
    '''
      Perform the training
    '''

    xx_train = {}
    y_train  = {}
    w_train  = {}
    h_train  = {}

    xx_val   = {}
    y_val    = {}
    w_val    = {}
    h_val    = {}

    model    = {}
    history  = {}

    # create cyclic list, f.e. for nfold = 5, this is:
    # [0,1,2,3,4,0,1,2,3,4]
    cyclic = list(range(self.nfolds)) * 2  

    for n in range(self.nfolds):

      train_indices = [cyclic[n + (i+1)] for i in range(self.nfolds - 2)  ]
      test_index    = cyclic[ n + self.nfolds - 1]

      print(f"For fold/net {n} we train on folds {train_indices} and test on fold {test_index}")

      #important to re-index after concatinating several pandas array! drop = True drops the old numbering 
      xx_train[n] = pd.concat([ xx_folds[i] for i in train_indices]).reset_index(drop=True)
      y_train[n]  = pd.concat([ y_folds[i]  for i in train_indices]).reset_index(drop=True)
      w_train[n]  = pd.concat([ w_folds[i]  for i in train_indices]).reset_index(drop=True)
      h_train[n]  = pd.concat([ h_folds[i]  for i in train_indices]).reset_index(drop=True)
     
      xx_val[n]   = xx_folds[test_index] .reset_index(drop=True)
      y_val [n]   = y_folds[test_index]  .reset_index(drop=True)
      w_val [n]   = w_folds[test_index]  .reset_index(drop=True)
      h_val [n]   = h_folds[test_index]  .reset_index(drop=True)


      # calculate the class weight (not the same as sample weight!
      # this corrects for class imbalance including the effect of the sample
      # weights
      eff_pop = np.zeros(6)

      for i in range(0,6):

        # effective population of classes 0-5 for fold n:
        eff_pop[i] = w_train[n][y_train[n]['is_signal'] == i].sum()

      total_pop = np.sum(eff_pop)
      class_w   = {i: total_pop / eff_pop[i] for i in range (0,6)} 

      #now, absorb the class weight into the sample weights as well :)

      #pdb.set_trace()
                                          #this gives the class 0, ... 5
      w_train[n]["total_w"] = [ class_w[  y_train[n]['is_signal'][ev]    ] * w for ev,w in enumerate(w_train[n]["total_w"]) ]
      w_val  [n]["total_w"] = [ class_w[  y_val  [n]['is_signal'][ev]    ] * w for ev,w in enumerate(w_val  [n]["total_w"]) ]
     

      #convert to numpy for training (and use one-hot for y)
      xx_train[n] = xx_train[n].to_numpy()
      xx_val[n]   = xx_val[n].to_numpy()


      # not one hot encoded + convert to numpy (used for penalty term)
      y_true_label_train = y_train[n]['is_signal'].to_numpy()
      y_true_label_val   = y_val[n]['is_signal'].to_numpy()

      # encode + convert to numpy
      y_train[n] = tf.one_hot(y_train[n]['is_signal'].to_numpy(), 6)
      y_train[n] = y_train[n].numpy()

      y_val[n] = tf.one_hot(y_val[n]['is_signal'].to_numpy(), 6)
      y_val[n] = y_val[n].numpy()

      # convert to numpy
      w_train[n] = w_train[n].to_numpy()
      w_val[n]   = w_val[n].to_numpy()

      h_train[n] = h_train[n].to_numpy()
      h_val[n]   = h_val[n].to_numpy()


      #clear model and redefine a new one for every fold! 

      from tensorflow.keras import backend as K
      K.clear_session()
      model[n] = self.defineModel()
      print(" ========> Received model")

      ##########################################
      # prepare data for custom train/val loop #
      ##########################################
      
      batch_size = self.batch_size
      epochs = self.epochs
      lambda_penalty_scalar = self.lambda_penalty_scalar 
      lambda_penalty_vector = self.lambda_penalty_vector 
      optimizer = tf.keras.optimizers.Adam(learning_rate=self.learning_rate)
     
      #define accuracy
      acc_train = tf.keras.metrics.CategoricalAccuracy()     
      acc_val   = tf.keras.metrics.CategoricalAccuracy()     
 
      #########################################
      # Implement custom train/val loop       #
      #########################################
       
      @tf.function
      def train_step(x, y, w, h):
          with tf.GradientTape() as tape:

              # Get score predictions (outputs)
              outputs = model[n](x, training=True)
              #main_loss = loss_main(y, outputs, sample_weight=w)
              main_loss = model[n].compiled_loss(y, outputs, regularization_losses=model[n].losses, sample_weight=w)
              #print(" ====> Get crossentropy loss", main_loss)
              penalty_scalar = 0.0 #dist_corr(outputs, h , y, scalar = True)  
              penalty_vector = dist_corr(outputs, h , y, w, vector = True)  
              total_loss = main_loss      +  \
                           lambda_penalty_scalar * penalty_scalar +  \
                           lambda_penalty_vector * penalty_vector
              #total_loss = lambda_penalty_vector * penalty_vector
              #print(" ====> Get total loss", total_loss)
              # Get accuracy and update it
              acc_train.update_state(y, outputs, sample_weight=w)              

          # calc gradients
          gradients = tape.gradient(total_loss, model[n].trainable_variables)

          # apply gradients
          optimizer.apply_gradients(zip(gradients, model[n].trainable_variables))
      
          return total_loss, main_loss, penalty_scalar, penalty_vector, gradients #, grad_norms, clipped_grad_norms
     


      @tf.function
      def val_step(x, y, w, h):
          with tf.GradientTape() as tape:

              # Get score predictions (outputs)
              outputs = model[n](x, training=False)
              #main_loss = loss_main(y, outputs, sample_weight=w)
              main_loss = model[n].compiled_loss(y, outputs, regularization_losses=model[n].losses, sample_weight=w)
              #print(" ====> Get crossentropy loss", main_loss)
              penalty_scalar = 0.0 #dist_corr(outputs, h , y, scalar = True)  
              penalty_vector = dist_corr(outputs, h , y, w, vector = True)  
              total_loss = main_loss      +   \
                           lambda_penalty_scalar * penalty_scalar +  \
                           lambda_penalty_vector * penalty_vector
              #total_loss = lambda_penalty_vector * penalty_vector
              #print(" ====> Get total loss", total_loss)
              # Get accuracy and update it
              acc_val.update_state(y, outputs, sample_weight=w)              
      
          return total_loss, main_loss, penalty_scalar, penalty_vector


      history[n] = {
        "xentropy_loss_train" :[], 
        "total_loss_train"    :[], 
        "penalty_scalar_train":[],
        "penalty_vector_train":[],
        "acc_train"           :[],
        "xentropy_loss_val"   :[],
        "total_loss_val"      :[],
        "penalty_scalar_val"  :[],
        "penalty_vector_val"  :[],
        "acc_val"             :[],
        }


      #create tf datasets
      tf_train_dataset    =  tf.data.Dataset.from_tensor_slices((xx_train[n], y_train[n], w_train[n], h_train[n]))
      tf_val_dataset      =  tf.data.Dataset.from_tensor_slices((xx_val[n], y_val[n], w_val[n], h_val[n]))        

      #set the counter for early stopping
      es_counter = 0
       
      for i, epoch in enumerate(range(self.epochs)):

        print(f" ========> At epoch {epoch+1}/{epochs}")

        start = time()

        # divide into batches 
        # for every epoch, shuffle events differently into the mini batches, to avoid pattern learning. 
        # It does not disturb the relatice order between x and y, i checked. It simply shuffles events.
        tf_train    =  tf_train_dataset.shuffle(buffer_size=1000).batch(self.batch_size).prefetch(tf.data.AUTOTUNE)
        #no need to shuffle validation :)
        tf_val      =  tf_val_dataset                            .batch(self.batch_size).prefetch(tf.data.AUTOTUNE)

        print("Time to shuffle:", time() - start)

        #prepare metrics for train and test
        xentropy_loss_train = tf.keras.metrics.Mean()
        penalty_scalar_train       = tf.keras.metrics.Mean()
        penalty_vector_train       = tf.keras.metrics.Mean()
        total_loss_train    = tf.keras.metrics.Mean()

        xentropy_loss_val   = tf.keras.metrics.Mean()
        penalty_scalar_val         = tf.keras.metrics.Mean()
        penalty_vector_val         = tf.keras.metrics.Mean()
        total_loss_val      = tf.keras.metrics.Mean()

        for batch, (x, y, w, h) in enumerate(tf_train):

          #print(f" ====> At batch {batch}/{len(tf_train)}")
          #total_loss, xentropy_loss, penalty, grad_norms, clipped_grad_norms = train_step(x, y, w, h)
          total_loss, xentropy_loss, penalty_scalar, penalty_vector, gradients = train_step(x, y, w, h)

          #print("Gradient norms before clipping:", gradients)
          #print("Gradient norms after clipping:", clipped_grad_norms)

          #tf.print(f"Step {step} Total Loss:", total_loss, "Penalty:", penalty)

          total_loss_train   .update_state(total_loss   )
          xentropy_loss_train.update_state(xentropy_loss)
          penalty_scalar_train      .update_state(penalty_scalar      )
          penalty_vector_train      .update_state(penalty_vector      )

        print("Time to train:", time() - start)
        print("Penalty:", penalty_vector)


        for batch, (x, y, w, h) in enumerate(tf_val):

          #print(f" ====> At batch {batch}/{len(tf_val)}")
          total_loss, xentropy_loss, penalty, penalty = val_step(x, y, w, h)
          #tf.print(f"Step {step} Total Loss:", total_loss, "Penalty:", penalty)

          total_loss_val   .update_state(total_loss   )
          xentropy_loss_val.update_state(xentropy_loss)
          penalty_scalar_val      .update_state(penalty_scalar      )
          penalty_vector_val      .update_state(penalty_vector      )

        print("Time to evaluate:", time() - start)

        #after one epoch, take the average of all metrics (loss, acc) by calling .result()

        #save validation acc/loss in variables, as we need them later for early stopping and model saving
        loss_val_now   = total_loss_val.result().numpy()
        acc_val_now    = acc_val       .result().numpy()
 

        history[n]["total_loss_train"   ].append( total_loss_train   .result().numpy())
        history[n]["xentropy_loss_train"].append( xentropy_loss_train.result().numpy())
        history[n]["penalty_scalar_train"      ].append( penalty_scalar_train      .result().numpy())     
        history[n]["penalty_vector_train"      ].append( penalty_vector_train      .result().numpy())     
        history[n]["acc_train"          ].append( acc_train          .result().numpy())     

        history[n]["total_loss_val"     ].append( loss_val_now                      )
        history[n]["xentropy_loss_val"  ].append( xentropy_loss_val.result().numpy())
        history[n]["penalty_scalar_val"        ].append( penalty_scalar_val      .result().numpy())     
        history[n]["penalty_vector_val"        ].append( penalty_vector_val      .result().numpy())     
        history[n]["acc_val"            ].append( acc_val_now                       )     


        #implement early stopping
        if (i == 0): 
          #first epoch
          best_loss = loss_val_now
        elif (loss_val_now >= best_loss):
          #if validation loss is bigger or equal (i.e. no improvement, increase the counter)
          es_counter += 1
          print("Increasing counter! Best loss is: ", best_loss, "while current loss is: ", loss_val_now)
        else: 
          # if we found an epoch that decreased the loss, reset counter to 0 and assign new best loss
          es_counter = 0
          best_loss = loss_val_now

        #print("Best loss is:", best_loss)
        #at every epoch, check if the counter exceeds the patience
        if es_counter > self.es_patience:
          break; 

        #save the model 
        filepath = '/'.join([
          self.outdir, 
          f'fold_{n}' + f'_saved-model-{epoch:04d}_val_loss_{total_loss_val.result().numpy():.4f}_val_acc_{acc_val.result().numpy():.4f}.h5'
        ])

        #first epoch assign as best acc
        if (i == 0): 
          best_acc = acc_val_now 
          model[n].save(filepath)

        #if current acc better than best acc
        if ( acc_val_now > best_acc): 
          #save model and set new best value:
          model[n].save(filepath)
          best_acc = acc_val_now 


        #reset accuracy for next epoch, to calculate fresh averages
        # loss mean is automatically reset (local object)
        acc_train.reset_states()
        acc_val  .reset_states()
      
 
        print("Epoch time:", time() - start)



      #plot each fold immediately

      self.plotMetric(history, "xentropy_loss_train","Training X-Entropy Loss"                                                     , "Loss"   , fold = n)
      self.plotMetric(history, "total_loss_train"   ,"Training Total Loss"                                                         , "Loss"   , fold = n)
      self.plotMetric(history, "penalty_scalar_train"      ,r"Training Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_scalar}" , "Penalty", fold = n)
      self.plotMetric(history, "penalty_vector_train"      ,r"Training Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_vector}" , "Penalty", fold = n)
      self.plotMetric(history, "acc_train"          ,"Training Accuracy"                                                           , "Acc."   , fold = n)
      self.plotMetric(history, "xentropy_loss_val","Validation X-Entropy Loss"                                                     , "Loss"   , fold = n)
      self.plotMetric(history, "total_loss_val"   ,"Validation Total Loss"                                                         , "Loss"   , fold = n)
      self.plotMetric(history, "penalty_scalar_val"      ,r"Validation Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_scalar}" , "Penalty", fold = n)
      self.plotMetric(history, "penalty_vector_val"      ,r"Validation Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_vector}" , "Penalty", fold = n)
      self.plotMetric(history, "acc_val"          ,"Validation Accuracy"                                                           , "Acc."   , fold = n)


      self.plotScore (model, xx_train, y_train, xx_val, y_val, 0,class_label, fold = n)
      self.plotScore (model, xx_train, y_train, xx_val, y_val, 1,class_label, fold = n)
      self.plotScore (model, xx_train, y_train, xx_val, y_val, 2,class_label, fold = n)
      self.plotScore (model, xx_train, y_train, xx_val, y_val, 3,class_label, fold = n)
      self.plotScore (model, xx_train, y_train, xx_val, y_val, 4,class_label, fold = n)
      self.plotScore (model, xx_train, y_train, xx_val, y_val, 5,class_label, fold = n)
      self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 0, fold = n)
      self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 1, fold = n)
      self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 2, fold = n)
      self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 3, fold = n)
      self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 4, fold = n)
      self.plotKSTest(model, xx_train, y_train, xx_val, y_val, 5, fold = n)




    #return model, history, xx_train, y_train, h_train, xx_val, y_val, h_val
    return history, model, xx_train, y_train, h_train, xx_val, y_val, h_val


  def plotMetric(self, history, history_key, title, ylabel, fold = None):
    '''
      Plot the loss/acc metric for training and validation sets
    '''

    # folds can have different length due to early stopping!!!

    ############
    # Training #
    ############


    if fold is None:

      # TOTAL LOSS
     
      average = []
      #max epochs (get max epochs)
      max_epochs = max(len(history[n][history_key]) for n in range(self.nfolds))
  
      for n in range(self.nfolds):
  
        toPlot = history[n][history_key]
        epochs = range(1, len(toPlot)+1)
        plt.plot(epochs, toPlot, self.colors[n], label=f'Fold {n}')
  
        #create a full array of nans of length max_epochs       
        padded = np.full(max_epochs, np.nan)
        #fill the first part with (the rest stays nan)
        padded[:len(toPlot)] = toPlot 
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
      self.saveFig(plt, history_key.replace(" ", "_"))
      self.saveFig(plt, history_key.replace(" ", "_"))
  
      #log
   
      plt.yscale('log')  
      #plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
      self.saveFig(plt, history_key.replace(" ", "_") + "_log")
      self.saveFig(plt, history_key.replace(" ", "_") + "_log")
      #plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))
      plt.clf()
      plt.close()


    else: 

      n = fold
      # TOTAL LOSS
  
      toPlot = history[n][history_key]
      epochs = range(1, len(toPlot)+1)
      plt.plot(epochs, toPlot, self.colors[n], label=f'Fold {n}')
  
      plt.subplots_adjust(left=0.2, right=0.95, top=0.85, bottom=0.20)
  
      plt.title(title)
      plt.xlabel('Epochs')
      plt.ylabel(ylabel)
      plt.legend()
      self.saveFig(plt, history_key.replace(" ", "_") + f"_fold{n}")
      self.saveFig(plt, history_key.replace(" ", "_") + f"_fold{n}")
  
      #log
   
      plt.yscale('log')  
      #plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
      self.saveFig(plt, history_key.replace(" ", "_") + "_log" + f"_fold{n}")
      self.saveFig(plt, history_key.replace(" ", "_") + "_log" + f"_fold{n}")
      #plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))
      plt.clf()
      plt.close()


    ##validation
    #average_loss = []
    #for n in range(self.nfolds):

    #  loss_val = history[n].history['val_loss']
    #  epochs = range(1, len(loss_val)+1)
    #  plt.plot(epochs, loss_val, self.colors[n], label=f'Fold {n}')

    #  #create a full array of nans of length max_epochs       
    #  padded = np.full(max_epochs, np.nan)
    #  #fill the first part with the loss (the rest stays nan)
    #  padded[:len(loss_val)] = loss_val
    #  #append the padded array to the average loss
    #  average_loss.append(padded)


    #average_loss = np.nanmean(np.vstack(average_loss), axis=0)
    #plt.plot(range(1, max_epochs + 1), average_loss, 'black' , label='Average loss', linestyle= 'dashed')
    #plt.title('Validation Loss')
    #plt.xlabel('Epochs')
    #plt.ylabel('Loss')
    #plt.legend()
    #self.saveFig(plt, 'validation_loss')
    #plt.clf()

    ##log
    #average_loss = []
    #for n in range(self.nfolds):

    #  loss_val = history[n].history['val_loss']
    #  epochs = range(1, len(loss_val)+1)
    #  plt.plot(epochs, loss_val, self.colors[n], label=f'Fold {n}')

    #  #create a full array of nans of length max_epochs       
    #  padded = np.full(max_epochs, np.nan)
    #  #fill the first part with the loss (the rest stays nan)
    #  padded[:len(loss_val)] = loss_val
    #  #append the padded array to the average loss
    #  average_loss.append(padded)


    #average_loss = np.nanmean(np.vstack(average_loss), axis=0)
    #plt.plot(range(1, max_epochs + 1), average_loss, 'black' , label='Average loss', linestyle= 'dashed')
    #plt.title('Validation Loss')
    #plt.xlabel('Epochs')
    #plt.ylabel('Loss')
    #plt.yscale('log')  
    #plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=None, numticks=10))
    #plt.gca().yaxis.set_major_formatter(LogFormatterExponent(base=10.0))
    #plt.legend()
    #self.saveFig(plt, 'validation_loss_log')
    #plt.clf()


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

  def plotScoreTauOnly(self, model, xx_train, y_train, h_train, xx_val, y_val, h_val, class_label):
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
        h_chan    = h_val[n] [ true_1d == chan ]
        #...and predict their score! (data is already scaled!)
        y_chan    = model[n].predict([x_chan, h_chan])[:,sig] 

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

   
  def plotScore(self, model, xx_train, y_train, xx_val, y_val, sig, class_label, fold = None):
    '''
      Plot the score of all channels class sig
    '''

    channels = [0,1,2,3,4,5]
    col ={0:'c',1:'g',2:'m',3:'r',4:'b',5:'k'}
    linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}

    hist_content = {}


    if fold is None:

      #plot average
      for chan in channels:
  
        dummy = []
  
        for n in range(self.nfolds):
       
          print("chan is ", chan, ", n is ", n, ", sig is ", sig) 
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
      plt.close()


    else:

      fig = plt.figure()
      #plot average
      for chan in channels:
  
        dummy = []
        n = fold 
       
        print("chan is ", chan, ", n is ", n, ", sig is ", sig) 
        #class predictions 1D
        #is of shape [1,2,5,4,3,2,4,2,1, ....] 
        true_1d = np.argmax(y_val[n], axis=1)
  
        #select only events where the true class is chan...
        x_chan    = xx_val[n][ true_1d == chan ]
        #...and predict their score! (data is already scaled!)
        y_chan    = model[n].predict(x_chan)[:,sig] 
  
        # Plot data into histo
        hist = plt.hist(y_chan, bins=np.arange(0,1.025,0.025), color=col[chan], alpha=0.5, label=class_label[chan], histtype='stepfilled',density = True,linestyle = linestyles[chan], linewidth = 1.5)

      plt.legend(loc='upper right')
      plt.title(f'Score for class ' + class_label[sig] + " (Average over folds)")
      plt.xlabel('score')
      plt.ylabel('events')
      #fig.savefig('outputs/score_' + str(sig) + '.pdf')
      #fig.savefig('outputs/score_' + str(sig) + '.png')
      self.saveFig(fig, f'score_' + str(sig) + f"_fold_{fold}" )
      plt.clf()
      plt.close()


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

  
  def plotKSTest(self, model, xx_train, y_train, xx_val, y_val,sig, fold = None):
    '''
      Plot the outcome of the Kolmogorov test
      Used to test the overfitting
    '''

    h1 = ROOT.TH1F(f'train_{sig}', f'train_{sig}', 30, 0, 1)
    h2 = ROOT.TH1F(f'val_{sig}', f'val_{sig}', 30, 0, 1)



    if fold is None:

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
      ks_value.AddText(f'Average KS score of class {sig} = {round(ks_score,3)}')
      ks_value.SetFillColor(0)
      ks_value.Draw('EP SAME')
  
      leg = ROOT.TLegend(.55,.82,.83,.88)
      leg.AddEntry(h1 ,'Training' ,'F' )
      leg.AddEntry(h2 ,'Validation' ,'EP' )
      leg.Draw("SAME")
  
  
      c1.SaveAs(self.outdir + f'/KS_test_{sig}.pdf')
      c1.SaveAs(self.outdir + f'/KS_test_{sig}.png')
      #print('KS score: ',ks_score, len(train_pred),len(test_pred))
      c1.Close()
      del c1
  

    else:

      n = fold

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
 

      c1=ROOT.TCanvas()
      if h1_dummy.Integral()!=0: h1_dummy.Scale(1./h1_dummy.Integral())
      if h2_dummy.Integral()!=0: h2_dummy.Scale(1./h2_dummy.Integral())
      c1.Draw()
      h1_dummy.GetXaxis().SetTitle("score")
      h1_dummy.GetYaxis().SetRangeUser(0, max([h2_dummy.GetMaximum(),h1_dummy.GetMaximum()])*1.6)
  
      h1_dummy.SetFillColor(ROOT.kGreen+2)
      h1_dummy.SetLineColor(ROOT.kGreen+2)
      h1_dummy.SetFillStyle(3345)
      h1_dummy.Draw('HIST')
      
      h2_dummy.SetLineColor(ROOT.kBlack)
      h2_dummy.SetMarkerStyle(8)
      h2_dummy.SetMarkerSize(0.5)
      h2_dummy.SetMarkerColor(ROOT.kBlack)
      h2_dummy.Draw('EP SAME')
  
      ks_score = h1_dummy.KolmogorovTest(h2_dummy)
      ks_value = ROOT.TPaveText(0.5, 0.76, 0.88, 0.80, 'nbNDC')
      ks_value.AddText(f'Average KS score of class {sig} = {round(ks_score,3)}')
      ks_value.SetFillColor(0)
      ks_value.Draw('EP SAME')
  
      leg = ROOT.TLegend(.55,.82,.83,.88)
      leg.AddEntry(h1_dummy ,'Training' ,'F' )
      leg.AddEntry(h2_dummy ,'Validation' ,'EP' )
      leg.Draw("SAME")
  
  
      c1.SaveAs(self.outdir + f'/KS_test_{sig}_fold_{n}.pdf')
      c1.SaveAs(self.outdir + f'/KS_test_{sig}_fold_{n}.png')
      #print('KS score: ',ks_score, len(train_pred),len(test_pred))
      c1.Close()
      del c1











  def process(self):
    print( '------------------- MVA Trainer --------------------')
 
    global xx_folds
    global y_folds
    global w_folds
    global h_folds
    global weight
    global class_label
   
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
    main_df, qt, xx_folds, y_folds, w_folds, h_folds = self.preprocessing(train_notnan, test_notnan)


    print('\n========> preparing the inputs') 

    xx_folds, y_folds, w_folds, h_folds  = self.prepareInputs(xx_folds, y_folds, w_folds, h_folds)

    # do the training
    print('\n========> training...') 
    #model, history, xx_train, y_train, h_train, xx_val, y_val, h_val = self.train(xx_folds, y_folds, w_folds, h_folds, weight)
    history, model, xx_train, y_train, h_train, xx_val, y_val, h_val = self.train(xx_folds, y_folds, w_folds, h_folds, weight, class_label)

    #save the output dictionaries
    with open(f'{self.outdir}/xx_val.pck', 'wb') as f:
      pickle.dump(xx_val,f)

    with open(f'{self.outdir}/xx_train.pck', 'wb') as f:
      pickle.dump(xx_train,f)

    with open(f'{self.outdir}/y_val.pck', 'wb') as f:
      pickle.dump(y_val,f)

    with open(f'{self.outdir}/y_train.pck', 'wb') as f:
      pickle.dump(y_train,f)

    #pdb.set_trace()

    # plotting
    print('\n========> plotting...' )
                             #dict key             # fig title                                                                   # y axis 
    self.plotMetric(history, "xentropy_loss_train","Training X-Entropy Loss"                                                     , "Loss"   )
    self.plotMetric(history, "total_loss_train"   ,"Training Total Loss"                                                         , "Loss"   )
    self.plotMetric(history, "penalty_scalar_train"      ,r"Training Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_scalar}" , "Penalty")
    self.plotMetric(history, "penalty_vector_train"      ,r"Training Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_vector}" , "Penalty")
    self.plotMetric(history, "acc_train"          ,"Training Accuracy"                                                           , "Acc.")

    self.plotMetric(history, "xentropy_loss_val","Validation X-Entropy Loss"                                                     , "Loss"   )
    self.plotMetric(history, "total_loss_val"   ,"Validation Total Loss"                                                         , "Loss"   )
    self.plotMetric(history, "penalty_scalar_val"      ,r"Validation Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_scalar}" , "Penalty")
    self.plotMetric(history, "penalty_vector_val"      ,r"Validation Distance Correlation with $\lambda = $" + f"{self.lambda_penalty_vector}" , "Penalty")
    self.plotMetric(history, "acc_val"          ,"Validation Accuracy"                                                           , "Acc.")



    #self.plotAccuracy(history)
    #self.plotScoreTauOnly(model, xx_train, y_train, h_train,  xx_val, y_val, h_val, class_label )
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 0,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 1,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 2,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 3,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 4,class_label)
    #self.plotScore(model, xx_train, y_train, xx_val, y_val, 5,class_label)
    ##self.plotScoreOneVsAll(model, xx_train, y_train, xx_val, y_val, 1,class_label)

    #pdb.set_trace()
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


    return history, model



if __name__ == '__main__':


  ROOT.gROOT.SetBatch(True)
  ##limiting cores/CPU/GPU
  #num_cores = 8 #one operation at a time and only one thread per operation  
  #num_CPU = 8 
  #num_GPU = 0 
  #config = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=num_cores, inter_op_parallelism_threads=num_cores, allow_soft_placement = False, device_count = {'CPU' : num_CPU, 'GPU' : num_GPU} )
  #session = tf.compat.v1.Session(config=config) 
  #tf.compat.v1.keras.backend.set_session(session)         

  tf.random.set_seed(1000)
  np.random.seed(1000)
  
  features = kin_var 
  epochs = 1000
  #batch_size = 128 #128 here
  batch_size = 4096 #128 here
  learning_rate = 0.0005
  lambda_penalty_scalar = 10.0
  lambda_penalty_vector = 100.0
  scaler_type = 'robust'
  es_patience = 4000
  do_reduce_lr = False
  dirname = 'test'
  baseline_selection = baseline_selection 
  nfolds = 10 
  frac_sb = 0.0
  frac_sf = 1.0

  trainer = Trainer(
      features           = features, 
      epochs             = epochs,
      batch_size         = batch_size,
      learning_rate      = learning_rate,
      lambda_penalty_scalar     = lambda_penalty_scalar,
      lambda_penalty_vector     = lambda_penalty_vector,
      scaler_type        = scaler_type,
      es_patience        = es_patience,
      do_reduce_lr       = do_reduce_lr,
      dirname            = dirname,
      baseline_selection = baseline_selection,
      nfolds             = nfolds,
      frac_sb            = frac_sb,
      frac_sf            = frac_sf
      )

  history, model = trainer.process()




