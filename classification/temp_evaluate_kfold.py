import os
import sys
from os import path
import uproot
import ROOT

#from root_pandas import read_root
from root_numpy import tree2array
import pickle
import numpy as np
import pandas as pd
from array import array

import sklearn as sk
from sklearn.model_selection import train_test_split

import tensorflow as tf
from datetime import datetime

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from helper import *
import argparse
from copy import deepcopy as dc
import re
import glob
import ast
import re

#########################################################

#####evaluates the BINARY classifier on data or mc#######

#########################################################



def limitCPU(n):
  '''
  limit CPU usage, not needed when assigned to batch 
  '''
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

class TrainingInfo(object):
  '''
  Get all the training info
  '''
  def __init__(self, training_path, model_label, fold):
    
    self.model_label = model_label
    self.indir = training_path
    self.fold = fold
    self.model = self.loadModel()
    self.qt = self.loadScaler()
    self.features = self.loadFeatures()


    if not path.exists(self.indir):
      raise RuntimeError(f'Model with label "{self.model_label}" was not to be found in "{self.indir}".')


  def loadModel(self):
    print( '\n ========> getting the model')
    
    model_filename = f'{self.indir}/fold_{self.fold}_saved-model-{self.model_label}.h5'
    model = tf.keras.models.load_model(model_filename)
    return model


  def loadScaler(self):
    print( '\n ========> getting the scaler')
    
    scaler_filename = f'{self.indir}/input_tranformation_weighted_fold_{self.fold}.pck'
    qt = pickle.load(open(scaler_filename, 'rb'))
    return qt


  def loadFeatures(self):
    print( '\n ========> getting the features')
    
    features_filename = f'{self.indir}/input_features.pck'
    features = pickle.load(open(features_filename, 'rb'))
    
    return features


class createDf(object):
  '''
  create df for evaluation on either mc or data
  '''

  def __init__(self,date_time, channel):
    self.date_time = date_time
    self.channel   = channel 

  def convertRootToDF(self, sample, training_info, treename, selection, var, channel):
  

    #import pdb; pdb.set_trace()
    #features are the same for all folds, so we can pick fold 0
    branches = dc(training_info[0].features) + [e for e in extra_vars if e not in dc(training_info[0].features) ]     
 
    #use uproot!!
    with uproot.open(sample) as f:
      tree = f[treename]
      #print(selection)
      selection = selection.replace("&&", "&")
      #print(selection)
      selection = selection.replace("||", "|")
      #print(selection)
      #import pdb; pdb.set_trace()
      df   = tree.arrays(branches, library = "pd", entry_start=None, entry_stop=None, cut=selection)

    pd.options.mode.chained_assignment = None   

    return df


  def createDataframe(self, training_info, files, selection,treename, var, channel):
    '''
      Function that returns a dataframe out of a list of samples
    '''

    df = pd.concat([self.convertRootToDF(ifile, training_info, treename, selection, var, channel) for ifile in files], sort=False)
    print("before removing nans i have:", len(df))
    # remove inf and nan

    df.replace([np.inf, -np.inf], np.nan, inplace=True)

    #collect all columns except the gen ones
    #non_nan = [col for col in df.keys() if ("gen_" not in col)]

    #select only rows in which we dont have nans in the "non_nan" columns. 
    #gen can contain nans, as not for every event we have tau or dsstar!
    #df.dropna(subset = non_nan, inplace = True) #inplace allows to directly overwrite the df and avoid df = df.dropna(...)

    print("after removing nans i have:", len(df))
    # re-index
    df = df.reset_index(drop=True)

    return df


  def getTrainingInfo(self, training_path):
    '''
      Get training information
    '''

    training_info = {}

    #get complete model names
    all_models = glob.glob(training_path + "/fold*.h5")

    #find number of folds
    pattern = r'fold_(\d+)'
    strings = [ re.search(pattern,model).group(1) for model in all_models     ] 
    ints    = [ int(re.search(pattern,model).group(1)) for model in all_models ]
    argmax  = ints.index(max(ints))
    nfolds  = strings[argmax]
    self.nfolds = int(nfolds) + 1 #f.e.: if the strings go from 0 to 4 -> means we have 5 folds

    #for each fold, find best model (the newest)
    models = {}

    for n in range(self.nfolds):
      all_models = glob.glob(training_path + f"/fold_{n}*.h5")
      pattern = rf'fold_{n}_saved-model-(\d+)'
      #take float here to avoid leading 0 problem
      strings   = [ re.search(pattern,model).group(1) for model in all_models]
      floats    = [ float(re.search(pattern,model).group(1)) for model in all_models]
      argmax    = floats.index(max(floats))
      nmodel    = strings[argmax]

      full_path = glob.glob(training_path + f'/fold_{n}_saved-model-{nmodel}*.h5')[0]
      pattern   = r'saved-model-(.+?)\.h5' 

      print(f"For fold {n} we pick model {re.search(pattern, full_path).group(1)}")

      training_info[n] = TrainingInfo(training_path = training_path, model_label = re.search(pattern, full_path).group(1), fold = n)
   

    return training_info


  def predictScore(self, training_info, df):
    '''
      Return score with scaled input features
    '''
    x = pd.DataFrame(df, columns=training_info.features)
    print(training_info.features)
    # apply the scaler
    xx = training_info.qt.transform(x[training_info.features])
    print(xx)

    # predict 
    score = training_info.model.predict(xx)

    return score


  def createFileWithAnalysisTree(self, training_path, files, selection, label, treename ,var, channel):
    '''
      Create tree that contains q2_new,cospiK and the score and store it in a root file
      This tree is going to be used for the rest of the analysis (see fitter)
      Note that both the baseline selection and category definition have to be applied
    '''

    # get the training information
    training_info = self.getTrainingInfo(training_path)

    # create dataframe
    print( '\n ========> creating the dataframe')
    
    df = self.createDataframe(training_info=training_info, files=files, selection=selection, treename=treename,var = var, channel = channel)

    # get the score for each fold!
    x_folds = {}
    score   = {}

    for n in range(self.nfolds):
      x_folds[n] = df[ df["event"] % self.nfolds == n ]

      print( '\n ========> predicting the score')
      score[n]   = self.predictScore(training_info=training_info[n], df=x_folds[n]) #use head function for debugging on few events
      print(score)     
      all_scores = []
    
      score0 = score[n][:,0] 
      score1 = score[n][:,1] 
      score2 = score[n][:,2] 
      score3 = score[n][:,3] 
      score4 = score[n][:,4] 
      score5 = score[n][:,5] 

      all_scores.append(score0)    
      all_scores.append(score1)    
      all_scores.append(score2)    
      all_scores.append(score3)    
      all_scores.append(score4)    
      all_scores.append(score5)    

      all_scores = np.array(all_scores)
      max_score = [ np.argmax( all_scores[:,i] ) for i in range(len(all_scores[0]))] 
      print("length of max score",len(max_score))
      print("length of all_scores",len(all_scores[:,0] ))

      x_folds[n]["score0"] = score0
      x_folds[n]["score1"] = score1
      x_folds[n]["score2"] = score2
      x_folds[n]["score3"] = score3
      x_folds[n]["score4"] = score4
      x_folds[n]["score5"] = score5
      x_folds[n]["class"] = max_score


    #concat all folds
    df_tot = pd.concat([ x_folds[n] for n in range(self.nfolds) ]) 

    root_filename = f"/scratch/pahwagne/score_trees/HOOK_CHANNEL/HOOK_CHANNEL_HOOK_MODELPATH_flatChunk_HOOK_CHUNK.root"

    with uproot.recreate(root_filename) as f:
        f["tree"] = df_tot 

    print(f' ========> {root_filename} created')

if __name__ == '__main__':

   # files
   files = [HOOK_FILES]

   # training path/model 
   training_path = f'/work/pahwagne/RDsTools/classification/outputs/test_HOOK_MODELPATH' 
   print("training path: ", training_path)

   print("=====> Using model HOOK_MODELPATH.") 

   #label under which we would like to store the new trees
   label = "HOOK_CHANNEL_pastNN"

   #treenames
   treename = "tree"

   force_overwrite=False

   features = []
   with open( training_path + "/settings.txt", "r") as f:
       for line in f:
           if line.strip().startswith("Features:"):
               match = re.search(r"Features:\s*(\[.*\])", line)
               if match:
                   features = ast.literal_eval(match.group(1))
                   break
   

   #append sf weights for data
   if   ("HOOK_CHANNEL" == "data" ) :  
     # data variables which we would like to save in the new trees (just take one file)
     test_file_data = ROOT.TFile(files[0])
     test_tree_data = test_file_data.Get("tree")
     var_data       = [branch.GetName() for branch in test_tree_data.GetListOfBranches()]
     extra_vars = [e for e in var_data if e not in features ]

   elif ("HOOK_CHANNEL" == "sig") :
     # sig variables which we would like to save in the new trees (just take one file)
     test_file_sig   = ROOT.TFile(files[0])
     test_tree_sig   = test_file_sig.Get("tree")
     var_sig         = [branch.GetName() for branch in test_tree_sig.GetListOfBranches()]
     extra_vars      = [e for e in var_sig  if e not in features ]

   else:
     # hb variables which we would like to save in the new trees (just take one file)
     test_file_hb   = ROOT.TFile(files[0])
     test_tree_hb   = test_file_hb.Get("tree")
     var_hb         = [branch.GetName() for branch in test_tree_hb.GetListOfBranches()]
     extra_vars = [e for e in var_hb   if e not in features ]

   print("Branches before clipping:", len(extra_vars))

   # clip further (these trees we only need for the final fit

   # remove the reco weighted columns (for now) 
   extra_vars = [e for e in extra_vars if "_weighted"   not in e ]

   # remove all the refitted quantities 
   extra_vars = [e for e in extra_vars if "fitted"   not in e ]

   # remove all the prescale info
   extra_vars = [e for e in extra_vars if "prescale" not in e ]

   # remove the track matching info
   extra_vars = [e for e in extra_vars if "matched"  not in e ]

   # remove all gen variables, except the signal info
   extra_vars = [e for e in extra_vars if "gen_"      not in e ]
   if "HOOK_CHANNEL" != "data": 
     extra_vars.append("gen_sig")
     extra_vars.append("gen_match_success")

   # remove all the PV info, except the custom one
   extra_vars = [e for e in extra_vars if "pv_general" not in e and "pv_dz" not in e]

   # remove all the isolation info, except the custom w.r.t PV
   extra_vars = [e for e in extra_vars if "iso_03_sv" not in e]
   extra_vars = [e for e in extra_vars if "iso_03_tv" not in e]
   extra_vars = [e for e in extra_vars if "iso_04_sv" not in e]
   extra_vars = [e for e in extra_vars if "iso_04_tv" not in e]

   # remove all the photon doubled variables 
   extra_vars = [e for e in extra_vars if "photon" not in e]
   extra_vars.append("photon_pt")
   extra_vars.append("photon_eta")

   #extra_vars.append("dsPhoton_m") #this is not removed as it has "P"
   extra_vars.append("bs_mass_corr_photon")
   extra_vars.append("ds_perp_photon")
   extra_vars.append("ds_mu_perp_photon")

   #remove all the btv isolation branches
   extra_vars = [e for e in extra_vars if "btv" not in e]

   #remove signed impact 
   extra_vars = [e for e in extra_vars if "signed_impact" not in e]

   #remove cosPlanesBs (what is this even!?)
   extra_vars = [e for e in extra_vars if "cosPlaneBs" not in e]

   #remove dxy 
   extra_vars = [e for e in extra_vars if "dxy" not in e]

   print("Branches after clipping:", len(extra_vars))

   print("========> Evaluating on HOOK_CHANNEL ...")

   #create df on which we want to evaluate
   tools = createDf("HOOK_MODELPATH", "HOOK_CHANNEL") 

   #create tree                                                    #list       #string
   FileWithScore = tools.createFileWithAnalysisTree(training_path, files, "HOOK_SELEC", label, treename, extra_vars, "HOOK_CHANNEL")







