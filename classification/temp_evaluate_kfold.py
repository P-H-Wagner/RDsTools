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

  def convertRootToDF(self, sample, training_info, treename, selection, weights, var, channel):
  

    #import pdb; pdb.set_trace()
    #add weights to extra_columns if any 
    #features are the same for all folds, so we can pick fold 0
    branches = dc(training_info[0].features) + [e for e in extra_vars if e not in dc(training_info[0].features) ]     
 
    if weights:
      branches += weights

    #if channel != "data":
    #  branches += [
    #   "gen_sig",
    #   "gen_match_success",
    #   "gen_mu_pt",
    #   "gen_tau_pt",
    #   "gen_ds_pt",
    #   "gen_dsStar_pt",
    #   "gen_bs_pt",
    #
    #   "gen_mu_eta",
    #   "gen_tau_eta",
    #   "gen_ds_eta",
    #   "gen_dsStar_eta",
    #   "gen_bs_eta",
    #
    #   "gen_mu_phi",
    #   "gen_tau_phi",
    #   "gen_ds_phi",
    #   "gen_dsStar_phi",
    #   "gen_bs_phi",
    #
    #   "gen_mu_pdgid",
    #   "gen_tau_pdgid",
    #   "gen_ds_pdgid",
    #   "gen_dsStar_pdgid",
    #   "gen_bs_pdgid",
    #
    #   "gen_ds_charge",
    #   "gen_bs_charge",

    #   "gen_mu_m",
    #   "gen_tau_m",
    #   "gen_ds_m",
    #   "gen_dsStar_m",
    #   "gen_bs_m",

    #   ]
      
    #evaluation on MC data
    #we can evaluate on the whole MC sample, since we checked that overtraining is not a problem :)


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

    #f = ROOT.TFile(sample)
    #tree_obj = f.Get(treename)
    #arr = tree2array(tree_obj,selection = selection, branches = branches)
    #df = pd.DataFrame(arr)
    #df = read_root(sample, treename, where=baseline_selection, warn_missing_tree=True, columns=training_info.features+extra_columns)

 
    pd.options.mode.chained_assignment = None   

    """
    if treename == 'tree_reco' and k == 'hb':
      print("save signal -1")
      #add signals to extra_columns
      
      df['tree_gen_sig'] = -1 #Hb

    
    if treename == 'tree_data' and k == 'left':

      print("save left sideband")
      #left sideband is assigned to a new sinal number, namely 5        

      df['tree_gen_sig'] = 5 #left

    if treename == 'tree_data' and k == 'right':

      print("save right sideband")
      #right sideband is assigned to a new sinal number, namely 6 

      df['tree_gen_sig'] = 6 #right
    """
    return df


  def createDataframe(self, training_info, files, selection, weights, treename, var, channel):
    '''
      Function that returns a dataframe out of a list of samples
    '''

    df = pd.concat([self.convertRootToDF(ifile, training_info, treename, selection, weights, var, channel) for ifile in files], sort=False)
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


  def createFileWithAnalysisTree(self, training_path, files, selection, weights, label, treename ,var, channel, constrained, nchunks):
    '''
      Create tree that contains q2_new,cospiK and the score and store it in a root file
      This tree is going to be used for the rest of the analysis (see fitter)
      Note that both the baseline selection and category definition have to be applied
    '''

    # get the training information
    training_info = self.getTrainingInfo(training_path)

    # create dataframe
    print( '\n ========> creating the dataframe')
    
    df = self.createDataframe(training_info=training_info, files=files, selection=selection, weights=weights, treename=treename,var = var, channel = channel)

    # get the score for each fold!
    x_folds = {}
    score   = {}

    for n in range(self.nfolds):
      x_folds[n] = df[ df["event"] % self.nfolds == n ]

      print( '\n ========> predicting the score')
      score[n]   = self.predictScore(training_info=training_info[n], df=x_folds[n]) #use head function for debugging on few events
      print(score)     
      import pdb
      pdb.set_trace()
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

    """
    # get other quantities to fill the branches with
    if channel != 'data':
     var.append("gen_sig")
     var.append("gen_match_success")

    var = dc(training_info.features) + var

    for name in var:
      globals()[name] = df[name]

    weight_val = {}
    if weights != None:
      for weight in weights:
        weight_val[weight] = df[weight]

    """


    #put all into one file
    if constrained: 
      name = "_cons"
    else:            
      name = "_unc"

    # create file
    #root_filename = f"/scratch/pahwagne/score_trees/{channel}/{self.channel}_{self.date_time}{name}.root"

    #concat all folds
    df_tot = pd.concat([ x_folds[n] for n in range(self.nfolds) ]) 

    #split df into n chunks for better post-processing 
    chunks = np.array_split(df_tot, nchunks)

    print("saving into directory", channel)

    for i, chunk in enumerate(chunks):
      root_filename = f"/scratch/pahwagne/score_trees/{channel}/{self.channel}_{self.date_time}_flatChunk_{i}.root"

      with uproot.recreate(root_filename) as f:
        f["tree"] = chunk 


    #else: 
    #  name = "_unc"
    #  #one file for every fold
    #  for n in range(self.nfolds):
   
    #    root_filename = f"/scratch/pahwagne/score_trees/{channel}/{self.channel}_{self.date_time}_fold_{n}.root"

    #    with uproot.recreate(root_filename) as f:
  
    #      f["tree"] = pd.concat([ x_folds[n] for n in range(self.nfolds) ])

    #"""
    #out_file = ROOT.TFile(root_filename, "RECREATE")

    ## create tree
    #print(' ========> creating tree')
    #tree = ROOT.TTree(f'tree', f'tree')

    ## initialise branches

    #for name in var:
    #  globals()["the_" + name] = array('d',[0])
    #the_score0 = array('d', [0])
    #the_score1 = array('d', [0])
    #the_score2 = array('d', [0])
    #the_class  = array('d', [0])
    ##the_score3 = array('d', [0])
    ##the_score4 = array('d', [0])
    ##the_score5 = array('d', [0])
    #
    #the_weight = {}
    #if weights != None:
    #  for weight in weights:
    #    the_weight[weight] = array('d', [0])

    #for name in var:
    #  tree.Branch(name, globals()["the_" + name], name + '/D')
    #tree.Branch('score0', the_score0, 'score0/D')
    #tree.Branch('score1', the_score1, 'score1/D')
    #tree.Branch('score2', the_score2, 'score2/D')
    #tree.Branch('class',  the_class,  'class/D')
    ##tree.Branch('score3', the_score3, 'score3/D')
    ##tree.Branch('score4', the_score4, 'score4/D')
    ##tree.Branch('score5', the_score5, 'score5/D')

    #if weights != None:
    #  for weight in weights:
    #    tree.Branch(weight, the_weight[weight], f'{weight}/D')

    ## fill the tree
    #for entry in range(len(score)):
   
    #  idx = np.argmax(score[entry,:])
 
    #  for name in var:
    #    dummy = globals()["the_" + name]
    #    dummy2 = globals()[name]
    #    dummy[0] = dummy2[entry]

    #  the_score0[0] = score0[entry]#[0]
    #  the_score1[0] = score1[entry]#[0]
    #  the_score2[0] = score2[entry]#[0]
    #  the_class[0]  = idx#[0]
    #  #the_score3[0] = score3[entry][0]
    #  #the_score4[0] = score4[entry][0]
    #  #the_score5[0] = score5[entry][0]

    #  if weights != None:
    #    for weight in weights:
    #      the_weight[weight][0] = weight_val[weight][entry]

    #  tree.Fill()
    #tree.Write()
    #out_file.Close()
    #"""
    print(f' ========> {root_filename} created')

if __name__ == '__main__':

   # parsing shell variable
   constrained  = os.getenv("constrained") 
   channel      = os.getenv("channel") 
   modelpath    = os.getenv("modelpath") 
   nchunks      = os.getenv("nchunks")
   prod         = os.getenv("prod")
   print("parsing command line... evaluate on: ", channel, nchunks)   

   
   nchunks = int(nchunks)
   if constrained == 'True': constrained = True
   else: constrained = False

 
   files = {}
    
   if constrained:
     if prod == "24":
     
       baseline_selection = base_wout_tv_24
      
       #bdt is evaluated on the skimmed datasets :) 
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_24}/"
       files["data"] = [base + f for f in os.listdir(base)]
      
       #hammer is evaluated on the skimmed datasets :) 
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/signal_default_29_04_2025_13_51_27/" 
       files["sig"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_24[0]}/"
       files["hb"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bs_cons_24[0]}/"
       files["bs"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{b0_cons_24[0]}/"
       files["b0"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bplus_cons_24[0]}/"
       files["bplus"] = [base + f for f in os.listdir(base)]
     
     elif prod == "25":
     
       baseline_selection = base_wout_tv_25
       # for now we only consider mu7
       baseline_selection += " && (mu7_ip4) && (mu_is_global) && (ds_vtx_cosine_xy_pv > 0.8)"
     
       #bdt is evaluated on the skimmed datasets :) 
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_25}/"
       files["data"] = [base + f for f in os.listdir(base)]
       #hammer is evaluated on the skimmed datasets :) 
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/signal_default_16_05_2025_10_36_51/" 
       files["sig"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_25[0]}/"
       files["hb"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bs_cons_25[0]}/"
       files["bs"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{b0_cons_25[0]}/"
       files["b0"] = [base + f for f in os.listdir(base)]
       
       base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bplus_cons_25[0]}/"
       files["bplus"] = [base + f for f in os.listdir(base)]

    
     #normally 
     ##this is only constrained data! (cut on fv only)
     #files["data"]  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_cons ]

     ##not hammered
     ##files["sig"]  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_cons  ]
     #direc          = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{sig_cons[0]}/"

     ##hammered
     #direc          = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/signal_default_13_03_2025_08_42_43/"

     #files["sig"]   = [os.path.join(direc, f) for f in os.listdir(direc)]

     #files["hb"]    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_cons   ]
     #files["b0"]    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_cons   ]
     #files["bs"]    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_cons   ]
     #files["bplus"] = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_cons]

   else:
     #this is only unconstrained data! (cut on fv only)
     files["data"]  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_unc ]
     files["sig"]   = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_unc  ]
     files["hb"]    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_unc   ]
     files["b0"]    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_unc   ]
     files["bs"]    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_unc   ]
     files["bplus"] = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_unc]


   training_path = f'/work/pahwagne/RDsTools/classification/outputs/test_{modelpath}' #mu and tau in same class  (adam) uncon.
   #date_time   = modelpath[5:] 
   #print("date_time:", date_time)

   print("training path: ", training_path)

   print(f"=====> Using model {modelpath}.") 

   weights=None

   #label under which we would like to store the new trees
   labels = {"sig": "sig_pastNN", "data": "data_pastNN", "hb": "hb_pastNN", "bs":"bs_pastNN", "b0": "b0_pastNN", "bplus":"bplus_pastNN" }

   #treenames
   treename = "tree"

   force_overwrite=False

   # sig variables which we would like to save in the new trees (just take one file)
   test_file_sig   = ROOT.TFile(files["sig"][0])
   test_tree_sig   = test_file_sig.Get("tree")
   var_sig         = [branch.GetName() for branch in test_tree_sig.GetListOfBranches()]

   # hb variables which we would like to save in the new trees (just take one file)
   test_file_hb   = ROOT.TFile(files["hb"][0])
   test_tree_hb   = test_file_hb.Get("tree")
   var_hb         = [branch.GetName() for branch in test_tree_hb.GetListOfBranches()]

   # data variables which we would like to save in the new trees (just take one file)
   test_file_data = ROOT.TFile(files["data"][0])
   test_tree_data = test_file_data.Get("tree")
   var_data       = [branch.GetName() for branch in test_tree_data.GetListOfBranches()]

   #read features from settings.txt file of the corresponding model

   features = []
   with open( training_path + "/settings.txt", "r") as f:
       for line in f:
           if line.strip().startswith("Features:"):
               match = re.search(r"Features:\s*(\[.*\])", line)
               if match:
                   features = ast.literal_eval(match.group(1))
                   break
   
   print(features)  


   #append sf weights for data
   if (channel == "data" and constrained) :  # throw away all the refitted quantities (not needed for final fit)
     extra_vars = [e for e in var_data if e not in features ]

   elif (channel == "sig" and constrained) :
     extra_vars = [e for e in var_sig  if e not in features ]

   else:
     extra_vars = [e for e in var_hb   if e not in features ]

   print(extra_vars)
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
   if channel != "data": 
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

   #HAMMER
   #if ((channel == "sig") and constrained) :

   #  print("Adding hammer variables") 
   #  extra_vars.append("central_w")
   #  for i in range(1,11):
   #    extra_vars.append(f"e{i}_up")
   #    extra_vars.append(f"e{i}_down")
       

   if prod == "24": 
     cut = base_wout_tv_24

   else: 
     cut = base_wout_tv_25 + " && (mu7_ip4) && (mu_is_global) && (ds_vtx_cosine_xy_pv > 0.8)"

   #selection for mc, exlude hb from the inclusive sample and include it explicitly with the hb only sample
   selections = {"sig":cut, "data":cut, "bs":cut ,"b0":cut , "bplus": cut , "hb": cut}

   print(f"========> Evaluating on {channel} ...")

   #create df on which we want to evaluate
   tools = createDf(modelpath, channel) 

   #create tree
   FileWithScore = tools.createFileWithAnalysisTree(training_path, files[channel], selections[channel], weights, labels[channel], treename, extra_vars,channel, constrained = constrained, nchunks = nchunks)

   #/work/pahwagne/ma22/programs/Kinematical_variable_programs/outputs/test_18Jan2023_09h36m59s is the last test we did of trainer.py with class_weight == class_weights






