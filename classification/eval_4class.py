import os
import sys
from os import path

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
  def __init__(self, training_path, model_label):
    
    self.model_label = model_label
    self.indir = training_path
    self.model = self.loadModel()
    self.qt = self.loadScaler()
    self.features = self.loadFeatures()

    if not path.exists(self.indir):
      raise RuntimeError(f'Model with label "{self.model_label}" was not to be found in "{self.indir}".')


  def loadModel(self):
    print( '\n ========> getting the model')
    
    model_filename = f'{self.indir}/saved-model-{self.model_label}.h5'
    model = tf.keras.models.load_model(model_filename)
    return model


  def loadScaler(self):
    print( '\n ========> getting the scaler')
    
    scaler_filename = f'{self.indir}/input_tranformation_weighted.pck'
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

  def convertRootToDF(self, sample, training_info, treename, baseline_selection, weights, channel):
   
    #add weights to extra_columns if any 
    if weights:
      extra_columns = weights

    if channel != 'data':
         
      #evaluation on MC data
      #we can evaluate on the whole MC sample, since we checked that overtraining is not a problem :)

      print("evaluate on MC ...")

      f = ROOT.TFile(sample)
      tree_obj = f.Get(treename)
      arr = tree2array(tree_obj,selection = baseline_selection, branches = training_info.features+extra_columns)
      df = pd.DataFrame(arr)
      #df = read_root(sample, treename, where=baseline_selection, warn_missing_tree=True, columns=training_info.features+extra_columns)
      
      pd.options.mode.chained_assignment = None   

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



    else:
     
      #evaluation on real data
      #we evaluate on the whole sample!
      print("evaluate on data ...")

      f = ROOT.TFile(sample)
      tree_obj = f.Get(treename)
      arr = tree2array(tree_obj,selection = baseline_selection, branches = training_info.features+extra_columns)
      df = pd.DataFrame(arr)

      #df = read_root(sample, treename, where=baseline_selection, warn_missing_tree=True, columns=training_info.features+extra_columns)
      pd.options.mode.chained_assignment = None   

    return df


  def createDataframe(self, training_info, files, minimal_selections, weights, treenames,key):
    '''
      Function that returns a dataframe out of a list of samples
    '''
    df = pd.concat([self.convertRootToDF(idt, training_info, name, selec, weights,k) for idt,name,selec,k in zip(files,treenames,minimal_selections,key)], sort=False)

    # remove inf and nan
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)

    # re-index
    df = df.reset_index(drop=True)

    return df


  def getTrainingInfo(self, training_path, model_label):
    '''
      Get training information
    '''
    training_info = TrainingInfo(training_path=training_path, model_label=model_label)

    return training_info


  def predictScore(self, training_info, df):
    '''
      Return score with scaled input features
    '''
    x = pd.DataFrame(df, columns=training_info.features)

    # apply the scaler
    xx = training_info.qt.transform(x[training_info.features])
    print(type(xx))
    # predict
    score = training_info.model.predict(xx)

    return score


  def createFileWithAnalysisTree(self, training_path, model_label, files, minimal_selections, weights, label, treenames,var, channel):
    '''
      Create tree that contains q2_new,cospiK and the score and store it in a root file
      This tree is going to be used for the rest of the analysis (see fitter)
      Note that both the baseline selection and category definition have to be applied
    '''
    # get the training information
    training_info = self.getTrainingInfo(training_path, model_label)

    # create dataframe
    print( '\n ========> creating the dataframe')
    df = self.createDataframe(training_info=training_info, files=files, minimal_selections=minimal_selections, weights=weights, treenames=treenames,channel = channel)

    # get the score
    print( '\n ========> predicting the score')
    score = self.predictScore(training_info=training_info, df=df )
    print(score)     
    score0 = score[:,0] 
    score1 = score[:,1] 
    score2 = score[:,2] 
    score3 = score[:,3] 

    # get other quantities to fill the branches with
 

    for name in var:
      globals()[name] = df[name]

    """
    q2_new = df['q2_new']
    cospiK = df['cospiK']
    m2_miss_new = df['m2_miss_new']
    costheta1 = df['costheta1']
    cosplanes = df['cosplanes']
    mu_pt = df['mu_pt']
    dr_m_ds = df['dr_m_ds']
    dr_m_pi = df['dr_m_pi']
    mDs = df['mDs']
    e_star_mu = df['e_star_mu']
    pt_miss = df['pt_miss']
    mu_eta = df['mu_eta']
    """

    if key[0] != 'data':
      tree_gen_sig = df['tree_gen_sig']
     

    weight_val = {}
    if weights != None:
      for weight in weights:
        weight_val[weight] = df[weight]

    # create file
    root_filename = f'./score_trees/{label}_' + datetime.now().strftime('%d%b%Y_%Hh%Mm%Ss') + '.root'
    out_file = ROOT.TFile(root_filename, 'RECREATE')

    # create tree
    print(' ========> creating tree')
    tree = ROOT.TTree(f'{label}_tree', f'{label}_tree')

    # initialise branches
    for name in var:
      globals()["the_" + name] = array('d',[0])
    the_score0 = array('d', [0])
    the_score1 = array('d', [0])
    the_score2 = array('d', [0])
    the_score3 = array('d', [0])
    the_score4 = array('d', [0])
    the_score5 = array('d', [0])


    """
    the_score = array('d', [0])
    the_q2_new = array('d', [0])
    the_cospiK = array('d', [0])
    the_m2_miss_new =array('d', [0])
    the_costheta1 = array('d', [0])
    the_cosplanes = array('d', [0])
    the_mu_pt = array('d', [0])
    the_dr_m_ds =array('d', [0])
    the_dr_m_pi = array('d', [0])
    the_mDs = array('d', [0])
    the_e_star_mu = array('d', [0])
    the_pt_miss= array('d', [0])
    the_mu_eta = array('d', [0])
    """
    if key[0] != 'data':
      the_tree_gen_sig = array('d', [0])

    
    the_weight = {}
    if weights != None:
      for weight in weights:
        the_weight[weight] = array('d', [0])

    for name in var:
      tree.Branch(name, globals()["the_" + name], name + '/D')
    tree.Branch('score0', the_score0, 'score0/D')
    tree.Branch('score1', the_score1, 'score1/D')
    tree.Branch('score2', the_score2, 'score2/D')
    tree.Branch('score3', the_score3, 'score3/D')
    tree.Branch('score4', the_score4, 'score4/D')
    tree.Branch('score5', the_score5, 'score5/D')

    """
    tree.Branch('score', the_score, 'score/D')
    tree.Branch('q2_new', the_q2_new, 'q2_new/D')
    tree.Branch('cospiK', the_cospiK, 'cospiK/D')
    tree.Branch('m2_miss_new',the_m2_miss_new,'m2_miss_new/D')
    tree.Branch('cosplanes',the_cosplanes,'cosplanes/D')
    tree.Branch('costheta1',the_costheta1,'costheta1/D')
    tree.Branch('mu_pt',the_mu_pt,'mu_pt/D')
    tree.Branch('dr_m_pi',the_dr_m_pi,'dr_m_pi/D')
    tree.Branch('dr_m_ds',the_dr_m_ds,'dr_m_ds/D')
    tree.Branch('mDs',the_mDs,'mDs/D')
    tree.Branch('e_star_mu',the_e_star_mu,'e_star_mu/D')
    tree.Branch('pt_miss',the_pt_miss,'pt_miss/D')
    tree.Branch('mu_eta',the_mu_eta,'mu_eta/D')
    """

    if key[0] != 'data':
      tree.Branch('tree_gen_sig',the_tree_gen_sig,'tree_gen_sig/D')
 
    if weights != None:
      for weight in weights:
        tree.Branch(weight, the_weight[weight], f'{weight}/D')

    # fill the tree
    if key[0] != 'data': 
      for entry in range(len(score)):
      
        for name in var:
          dummy = globals()["the_" + name]
          dummy2 = globals()[name]
          dummy[0] = dummy2[entry]

        the_tree_gen_sig[0] = tree_gen_sig[entry]
	the_score0[0] = score0[entry][0]
        the_score1[0] = score1[entry][0]
        the_score2[0] = score2[entry][0]
        the_score3[0] = score3[entry][0]
        the_score4[0] = score4[entry][0]
        the_score5[0] = score5[entry][0]

        """ 
        the_score[0] = score[entry][0]
        the_q2_new[0] = q2_new[entry]
        the_cospiK[0] = cospiK[entry]
        the_m2_miss_new[0] = m2_miss_new[entry]
        the_cosplanes[0] = cosplanes[entry]
        the_costheta1[0] = costheta1[entry]
        the_mu_pt[0] = mu_pt[entry]
        the_dr_m_pi[0] = dr_m_pi[entry]
        the_dr_m_ds[0] = dr_m_ds[entry]
        the_tree_gen_sig[0] = tree_gen_sig[entry]
        the_mDs[0] = mDs[entry]
        the_e_star_mu[0] = e_star_mu[entry]
        the_pt_miss[0] = pt_miss[entry]
        the_mu_eta[0] = mu_eta[entry]
        """
        if weights != None:
          for weight in weights:
            the_weight[weight][0] = weight_val[weight][entry]

        tree.Fill()
      tree.Write()
      out_file.Close()

    else:
      for entry in range(len(score)):

        for name in var:
          dummy = globals()["the_" + name]
          dummy2 = globals()[name]
          dummy[0] = dummy2[entry]

	the_score0[0] = score0[entry][0]
        the_score1[0] = score1[entry][0]
        the_score2[0] = score2[entry][0]
        the_score3[0] = score3[entry][0]
        the_score4[0] = score4[entry][0]
        the_score5[0] = score5[entry][0]

        """
        the_score[0] = score[entry][0]
        the_q2_new[0] = q2_new[entry]
        the_cospiK[0] = cospiK[entry]
        the_m2_miss_new[0] = m2_miss_new[entry]
        the_cosplanes[0] = cosplanes[entry]
        the_costheta1[0] = costheta1[entry]
        the_mu_pt[0] = mu_pt[entry]
        the_dr_m_pi[0] = dr_m_pi[entry]
        the_dr_m_ds[0] = dr_m_ds[entry]
        the_mDs[0] = mDs[entry]
        the_e_star_mu[0] = e_star_mu[entry]
        the_pt_miss[0] = pt_miss[entry]
        the_mu_eta[0] = mu_eta[entry]
        """
        if weights != None:
          for weight in weights:
            the_weight[weight][0] = weight_val[weight][entry]

        tree.Fill()
      tree.Write()
      out_file.Close()


    print(f' ========> {root_filename} created')

if __name__ == '__main__':


  
   #this is only unconstrained data! (cut on fv and tv)
   files_data  = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in data_unc ]
   files_sig   = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in sig_unc  ]
   files_hb    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in hb_unc   ]
   files_b0    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in b0_unc   ]
   files_bs    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in bs_unc   ]
   files_bplus = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_{f}.root" for f in bplus_unc]
   #this is only unconstrained data! (cut on fv only)
   file_data  = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_unc ]
   file_sig   = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_unc  ]
   file_hb    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_unc   ]
   file_b0    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_unc   ]
   file_bs    = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_unc   ]
   file_bplus = ["/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_unc]

 
   #specify model
   training_path = '/work/pahwagne/RDsTools/classification/outputs/test_30Jul2024_14h16m55s'  
   model_label = '0020_val_loss_0.3508_val_acc_0.8379' 

   weights=None

   #label under which we would like to store the new trees
   label_mc ='pastNN'
   label_data = 'pastNN_big_data_multi'

   #treenames
   tree_name = 'tree'

   force_overwrite=False

   #apply save selection as before training #TODO automize this 
   selection = ma_cut

   # mc variables which we would like to save in the new trees (just take one file)
   test_file_mc   = ROOT.TFile(files_sig[0])
   test_tree_mc   = test_file_mc.Get("tree")
   var_mc         = [branch.GetName() for branch in test_tree_mc.GetListOfBranches()]

   # data variables which we would like to save in the new trees (just take one file)
   test_file_data = ROOT.TFile(files_data[0])
   test_tree_data = test_file_data.Get("tree")
   var_data       = [branch.GetName() for branch in test_tree_mc.GetListOfBranches()]

   # get consants for sb method
   with open('massfit.pickle', 'rb') as handle:
     massfit = pickle.load(handle)

   sigma = massfit["sigma"]

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


   #selection for mc, exlude hb from the inclusive sample and include it explicitly with the hb only sample
   minimal_selections_mc = [minimal_selection + "& tree_gen_sig != -1", minimal_selection + SBL, minimal_selection + SBR,minimal_selection]
   #selection for data
   minimal_selections_data = [minimal_selection]
  
   key_mc = ["mc","left","right","hb"]
   key_data = ["data"]

   print("DO YOU WANT TO EVALUATE ON DATA ('data') OR ON MC ('mc')?")
   evaluation_type = 'data'
   print(f"========> you selected {evaluation_type}")

   #create df on which we want to evaluate
   tools = createDf() 

   if evaluation_type == 'mc':
      label = label_mc
      files = files_mc
      treenames = treenames_mc
      minimal_selections = minimal_selections_mc   
      key = key_mc

   else:
      label = label_data
      files = files_data
      treenames = treenames_data
      minimal_selections = minimal_selections_data   
      key = key_data

   #create tree
   FileWithScore = tools.createFileWithAnalysisTree(files, training_path, model_label, weights, label, minimal_selections, force_overwrite, tree_name, var, channel)

   #/work/pahwagne/ma22/programs/Kinematical_variable_programs/outputs/test_18Jan2023_09h36m59s is the last test we did of trainer.py with class_weight == class_weights







