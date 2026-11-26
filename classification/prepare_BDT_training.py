import pdb
import pickle
import uproot
import xgboost as xgb
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import sys
import os
import glob
import time
import ROOT
import json 
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, confusion_matrix
from datetime import datetime
import argparse
import yaml
import seaborn

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.EnableImplicitMT()

###########################
# Parse arguments         #
###########################

parser = argparse.ArgumentParser()
parser.add_argument('-dt', '--datetime'  , required=True)
parser.add_argument("--sel"              , required = True, help = "specify baseline selection (choose one in helper.py)")
parser.add_argument('-d' , '--debug'     , action='store_true' )
args = parser.parse_args()


if args.sel not in baselines.keys() :
  raise ValueError ("Error: Not a valid key for --sel, please use one defined in helper.py")

debug = False
if args.debug: debug = args.debug

scratch_dir = f"/scratch/pahwagne/{args.datetime}"

with open("/work/pahwagne/RDsTools/hammercpp/development_branch/weights/average_weights.yaml","r") as f:
  averages = yaml.safe_load(f)

#date time
now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")
path = f"./bdt_outputs/bdt_{dt}/"
os.makedirs(path)

#number of folds
nfolds = 10

###########################
# Define features         #
###########################

# Define features and target
features = [

# helicity
#'cosPhiDs_lhcb_alt',
'cosPiK1',
'cosMuW_lhcb_alt',
#'cosMuW_reco_1',
#'cosMuW_reco_2',

# vertex info
'fv_chi2', # dont take prob, they can hit the prec. limit when small!
'tv_chi2', # "
'sv_chi2', # "

# displacement
#'lxy_bs_sig', #contains nans
#'lxy_ds_sig', #contains nans

#'dxy_mu_sig_pv',
#'dxy_mu_sig_sv',

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
#'kk_m',
#'phiPi_m',
#'dsMu_m',

# delta R
#'kk_deltaR',
#'phiPi_deltaR',
#'dsMu_deltaR',

# q2 and co
#'q2_coll',
#'e_star_coll',

#'log(bs_boost_lhcb_alt)',
'bs_pt_lhcb_alt',
'e_star_lhcb_alt',
'm2_miss_lhcb_alt',
'pt_miss_lhcb_alt',

#'log(bs_boost_reco_1)',
'bs_pt_reco_1',
'e_star_reco_1',
'q2_reco_1',
'pt_miss_reco_1',
'm2_miss_reco_1',

#'log(bs_boost_reco_2)',
'bs_pt_reco_2',
'e_star_reco_2',
'q2_reco_2',
'pt_miss_reco_2',
'm2_miss_reco_2',

#'disc_negativity',
#'ds_vtx_cosine_xy_pv',
#'kappa',

#'signed_decay_ip3d_mu_ds_sv',
#'signed_decay_ip3d_mu_bs_sv',

#'rel_iso_03',
'rel_iso_03_pv',
#'rel_iso_03_sv',
#'rel_iso_03_ds_sv',
#'rel_iso_03_tv',

"e_gamma",

#"photon_pt",
#"dsPhoton_m",
#"dsPhotonMu_m",

"bs_mass_corr",
#"bs_mass_corr",
#"bs_mass_corr_photon)",
"ds_perp",
#"ds_perp_photon",
#"ds_mu_perp",
#"ds_mu_perp_photon"


]

#needed for filtering, selection, etc.
more_vars = ["lxy_ds","mu_id_medium","fv_prob","mu_is_global",
"ds_vtx_cosine_xyz_pv",
"mu7_ip4",
"mu8_ip3",
"mu8_ip5",
"mu8_ip6",
"mu9_ip4",
"mu9_ip5",
"mu9_ip6",
"mu12_ip6",
"k1_charge", "k2_charge", "pi_charge", "mu_charge"
]


###########################
# Define helper functions #
###########################
def getDf(path, branches=None, selection=None, debug=None):

    import ROOT
    import glob
    import pandas as pd

    print(f"Start loading df for path {path}")

    if isinstance(path, list):
        files = []
        for p in path:
            files += glob.glob(p)
    else:
        files = glob.glob(path)

    if debug:
        print("Debug mode: using only 3 files")
        files = files[:3]

    rdf = ROOT.RDataFrame("tree", files)

    if selection:
        print(selection)
        rdf = rdf.Filter(selection)  # MUST be C++ syntax

    if branches:
        rdf = rdf.AsNumpy(branches)

    df = pd.DataFrame(rdf)

    nmax = 1000000
    if len(df) > nmax: 
      df = df.sample(n=nmax, random_state=42, replace = False)

    return df

#def getDf(path, branches = None, selection = None, debug = None):
#
#   
#
#  print(f"Start loading df for path {path}")
#
#  #make it python compatible
#  selection = selection.replace("||"," or ")
#  selection = selection.replace("&&"," and ")
#
#  if isinstance(path,list):
#    files = []
#    for p in path:
#      files += glob.glob(p)
#  else:
#    files = glob.glob(path)
#
#
#  if debug:
#    print("Loading less files for debugging...")
#    files = files[:3]
#    print(files)
#
#  pd_list = [] 
#
#  total = int(len(files))
#
#  for i,name in enumerate(files):
#
#    tree = uproot.open(name)["tree"]
#
#    if selection:
#      #arr = tree.arrays(
#      #    expressions=branches,
#      #    cut=selection,
#      #    library="pd"
#      #)
#     
#      arr = tree.arrays(branches, library="pd")
#      arr = arr.query(selection, engine="python") 
# 
#    else:
#      arr = tree.arrays(
#          expressions=branches,
#          library="pd"
#      )
#
#    pd_list.append(arr)
#
#    progress = int((i + 1) * 100 / total)
#    last_printed = -1
#
#    if (progress % 10 == 0 ) and (progress != last_printed): 
#      print(f" ... processed {progress} %")
#      last_printed = progress
#
#  df = pd.concat(pd_list) 
# 
#  return df

def filterDf(df, selection):

  #only used when using uproot (slow af)

  #change to python expr
  selection = selection.replace("&&", " and ")
  selection = selection.replace("||", " or ")
  selection = selection.replace("!(", " ~( ") #expressions like !(mu7_ip4) -> ~(mu7_ip4)

  print(f"Filtering files with selection {selection}")
  df_sel = df.query(selection).copy()

  return df_sel


def getRdf(path, debug = None):

  chain = ROOT.TChain("tree")

  if debug:
    print("Loading less files for debugging...")
    for f in glob.glob(path)[:100]: 
      print(f"Adding file {f}")
      chain.Add(f) 
  
  else:
    chain.Add(path)

  rdf = ROOT.RDataFrame(chain)

  print(f"Rdf is well defined here with {rdf.Count().GetValue()} events")

  #need to return both, chain and rdf to keep rdf alive!! Otherwise segfault

  return rdf, chain

#######################################
# Plot ROC curve and prob. histo      #
#######################################

def plotRoc(model, X, y, roc_type = "train"):

  # convert to DMatrix format
  X = xgb.DMatrix(X)  

  #convert true values to one hot format
  label_binarizer = LabelBinarizer().fit(y)
  y_onehot = label_binarizer.transform(y) # shape: [0,0,1, ...,0]

  # Make predictions on train/test X
  y_pred = model.predict(X) # this outputs probabilites [0.12, 0.13, 0.22, ... ]

  #loop over classes and plot ROC
  classes = 6

  plt.figure()
  fig, axes = plt.subplots(1, 2, figsize=(12, 6))
  axes[0].set_xlabel('False Positive Rate')
  axes[0].set_ylabel('True Positive Rate')
  axes[0].set_title('Receiver Operating Characteristic (ROC) Curve')
  axes[0].plot([0, 1], [0, 1], linestyle='--', color='gray')  # Diagonal line

  for i in range(classes):

    fpr, tpr, _ = roc_curve(y_onehot[:,i], y_pred[:,i])
    roc_auc = auc(fpr,tpr)
    axes[0].plot(fpr, tpr, label=f'ROC class {i} (AUC = {round(roc_auc,2)})', color='blue')

  axes[0].legend()
  axes[0].grid()
    
  plt.savefig(f"{path}/roc_{roc_type}.pdf")
  plt.clf()  
  plt.close()  


def plotFeatureImportance(model):

  for n in range(nfolds):
    fig, ax = plt.subplots(figsize=(10, 12))  # Increase height
    xgb.plot_importance(model[n], ax=ax)
    plt.tight_layout()  # Adjust layout to fit labels
    plt.savefig(f"{path}/feature_importance_fold_{n}.pdf")
    plt.clf()  
    plt.close()  

#######################################
# Create weight signflip weight histo #
#######################################
  
def predictAndGetWeight(model, data, df, X_df, bdt_bins, binned_weights):

  #first, lets predict the prob for this df
  X_df = xgb.DMatrix(X_df)
  df['bdt_prob'] = model.predict(X_df)
  
  # apply the weights to our df 
  # bdt_bins are the edges, so we have to subtract 1
  masks = [(df['bdt_prob'] > bdt_bins[i]) & (df['bdt_prob'] <= bdt_bins[i+1]) for i in range(len(bdt_bins)-1)]
  
  sf_weights = []
  for i in range(len(masks)):
    # now scale all masks by the weight
    sf_weights.append(binned_weights[i] * masks[i])
  
  sf_weights = sum(sf_weights)

  #append the column
  df["sf_weights"] = sf_weights

  return df 

#######################################
# Plot loss function                  #
#######################################


def plotLoss(history, flag = ""):
  '''
    Plot the loss for training and validation sets
  '''

  loss_train = history["train"]["mlogloss"]
  loss_val   = history["eval" ]["mlogloss"]
  epochs     = range(1, len(loss_train)+1)

  plt.figure()
  plt.plot(epochs, loss_train, 'g', label='Training loss')
  plt.plot(epochs, loss_val, 'b', label='Validation loss')

  plt.xlabel('Rounds')
  plt.ylabel('Loss')
  plt.legend()
  plt.savefig(path + f'/loss{flag}.pdf')
  plt.clf()
  plt.close()

def plotHist(df_wrong,df_correct, var, bins, start,stop, flag = "", region = ""):

  # Define binning
  bins = np.linspace(start, stop, bins + 1)
  
  # Normalize histograms
  norm_factor_wrong    = len(df_wrong)   if len(df_wrong)   > 0 else 1
  norm_factor_correct  = len(df_correct) if len(df_correct) > 0 else 1
  norm_factor_weighted = sum(df_wrong['sf_weights']) if sum(df_wrong['sf_weights']) > 0 else 1
  
  # Create the figure and axis
  plt.figure(figsize=(8, 6))
  
  # Plot the first normalized histogram (df_wrong)
  plt.hist(df_wrong[var], bins=bins, label=f'Wrong sign', histtype='step', linewidth=2, density=True, color = "red")
  
  # Plot the second normalized histogram (df_correct)
  plt.hist(df_correct[var], bins=bins, label=f'Correct sign', histtype='step', linewidth=2, density=True, color = "green")
  
  # Plot the third normalized histogram (df_wrong with per event weights)
  plt.hist(df_wrong[var], bins=bins, weights=df_wrong['sf_weights']/norm_factor_weighted, label=f'Wrong sign weighted', histtype='step', linewidth=2, density=True, color = "blue")
  
  label = var
  # Labels and legend
  if "lhcb_alt" in var: label = var.replace("lhcb_alt", "xyz")
  if "reco"     in var: label = var.replace("reco", "math")
  plt.xlabel(label)
  plt.ylabel('a.u.')
  plt.title(f'{region} region', loc='left')
  plt.legend(loc='lower right', bbox_to_anchor=(1, 1), borderaxespad=0.)
  plt.savefig( dt + f"/plots{flag}/{var}_{region}.pdf")
  plt.clf()  
  plt.close()  


def plotConf(model, X_test,y_test):

  #convert into dmatrix and predict
  y_pred = model.predict(xgb.DMatrix(X_test))
 
  #rather than probabilities, we want argmax (one hot predictions)
  y_pred = np.argmax(y_pred, axis = 1)

  #get confusion matrix
  cm = confusion_matrix(y_test, y_pred)

  #normalize
  cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]


  plt.figure(figsize=(8, 6))
  seaborn.heatmap(cm, annot=True, fmt=".2f" , xticklabels=class_labels, yticklabels=class_labels, cmap="Blues")
  plt.xlabel("Predicted Label")
  plt.ylabel("True Label")
  plt.title("Confusion matrix")
  plt.tight_layout()
  plt.savefig(f"{path}/conf.pdf")
  plt.clf()  
  plt.close()  


def plotScore(model, X_test, y_test, flag = ""):

  #convert into dmatrix and predict
  y_pred = model.predict(xgb.DMatrix(X_test))

  classes = 6
  for i in range(classes):

    #only select events where true label is i
    scores = y_pred[y_test == i]


    #for these events, plot the score of all classes

    col ={0:'c',1:'g',2:'m',3:'r',4:'b',5:'k'}

    plt.figure()
    for j in range(classes):
      plt.hist(scores[:,j], bins = 31, range = (0,1), alpha = 0.5, color = col[j], histtype='stepfilled', label = f"class {j}", density=True)
    plt.legend()
    plt.savefig(f"{path}/scores_for_class_{i}{flag}.pdf")
    plt.clf()  
    plt.close()  


    plt.figure()
    #also get argmax for these events
    pred_class = np.argmax(scores, axis = 1)
    counts     = np.bincount(pred_class)

    #plot distro into histogram
    plt.bar(range(len(counts)), counts, color = col[i], alpha = 0.5, label = f"class {i}")
    plt.legend()

    plt.savefig(f"{path}/class_pred_for_class{i}.pdf")
    plt.clf()  
    plt.close()  



def doKsTest(model, X_train,X_test,y_train,y_test, single = None):

  if single is not None:
    #do a KS test only for fold <single> as a cross check
    folds = [single] 
    print(folds)
  else:
    folds = range(nfolds)
    print("I enter here !")
 
  #these will hold the average histo for each class (average over folds)
  h1_av = {}
  h2_av = {}


  classes = 6
  for i in range(classes):
    h1_av[i] = ROOT.TH1F(f'train_{i}_av', f'train_{i}_av', 30, 0, 1)
    h2_av[i] = ROOT.TH1F(f'test_{i}_av' , f'test_{i}_av' , 30, 0, 1)


  for n in folds:
  
    #convert into dmatrix and predict
    y_pred_train = model[n].predict(xgb.DMatrix(X_train[n]))#.head(50000)))
    y_pred_test  = model[n].predict(xgb.DMatrix(X_test [n]))#.head(50000)))
  
    #y_train = y_train.head(50000)
    #y_test = y_test.head(50000)
  
    for i in range(classes):
      #select only events where the true class is class {i}
      scores_train = y_pred_train[y_train[n] == i]
      scores_test  = y_pred_test [y_test[n]  == i]
  
      #only keep class i and save it
      scores_train = scores_train[:,i]
      scores_test  = scores_test [:,i]
  
      #np.savetxt(f"{path}/scoretrain_{i}.csv", scores_train, delimiter = ",")
      #np.savetxt(f"{path}/scoretest_{i}.csv" , scores_test , delimiter = ",")
  
      #fill root histos and do KS test
      h1 = ROOT.TH1F(f'train_{i}', f'train_{i}', 30, 0, 1)
      h2 = ROOT.TH1F(f'test_{i}' , f'test_{i}' , 30, 0, 1)
  
      for s_train, s_test in zip(scores_train, scores_test):
        h1.Fill(s_train)
        h2.Fill(s_test )
  
      c1=ROOT.TCanvas()
      if h1.Integral()!=0: h1.Scale(1./h1.Integral())
      if h2.Integral()!=0: h2.Scale(1./h2.Integral())
 
      # add to average
      h1_av[i].Add(h1)
      h2_av[i].Add(h2)
 
      c1.Draw()
  
      h1.GetXaxis().SetTitle("score")
      h1.GetYaxis().SetRangeUser(0, max([h2.GetMaximum(),h1.GetMaximum()])*1.6)
      h1.SetFillColor(ROOT.kGreen+2)
      h1.SetLineColor(ROOT.kGreen+2)
      h1.SetFillStyle(3345)
      h1.Draw('HIST')
      h1.Draw('E SAME')
  
      h2.SetLineColor(ROOT.kBlack)
      h2.SetMarkerStyle(8)
      h2.SetMarkerSize(0.5)
      h2.SetMarkerColor(ROOT.kBlack)
      h2.Draw('EP SAME')
  
      ks_score = h1.KolmogorovTest(h2)
      ks_value = ROOT.TPaveText(0.5, 0.76, 0.88, 0.80, 'nbNDC')
      ks_value.AddText(f'KS score of class {i} = {round(ks_score,3)}')
      ks_value.SetFillColor(0)
      ks_value.Draw('EP SAME')
      leg = ROOT.TLegend(.55,.82,.83,.88)
  
      leg.AddEntry(h1 ,'Training' ,'F' )
      leg.AddEntry(h2 ,'Validation' ,'EP' )
      leg.Draw("SAME")
      c1.SaveAs(f"{path}/KS_test_class_{i}_fold_{n}.pdf")
      c1.Close()          
      del c1  
  

  if single is None:     
    #after the fold loop we average
    for i in range(classes):
      h1_av[i].Scale(1/nfolds)
      h2_av[i].Scale(1/nfolds)
  
      #lets also draw the average
      c1=ROOT.TCanvas()
      if h1_av[i].Integral()!=0: h1_av[i].Scale(1./h1_av[i].Integral())
      if h2_av[i].Integral()!=0: h2_av[i].Scale(1./h2_av[i].Integral())
      c1.Draw()
      h1_av[i].GetXaxis().SetTitle("score")
      h1_av[i].GetYaxis().SetRangeUser(0, max([h2_av[i].GetMaximum(),h1_av[i].GetMaximum()])*1.6)
      h1_av[i].SetFillColor(ROOT.kGreen+2)
      h1_av[i].SetLineColor(ROOT.kGreen+2)
      h1_av[i].SetFillStyle(3345)
      h1_av[i].Draw('HIST')
      h2_av[i].SetLineColor(ROOT.kBlack)
      h2_av[i].SetMarkerStyle(8)
      h2_av[i].SetMarkerSize(0.5)
      h2_av[i].SetMarkerColor(ROOT.kBlack)
      h2_av[i].Draw('EP SAME')
      ks_score = h1_av[i].KolmogorovTest(h2_av[i])
      ks_value = ROOT.TPaveText(0.5, 0.76, 0.88, 0.80, 'nbNDC')
      ks_value.AddText(f'Average KS score of class {i} = {round(ks_score,3)}')
      ks_value.SetFillColor(0)
      ks_value.Draw('EP SAME')
      leg = ROOT.TLegend(.55,.82,.83,.88)
      leg.AddEntry(h1_av[i] ,'Training' ,'F' )
      leg.AddEntry(h2_av[i] ,'Validation' ,'EP' )
      leg.Draw("SAME")
      c1.SaveAs(f"{path}/KS_test_class_{i}_average.pdf")
      #print('KS score: ',ks_score, len(train_pred),len(test_pred))
      c1.Close()          
      del c1  
   

def getRegions(sigma):

  #signal region
  mlow   = dsMass_ - nSignalRegion*sigma
  mhigh  = dsMass_ + nSignalRegion*sigma

  #sideband start
  mlow2  = dsMass_ - nSidebands*sigma
  mhigh2 = dsMass_ + nSidebands*sigma

  #sideband stops
  mlow3  = mlow2  - sbWidth*sigma
  mhigh3 = mhigh2 + sbWidth*sigma

  signalRegion      = f"&& ({mlow} < phiPi_m) && (phiPi_m < {mhigh})"
  anti_signalRegion = f"&& ((({mlow3} < phiPi_m) && (phiPi_m < {mlow2})) || (({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3}))) "
  leftSB            = f"&& ({mlow3} < phiPi_m) && (phiPi_m < {mlow2})"
  rightSB           = f"&& ({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3})"

  return mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, signalRegion, anti_signalRegion, leftSB, rightSB


#######################################
# Define regions and selections       #
#######################################

#define sidebands for now
sigma = 0.009 #GeV, eyeballing


mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, signalRegion, anti_signalRegion, left_sb_cut, right_sb_cut = getRegions(sigma)

# sign conditions
wrong_sign   = "&& ((mu_charge*pi_charge > 0) && (k1_charge*k2_charge < 0))" # only flip pi-mu charge!
correct_sign = "&& ((k1_charge*k2_charge < 0) && (mu_charge*pi_charge < 0))"

# NEW: remove D+ window from training!
d_plus_window = 0.06
lambda_window = 0.06
 
k1_mu_charge = f"(k1_charge * mu_charge)"
k2_mu_charge = f"(k2_charge * mu_charge)"

d_plus_k1_wrong     = f"fabs(kPiPi_m_k1_wrong - 1.869)"
d_plus_k2_wrong     = f"fabs(kPiPi_m_k2_wrong - 1.869)"
lambda_k1_wrong     = f"fabs(pKPi_m_k1_wrong - 2.28646)"
lambda_k2_wrong     = f"fabs(pKPi_m_k2_wrong - 2.28646)"

correct_sign += f"&& (({k2_mu_charge} < 0) * ({d_plus_k2_wrong} > {d_plus_window}) + ({k1_mu_charge} < 0) * ({d_plus_k1_wrong} > {d_plus_window})) "
wrong_sign   += f"&& (({k2_mu_charge} > 0) * ({d_plus_k2_wrong} > {d_plus_window}) + ({k1_mu_charge} > 0) * ({d_plus_k1_wrong} > {d_plus_window})) "

correct_sign += f"&& (({k2_mu_charge} < 0) * ({lambda_k2_wrong} > {lambda_window}) + ({k1_mu_charge} < 0) * ({lambda_k1_wrong} > {lambda_window})) "
wrong_sign   += f"&& (({k2_mu_charge} > 0) * ({lambda_k2_wrong} > {lambda_window}) + ({k1_mu_charge} > 0) * ({lambda_k1_wrong} > {lambda_window})) "



#high mass
highMass    = f"&& (dsMu_m > {bsMass_})"
#lowMass     = f"&& (dsMu_m < {bsMass_})"
lowMass     = f""


trigger_mc   = ""
trigger_data = ""

# Lets collect everything which is not signal and not combinatorial into hb (for now)
hbSelec = " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) && (gen_match_success == 1) && (gen_same_mother == 1)"

#Sign flip selection for data

start_time = time.time()

###########################
# define classes          #
###########################
classes = 6

#class selections
sel = {}

base = baselines[args.sel]
base += f"&& (phiPi_m > 1.9)" #exclude D+ region

sel[0] = base + trigger_mc   + "&& (gen_sig == 1)"  + lowMass + correct_sign 
sel[1] = base + trigger_mc   + "&& (gen_sig == 11)" + lowMass + correct_sign 
sel[2] = base + trigger_mc   + "&& (gen_sig == 0)"  + lowMass + correct_sign 
sel[3] = base + trigger_mc   + "&& (gen_sig == 10)" + lowMass + correct_sign 
sel[4] = base + trigger_mc   + hbSelec              + lowMass + correct_sign 
sel[5] = base + trigger_data                        + lowMass + wrong_sign 


###########################
# Load data and mc sets   #
###########################

sig_path   = [
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{sig_hammer_flatNano}/*",
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{sig_hammer_flatNano2}/*",
]

#hb_path    = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/{hb_flatNanos[0]}_D_only/*"
hb_path    = [
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/{hb_flatNanos[0]}_fullBPark/*",
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/{hb_flatNanos[1]}_fullBPark/*"
]

data_path  = [ 
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/20260227_105644_D_only/*"    ,
 f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_model}/data_bph1_{bdt_model}/*"

 # PART 1
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/20260227_105644_WS_only_fullBPark/*"    ,
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/25_05_2026_11_46_22_WS_only_fullBPark/*",
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/24_05_2026_17_46_53_WS_only_fullBPark/*",
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/22_05_2026_11_24_52_WS_only_fullBPark/*",

 # PART 2
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/20260227_104615_fullBPark/*"    ,
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/29_05_2026_07_05_38_fullBPark/*",
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/27_05_2026_22_31_02_fullBPark/*",
 #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/27_05_2026_17_18_28_fullBPark/*",



]

print("bonjour")

#filter df and define the seperate classes
dfs = {}
for i in range(4):
  dfs[i]  = getDf(sig_path , branches = features + ["event", "central_w", "trigger_sf"] + more_vars + ["gen_sig", "gen_match_success"], selection = sel[i], debug = debug)
  #dfs[i]  = getDf(sig_path , branches = features + ["event", "central_w"]  , selection = sel[i], debug = debug)
  dfs[i]["is_signal"] = i

dfs[4]    = getDf(hb_path  , branches = features + ["event", "trigger_sf"] + more_vars + ["gen_sig", "gen_match_success", "gen_same_mother"]             , selection = sel[4], debug = debug)
#dfs[4]    = getDf(hb_path  , branches = features + ["event"]               , selection = sel[4], debug = debug)
dfs[4]["is_signal"] = 4

dfs[5]    = getDf(data_path, branches = features + ["event", "sf_weights"] + more_vars, selection = sel[5], debug = debug)
dfs[5]["is_signal"] = 5

for i in range(6):
  print(f"=====> We have {len(dfs[i])} events for training in class {i}")


##########################
# Define sample weights  #     
##########################

#collect the effectvie populations of all classes
eff_pop = {}

names = {0:"dstau", 1:"dsstartau", 2:"dsmu", 3:"dsstarmu", 4:"hb", 5:"comb"}

#the efffective population of signal class i is just the sum of hammer weights
for i in range(4): 
  print(f"population of class {names[i]}")
  #eff_pop[i] = (dfs[i]["central_w"] / averages[f"central_w_{names[i]}"]).sum()
  eff_pop[i] = (dfs[i]["central_w"] / averages[f"central_w_{names[i]}"] * dfs[i]["trigger_sf"]).sum()

#the efffective population of hb class is just the population (no weight)
eff_pop[4] = (dfs[4]["trigger_sf"]).sum() 
#eff_pop[4] = len(dfs[4]) 

#the efffective population of data class is just the bdt signflip weight
eff_pop[5] = (dfs[5]["sf_weights"]).sum() 
#eff_pop[5] = len(dfs[5]) 

total_pop = sum(eff_pop.values())

#define total weights

# the total weight of class i is given by:
# hammer weights * class weight (signal)
# class weight (hb)
# resp. sf weight * class weight (data)

for i in range(4): 
  dfs[i]["total_w"] = (total_pop / eff_pop[i]) * dfs[i]["central_w"] / averages[f"central_w_{names[i]}"] * dfs[i]["trigger_sf"] 
  #dfs[i]["total_w"] = (total_pop / eff_pop[i]) * dfs[i]["central_w"] / averages[f"central_w_{names[i]}"] 

dfs[4]["total_w"]   = (total_pop / eff_pop[4]) * dfs[4]["trigger_sf"]  
#dfs[4]["total_w"]   = (total_pop / eff_pop[4]) 
dfs[5]["total_w"]   = (total_pop / eff_pop[5]) * dfs[5]["sf_weights"]
 

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Loading Time: {elapsed_time:.4f} seconds")

###########################
# Prepare train/test df   #
###########################

df = pd.concat(dfs, ignore_index = True)
# features must be complete in every class; only non-feature columns may be missing
assert not df[features].isna().any().any(), \
  df[features].isna().sum()[lambda s: s > 0]


#stack all the different pandas for training/testing
#replace all nans with 1.0 (f.e. if not sharing the same columns)
df = pd.concat(dfs, ignore_index = True).fillna(1.0)
#shuffle rows
df = df.sample(frac=1, random_state=1968).reset_index(drop=True)

#split data into folds

# create cyclic list, f.e. for nfold = 5, this is:
# [0,1,2,3,4,0,1,2,3,4]
cyclic = list(range(nfolds)) * 2

X_folds = {}
y_folds = {}
w_folds = {}
model   = {}
history = {n:{} for n in range(nfolds)}

for n in range(nfolds):
  
  #trick to split into folds using event number! :D
  X_folds[n] = df[features   ] [ df["event"] % nfolds == n ]
  y_folds[n] = df["is_signal"] [ df["event"] % nfolds == n ]
  w_folds[n] = df["total_w"  ] [ df["event"] % nfolds == n ]

  #reset the pandas index
  X_folds[n] = X_folds[n].reset_index(drop=True)
  y_folds[n] = y_folds[n].reset_index(drop=True)
  w_folds[n] = w_folds[n].reset_index(drop=True)

for n in range(nfolds):

  #get train and test indices for fold n
  train_indices = [cyclic[n + (i+1)] for i in range(nfolds - 2)  ]
  test_index    = cyclic[ n + nfolds - 1]
  print(f"For fold/net {n} we train on folds {train_indices} and test on fold {test_index}")
  
  #important to re-index after concatinating several pandas array! drop = True drops the old numbering 
  X_train  = pd.concat([ X_folds[i] for i in train_indices]).reset_index(drop=True)
  y_train  = pd.concat([ y_folds[i] for i in train_indices]).reset_index(drop=True)
  w_train  = pd.concat([ w_folds[i] for i in train_indices]).reset_index(drop=True)

  X_val    = X_folds[test_index] .reset_index(drop=True)
  y_val    = y_folds[test_index] .reset_index(drop=True)
  w_val    = w_folds[test_index] .reset_index(drop=True)


  
  with open(f"{scratch_dir}/X_train_{n}.pck", "wb") as f:
    pickle.dump(X_train, f)
  with open(f"{scratch_dir}/y_train_{n}.pck", "wb") as f:
    pickle.dump(y_train, f)
  with open(f"{scratch_dir}/w_train_{n}.pck", "wb") as f:
    pickle.dump(w_train, f)
  with open(f"{scratch_dir}/X_val_{n}.pck", "wb") as f:
    pickle.dump(X_val, f)
  with open(f"{scratch_dir}/y_val_{n}.pck", "wb") as f:
    pickle.dump(y_val, f)
  with open(f"{scratch_dir}/w_val_{n}.pck", "wb") as f:
    pickle.dump(w_val, f)
  with open(f"{scratch_dir}/features.pck", "wb") as f:
    pickle.dump(features, f)

  # use .train() rather than .fit() -> allows more complex handling

