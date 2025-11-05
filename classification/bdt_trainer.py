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
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
from datetime import datetime
import argparse
from root_numpy import tree2array

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

################################################
# Logic: We train on wrong  vs right sign data #
# in the left data SB, interpret the score as  #
# a weight for the wrong sign data.            #
# Evaluation happens on the right SB as a test #
################################################

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--prod", required = True, help = "Specify '24' or '25' to specify the data production")
parser.add_argument("-t", "--trigger", required = True, help = "Specify '7' or '9' to specify trigger menu")
parser.add_argument('-d', '--debug'  , action='store_true' )
args = parser.parse_args()

#ipython -i bdt_trainer.py -- -p 25 -t 7 


if args.prod not in ["24", "25"]:
  raise ValueError ("Error: Not a valid key for --prod, please use '24' or '25'")
else: prod = args.prod

if args.trigger not in ["7", "9"]:
  raise ValueError ("Error: Not a valid key for --trigger, please use '7' or '9'")
else: trig = args.trigger


debug = False
if args.debug: debug = args.debug


# Load datasets (only part 1)
chain = ROOT.TChain("tree")

if prod == "24":

  data_path = []
  for f in data_cons_24:
    data_path.append(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/*") # old prod
  #files_data = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/20240724_170443/*.root" #old prod
  trigger = f""

else:

  #files_data = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{data_cons_25[0]}/*.root" # new prod
  data_path = []
  for f in data_cons_25:
    data_path.append(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/*") # new prod
  #data_path.append(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/data_31Aug2025_20h37m02s/*") # new prod

  if trig == "7":
    trigger = f"&& (mu7_ip4)  "
  else: 
    trigger = f"&& (mu9_ip6) && (!mu7_ip4)  "

print(f"Training on files: {data_path}")

#define sidebands for now
sigma = 0.009 #GeV, eyeballing

maxevents = 500

#date time
now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")

#######################################
# Load the dataframe                  #
#######################################
from concurrent.futures import ProcessPoolExecutor
def load_one_file(name, branches, selection):

    f = ROOT.TFile(name)
    tree_obj = f.Get("tree")
    arr = tree2array(tree_obj, selection=selection, branches=branches)
    return pd.DataFrame(arr)

# top-level helper
def load_file_for_pool(fn_br_sel):
    fn, branches, selection = fn_br_sel
    return load_one_file(fn, branches, selection)

def getDf(path, branches = None, selection = None, debug = None):

  print("Start loading df...")

  if isinstance(path,list):
    files = []
    for p in path:
      files += glob.glob(p)
  else:
    files = glob.glob(path)

  if debug:
    print("Loading less files for debugging...")
    files = files[:3]
    print(files)

  pd_list = []


  total_len = len(files)


  # prepare arguments for each worker
  args = [(fn, branches, selection) for fn in files]

  with ProcessPoolExecutor(max_workers=8) as exe:
      results = list(exe.map(load_file_for_pool, args))

  df = pd.concat(results, ignore_index=True)

  #for i,name in enumerate(files):

  #  if ( i% 10 == 0 and i != 0): 
  #    print(f"====> {i}/{total_len} files loaded")
  #    #concat now
  #    big = pd.concat(pd_list, ignore_index=True)
  #    pd_list = [big]


  #  f = ROOT.TFile(name)
  #  tree_obj = f.Get("tree")

  #  if branches is None:
  #      branches = [br.GetName() for br in tree_obj.GetListOfBranches()]
  #
  #
  #  #arr = tree2array(tree_obj,selection = selection, branches = branches) #event will be removed later!
  #  arr = tree2array(tree_obj, branches = branches) #event will be removed later!
  #  pd_list.append(pd.DataFrame(arr).query(selection))

  #chain = ROOT.TChain("tree")
  #for i,name in enumerate(files):

  #  chain.Add(name)

  #  if ( i% 100 == 0 and i != 0): 
  #    print(f"====> {i}/{total_len} files added")

  #    arr = tree2array(chain, branches = branches) #event will be removed later!
  #    pd_list.append(pd.DataFrame(arr).query(selection))



  #df = pd.concat(pd_list)

  return df


#######################################
# Define signal regions and sidebands #
#######################################
 
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

  signalRegion      = f"({mlow} < phiPi_m) && (phiPi_m < {mhigh})"
  anti_signalRegion = f"((({mlow3} < phiPi_m) && (phiPi_m < {mlow2})) || (({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3}))) "
  leftSB            = f"({mlow3} < phiPi_m) && (phiPi_m < {mlow2})"
  rightSB           = f"({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3})"

  return mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, signalRegion, anti_signalRegion, leftSB, rightSB

mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, signalRegion, anti_signalRegion, left_sb_cut, right_sb_cut = getRegions(sigma)

# sign conditions
#wrong_sign   = "&& ((mu_charge*pi_charge > 0) && (k1_charge*k2_charge < 0))" # only flip pi-mu charge!
wrong_sign   = "&& ((mu_charge*pi_charge > 0))" # flip one of both charges 
#wrong_sign   = "&& ((mu_charge*pi_charge > 0))" # flip one of both charges 
correct_sign = "&& ((k1_charge*k2_charge < 0) && (mu_charge*pi_charge < 0))"

#high mass          ##high mass region        #low mass region + sidebands
high_mass    = f"(  (dsMu_m > {bsMass_})  || ((dsMu_m < {bsMass_}) && ((({mlow3} < phiPi_m) && (phiPi_m < {mlow2})) || (({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3})))) )"
low_mass     = f"&& (dsMu_m < {bsMass_})"

start_time = time.time()

#######################################
# Defining training features          #
#######################################


features = [
    #"phiPi_deltaR", 
    #"kk_deltaR", 
    "dsMu_deltaR", 
    #"bs_pt_lhcb_alt", 
    #"bs_eta_lhcb_alt", 
    #"bs_phi_lhcb_alt", 
    #"bs_pt_coll", 
    #"bs_phi_lhcb_alt", 
    #"q2_lhcb_alt", 
    #"q2_coll", 
    #"dsMu_m", 
    #"phiPi_m", 
    #"q2_lhcb_alt", 
    #"pi_pt",
    #"pi_eta",
    #"mu_pt",
    #"mu_eta",
    #"k1_pt",
    #"k2_pt",
    #"pv_prob",
    #"sv_prob",
    #"tv_prob",
    #"pt_miss_coll", #affects q2 coll 
    #"cosPiK1", 
    #"cosMuW_lhcb_alt",
    #"e_gamma"
    #"e_star_coll",
    #"mu_pt",
    #"tv_prob",
]
branches = [
    "q2_coll",
    "q2_lhcb_alt",
    "bs_pt_lhcb_alt",
    "pi_pt",
    "kk_deltaR",
    "phiPi_deltaR",
    "dsMu_deltaR",
    "cosPiK1",
    "pt_miss_coll",
    "phiPi_m",
    "mu_pt",
    "mu_eta",
    "pi_eta",
    "k1_pt",
    "k2_pt",
    "k1_eta",
    "k2_eta",
    "bs_eta_lhcb_alt",
    "bs_phi_lhcb_alt",
    "bs_pt_coll",
    "dsMu_m",
    "pv_prob",
    "sv_prob",
    "tv_prob",
    "cosMuW_lhcb_alt",
    "e_gamma",
    "event",
    "mu7_ip4",
    "mu9_ip6",
    "mu_charge",
    "pi_charge",
    "k1_charge",
    "k2_charge",
    "e_star_coll"
]

######################################
# wrong sign data for both sidebands #
######################################

dfs = {}

dfs["df_left_wrong"  ] = getDf(data_path, branches = branches, selection = left_sb_cut  + wrong_sign + low_mass + trigger, debug = debug  ) 
dfs["df_right_wrong" ] = getDf(data_path, branches = branches, selection = right_sb_cut + wrong_sign + low_mass + trigger, debug = debug  ) 
#dfs["df_left_wrong"  ] = getDf(data_path, branches = branches, selection = left_sb_cut  + wrong_sign + trigger, debug = debug  ) 
#dfs["df_right_wrong" ] = getDf(data_path, branches = branches, selection = right_sb_cut + wrong_sign + trigger, debug = debug  ) 
dfs["df_high_wrong"  ] = getDf(data_path, branches = branches, selection = high_mass    + wrong_sign + trigger           , debug = debug  ) 

# get training events
neg_events_train_left  = len(dfs["df_left_wrong"])
neg_events_train_right = len(dfs["df_right_wrong"])
neg_events_train_high  = len(dfs["df_high_wrong"])
print(f"=====> We have {neg_events_train_left} wrong sign events for training in left SB")

########################################
# correct sign data for both sidebands #
########################################

dfs["df_left_correct"  ] = getDf(data_path, branches = branches, selection = left_sb_cut  + correct_sign + low_mass + trigger, debug = debug  ) 
dfs["df_right_correct" ] = getDf(data_path, branches = branches, selection = right_sb_cut + correct_sign + low_mass + trigger, debug = debug  ) 
#dfs["df_left_correct"  ] = getDf(data_path, branches = branches, selection = left_sb_cut  + correct_sign + trigger, debug = debug  ) 
#dfs["df_right_correct" ] = getDf(data_path, branches = branches, selection = right_sb_cut + correct_sign + trigger, debug = debug  ) 
dfs["df_high_correct"  ] = getDf(data_path, branches = branches, selection = high_mass    + correct_sign + trigger           , debug = debug  ) 

# get training events
pos_events_train_left  = len(dfs["df_left_correct"])
pos_events_train_right = len(dfs["df_right_correct"])
pos_events_train_high  = len(dfs["df_high_correct"])

print(f"=====> We have {pos_events_train_left} right sign events for training in left SB")

end_time = time.time()
elapsed_time = end_time - start_time
print(f"File loading time: {elapsed_time:.4f} seconds")

########################################
# Balance classes and define weights   #
########################################

#balance classes: (negative is 0 and pos is 1)
weight        = pos_events_train_left / neg_events_train_left
weight_double = (pos_events_train_left + pos_events_train_right) / (neg_events_train_left + neg_events_train_right) 
weight_high   = (pos_events_train_high) / (neg_events_train_high) 

#######################################
# Prepare datasets for training       #
#######################################

dfs["df_left_wrong"   ]["weights"]        = weight 
dfs["df_left_correct" ]["weights"]        = 1.0
dfs["df_left_wrong"   ]["weights_double"] = weight_double
dfs["df_left_correct" ]["weights_double"] = 1.0
dfs["df_right_wrong"  ]["weights_double"] = weight_double 
dfs["df_right_correct"]["weights_double"] = 1.0
dfs["df_high_wrong"   ]["weights_high"]   = weight_high
dfs["df_high_correct" ]["weights_high"]   = 1.0

# train it on the left sideband only
dfs["df_left_wrong"   ]["target"] = 0
dfs["df_left_correct" ]["target"] = 1
dfs["df_right_wrong"  ]["target"] = 0
dfs["df_right_correct"]["target"] = 1
dfs["df_high_wrong"   ]["target"] = 0
dfs["df_high_correct" ]["target"] = 1

#concatenate to a main df
data        = pd.concat([dfs["df_left_wrong"], dfs["df_left_correct"]                                                 ], ignore_index=True)
data_double = pd.concat([dfs["df_left_wrong"], dfs["df_left_correct"], dfs["df_right_correct"], dfs["df_right_wrong"] ], ignore_index=True)
data_high   = pd.concat([dfs["df_high_wrong"], dfs["df_high_correct"]                                                 ], ignore_index=True)

X        = data       [features + ["weights"]       ]
X_double = data_double[features + ["weights_double"]]
X_high   = data_high  [features + ["weights_high"  ]]

X_right_wrong   = dfs["df_right_wrong"  ][features]
X_right_correct = dfs["df_right_correct"][features]
X_left_wrong    = dfs["df_left_wrong"   ][features]
X_left_correct  = dfs["df_left_correct" ][features]
X_high_wrong    = dfs["df_high_wrong"   ][features]
X_high_correct  = dfs["df_high_correct" ][features]

y        = data       ['target']  #  column (0 or 1)
y_double = data_double['target']  #  column (0 or 1)
y_high   = data_high  ['target']  #  column (0 or 1)

# Split data into training and testing sets
X_train, X_test, y_train, y_test                             = train_test_split(X,        y,        test_size=0.2, random_state=42)
X_double_train, X_double_test, y_double_train, y_double_test = train_test_split(X_double, y_double, test_size=0.2, random_state=42)
X_high_train, X_high_test, y_high_train, y_high_test         = train_test_split(X_high,   y_high, test_size=0.2, random_state=42)

weights        = X_train["weights"]
weights_double = X_double_train["weights_double"]
weights_high = X_high_train["weights_high"]

bdt_bins        = [0.00, 0.35] + list(np.linspace(0.4, 0.8, 20).tolist()) + [0.85, 1.0]
bdt_bins_double = [0.00, 0.35] + list(np.linspace(0.4, 0.8, 20).tolist()) + [0.85, 1.0]
bdt_bins_high = [0.00, 0.35] + list(np.linspace(0.4, 0.8, 20).tolist()) + [0.85, 1.0]


#now only keep going with features
X_train        = X_train[features]
X_test         = X_test [features]
X              = X      [features]

X_double_train = X_double_train[features]
X_double_test  = X_double_test [features]
X_double       = X_double      [features]

X_high_train = X_high_train[features]
X_high_test  = X_high_test [features]
X_high       = X_high      [features]

# use .train() rather than .fit() -> allows more complex handling

# prepare input as Dmatrix
dtrain         = xgb.DMatrix(X_train[features], label=y_train, weight=weights)
dtest          = xgb.DMatrix(X_test [features], label=y_test)

dtrain_double  = xgb.DMatrix(X_double_train[features], label=y_double_train, weight=weights_double)
dtest_double   = xgb.DMatrix(X_double_test[features] , label=y_double_test)

dtrain_high  = xgb.DMatrix(X_high_train[features], label=y_high_train, weight=weights_high)
dtest_high   = xgb.DMatrix(X_high_test[features] , label=y_high_test)

# define the model here!
params = {
  "objective"  : "binary:logistic",
  "eval_metric": "logloss",
  "max_depth"  : 3, 
  "eta"        : 0.0005,
  "tree_method": "gpu_hist",
 
}

evals        = [(dtrain, "train"),        (dtest, "eval")]
evals_double = [(dtrain_double, "train"), (dtest_double, "eval")]
evals_high = [(dtrain_high, "train"), (dtest_high, "eval")]

#to save history
history= {}
history_double= {}
history_high= {}

rounds = 20000 
es = 30

model = xgb.train(
    params,
    dtrain,
    num_boost_round=rounds,
    evals=evals,
    evals_result=history,
    early_stopping_rounds=es,
    verbose_eval=True
)

rounds_double = 20000 
es_double = 30

model_double = xgb.train(
    params,
    dtrain_double,
    num_boost_round=rounds_double,
    evals=evals_double,
    evals_result=history_double,
    early_stopping_rounds=es_double,
    verbose_eval=True
)

rounds_high = 20000 
es_high = 30

model_high = xgb.train(
    params,
    dtrain_high,
    num_boost_round=rounds_high,
    evals=evals_high,
    evals_result=history_high,
    early_stopping_rounds=es_high,
    verbose_eval=True
)



if not os.path.exists(dt):
  os.makedirs(dt)
  os.makedirs(dt + "/plots")
  os.makedirs(dt + "/plots_double")
  os.makedirs(dt + "/plots_high")

with open( dt + f"/info.txt", "w") as f:
  f.write( f" These plots use the following params: {params}, with {rounds} rounds and early stopping after {es}\n")

with open( dt + f"/info_double.txt", "w") as f:
  f.write( f" These plots use the following params: {params}, with {rounds_double} rounds and early stopping after {es_double}\n")

with open( dt + f"/info_high.txt", "w") as f:
  f.write( f" These plots use the following params: {params}, with {rounds_high} rounds and early stopping after {es_high}\n")



model.save_model       ( dt + '/bdt_model.json') 
model_double.save_model( dt + '/bdt_model_double.json') 
model_high.save_model( dt + '/bdt_model_high.json') 

print("=====> training finished")

#######################################
# Plot ROC curve and prob. histo      #
#######################################


def plotRoc(model, data, X, X_tt, y_tt, bdt_bins, flag = "", roc_type = "train"):

  # X_tt = X_train or X_test
  # y_tt = y_train or y_test

  # convert to DMatrix format
  X_tt = xgb.DMatrix(X_tt)  
  X    = xgb.DMatrix(X)     

  # Make predictions on train/test X
  y_prob  = model.predict(X_tt) # this already outputs probabilites
  y_pred = (y_prob > 0.5).astype(int) #convert to 0 or 1 

  # Make predictions on complete dataframe
  data['bdt_prob'] = model.predict(X) # this already outputs probabilites

  # Evaluate model
  accuracy = accuracy_score(y_tt, y_pred)
  auc      = roc_auc_score (y_tt, y_prob)
  
  #print("=====> evaluation finished")
  
  print(f'Accuracy of {flag}: {accuracy:.4f}')
  print(f'AUC Score {flag}: {auc:.4f}')
  
  fpr, tpr, _ = roc_curve(y_tt, y_prob)
 
 
  # Plot ROC curve
  plt.figure()
  fig, axes = plt.subplots(1, 2, figsize=(12, 6))
  
  # ROC Curve
  axes[0].plot(fpr, tpr, label=f'ROC Curve (AUC = {auc:.4f})', color='blue')
  axes[0].plot([0, 1], [0, 1], linestyle='--', color='gray')  # Diagonal line
  axes[0].set_xlabel('False Positive Rate')
  axes[0].set_ylabel('True Positive Rate')
  axes[0].set_title('Receiver Operating Characteristic (ROC) Curve')
  axes[0].legend()
  axes[0].grid()
  
  
  
  # Prediction distribution
  aa = axes[1].hist(data.bdt_prob[data.target==0], bins=bdt_bins, alpha=0.5, label='Wrong sign', color='red'  , density=True)
  bb = axes[1].hist(data.bdt_prob[data.target==1], bins=bdt_bins, alpha=0.5, label='Correct sign' , color='green', density=True)

  # get weights from histogram ratio
  binned_weights = np.zeros_like(aa[0])
  
  for i in range(len(aa[0])):
    if aa[0][i] != 0:
      binned_weights[i] = (bb[0][i]/np.sum(bb[0])) / (aa[0][i]/np.sum(aa[0]))
    else: 
      print("ALERT! There are empty bins for the wrong signed histo, weight are zero!!")

  # axes[1].set_yscale("log")  # Set y-axis to logarithmic scale
  axes[1].set_xlabel('Predicted Probability')
  axes[1].set_ylabel('Frequency')
  axes[1].set_title('Prediction Distribution')
  axes[1].legend()
  axes[1].grid()
  
  plt.tight_layout()
  print("here")
  plt.savefig( dt + f'/plots{flag}/roc_{roc_type}.pdf')
  plt.clf()
  return data, binned_weights


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

  loss_train = history["train"]["logloss"]
  loss_val   = history["eval"]["logloss"]
  epochs     = range(1, len(loss_train)+1)

  plt.figure()
  plt.plot(epochs, loss_train, 'g', label='Training loss')
  plt.plot(epochs, loss_val, 'b', label='Validation loss')

  plt.xlabel('Rounds')
  plt.ylabel('Loss')
  plt.legend()
  plt.savefig(dt + f'/plots{flag}/loss.pdf')
  plt.clf()


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
  plt.close()
  return 


def plotKS(model, X_train, y_train, X_test, y_test, flag = ""):

  # we only have class 0 and class 1, plot them on one.

  y_pred_train = model.predict(xgb.DMatrix(X_train))#.head(50000)))
  y_pred_test  = model.predict(xgb.DMatrix(X_test ))#.head(50000)))


  for i in [0,1]:

    #select only events where the true class is class {i}
    scores_train = y_pred_train[y_train == i]
    scores_test  = y_pred_test [y_test  == i]

    if i == 0:
      scores_train = 1 - scores_train # get score for class 0
      scores_test  = 1 - scores_test # get score for class 0

    #np.savetxt(f"{path}/scoretrain_{i}.csv", scores_train, delimiter = ",")
    #np.savetxt(f"{path}/scoretest_{i}.csv" , scores_test , delimiter = ",")
    #fill root histos and do KS test

    c1=ROOT.TCanvas()
    h1 = ROOT.TH1F(f'train_{i}', f'train_{i}', 30, 0, 1)
    h2 = ROOT.TH1F(f'test_{i}' , f'test_{i}' , 30, 0, 1)

    for s_train, s_test in zip(scores_train, scores_test):
      h1.Fill(s_train)
      h2.Fill(s_test )

    if h1.Integral()!=0: h1.Scale(1./h1.Integral())
    if h2.Integral()!=0: h2.Scale(1./h2.Integral())

    # add to average
    h1.GetXaxis().SetTitle("score")
    h1.GetYaxis().SetRangeUser(0, max([h2.GetMaximum(),h1.GetMaximum()])*1.6)

    if i == 0:
      h1.SetFillColor(ROOT.kBlue)
      h1.SetLineColor(ROOT.kBlue)
    else:
      h1.SetFillColor(ROOT.kRed)
      h1.SetLineColor(ROOT.kRed)

    h1.SetFillStyle(3345)
    h1.Draw('HIST E')
    h1.Draw('E SAME')
    h2.SetLineColor(ROOT.kBlack)
    h2.SetMarkerStyle(8)
    h2.SetMarkerSize(0.5)
    h2.SetMarkerColor(ROOT.kBlack)
    h2.Draw('EP SAME')

    ks_score = h1.KolmogorovTest(h2)
    xpos1 = 0.11 
    xpos2 = 0.5
    ks_value = ROOT.TPaveText(xpos1, 0.76, xpos2, 0.80, 'nbNDC')
    ks_value.AddText(f'KS score of class {i} = {round(ks_score,3)}')
    ks_value.SetFillColor(0)
    ks_value.Draw('EP SAME')

    xpos1 = 0.11 
    xpos2 = 0.5
    leg = ROOT.TLegend(xpos1,.82,xpos2,.88)
    leg.AddEntry(h1 ,'Training' ,'F' )
    leg.AddEntry(h2 ,'Validation' ,'EP' )
    leg.Draw("SAME")

    c1.SaveAs( dt + f"/plots{flag}/KS_test_class{i}.pdf")
    c1.Close()
    del c1
 


# plot roc and get the binned weights from hist ratio 
data,binned_weights = plotRoc(model, data,X, X_train, y_train, bdt_bins, roc_type  = "train")
data,_              = plotRoc(model, data,X, X_test,  y_test,  bdt_bins, roc_type  = "test")

data_double,binned_weights_double = plotRoc(model_double, data_double,X_double, X_double_train, y_double_train, bdt_bins_double, flag = "_double" , roc_type = "train")
data_double,_                     = plotRoc(model_double, data_double,X_double, X_double_test,  y_double_test,  bdt_bins_double, flag = "_double" , roc_type = "test")

data_high,binned_weights_high = plotRoc(model_high, data_high,X_high, X_high_train, y_high_train, bdt_bins_high, flag = "_high" , roc_type = "train")
data_high,_                     = plotRoc(model_high, data_high,X_high, X_high_test,  y_high_test,  bdt_bins_high, flag = "_high" , roc_type = "test")

with open( dt + "/bdt_tools.json", "w") as f: 
  json.dump({"binned_weights": binned_weights.tolist(),        "bdt_bins": bdt_bins, "features": features}, f)
with open( dt + "/bdt_tools_double.json", "w") as f: 
  json.dump({"binned_weights": binned_weights_double.tolist(), "bdt_bins": bdt_bins_double, "features": features}, f)
with open( dt + "/bdt_tools_high.json", "w") as f: 
  json.dump({"binned_weights": binned_weights_high.tolist(), "bdt_bins": bdt_bins_high, "features": features}, f)


#plot loss
plotLoss(history)
plotLoss(history_double, flag = "_double")
plotLoss(history_high, flag = "_high")


#plot KS test between train and test
plotKS(model, X_train, y_train, X_test, y_test, flag = "")
plotKS(model_double, X_double_train, y_double_train, X_double_test, y_double_test, flag = "_double")
plotKS(model_high, X_high_train, y_high_train, X_high_test, y_high_test, flag = "_high")



# predict for all the partial df and apply weight column
predictAndGetWeight(model, data, dfs["df_right_wrong"]  , X_right_wrong  , bdt_bins, binned_weights)
predictAndGetWeight(model, data, dfs["df_right_correct"], X_right_correct, bdt_bins, binned_weights)
predictAndGetWeight(model, data, dfs["df_left_wrong"]   , X_left_wrong   , bdt_bins, binned_weights)
predictAndGetWeight(model, data, dfs["df_left_correct"] , X_left_correct , bdt_bins, binned_weights)
predictAndGetWeight(model, data, dfs["df_high_wrong"]   , X_high_wrong   , bdt_bins, binned_weights)
predictAndGetWeight(model, data, dfs["df_high_correct"] , X_high_correct , bdt_bins, binned_weights)

plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "q2_coll", 20, 0 ,12        , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "q2_lhcb_alt", 20, 0 ,12    , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "bs_pt_lhcb_alt", 20, 0 ,60 , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pi_pt", 20, 0 ,15          , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "kk_deltaR", 20, 0 ,0.5     , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "phiPi_deltaR", 20, 0 ,0.5  , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "dsMu_deltaR", 20, 0 ,1     , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "cosPiK1", 20, -1 ,1        , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "cosMuW_lhcb_alt", 20, -1 ,1         , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pt_miss_coll", 20, 0 ,30   , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "phiPi_m", 20, 1.968, 2.028 , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "mu_pt", 20,  0, 15         , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "mu_eta", 25,  -2.4, 2.4    , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pi_eta", 25,  -2.4, 2.4    , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_pt", 20,  0, 15         , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_pt", 20,  0, 15         , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_eta", 25,  -2.4, 2.4    , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_eta", 25,  -2.4, 2.4    , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "e_star_coll", 25,0, 3      , region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "dsMu_m", 25,0, 8           , region = "Right SB")

plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "q2_coll", 20, 0 ,12          , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "q2_lhcb_alt", 20, 0 ,12      , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "bs_pt_lhcb_alt", 20, 0 ,60   , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pi_pt", 20, 0 ,15            , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "kk_deltaR", 20, 0 ,0.5       , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "phiPi_deltaR", 20, 0 ,0.5    , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "dsMu_deltaR", 20, 0 ,1       , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "cosPiK1", 20, -1 ,1          , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "cosMuW_lhcb_alt", 20, -1 ,1           , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pt_miss_coll", 20, 0 ,30     , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "phiPi_m", 20,  1.91, 1.968   , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "mu_pt", 20,  0, 15           , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "mu_eta", 25,  -2.4, 2.4      , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pi_eta", 25,  -2.4, 2.4      , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k1_pt", 20,  0, 15           , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k2_pt", 20,  0, 15           , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k1_eta", 25,  -2.4, 2.4      , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k2_eta", 25,  -2.4, 2.4      , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "e_star_coll", 25,0, 3        , region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "dsMu_m", 25,0, 8             , region = "Left SB")


plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "q2_coll", 20, -12 ,12        , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "q2_lhcb_alt", 20, 0 ,12      , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "bs_pt_lhcb_alt", 20, 0 ,60   , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pi_pt", 20, 0 ,15            , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "kk_deltaR", 20, 0 ,0.5       , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "phiPi_deltaR", 20, 0 ,0.5    , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "dsMu_deltaR", 20, 0 ,1       , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "cosPiK1", 20, -1 ,1          , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "cosMuW_lhcb_alt", 20, -1 ,1           , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pt_miss_coll", 20, 0 ,30     , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "phiPi_m", 20,  1.91, 2.028   , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "mu_pt", 20,  0, 15           , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "mu_eta", 25,  -2.4, 2.4      , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pi_eta", 25,  -2.4, 2.4      , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k1_pt", 20,  0, 15           , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k2_pt", 20,  0, 15           , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k1_eta", 25,  -2.4, 2.4      , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k2_eta", 25,  -2.4, 2.4      , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "e_star_coll", 25,0, 3        , region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "dsMu_m", 25,0, 8             , region = "High mass")

predictAndGetWeight(model_double, data_double, dfs["df_right_wrong"]  , X_right_wrong  , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, dfs["df_right_correct"], X_right_correct, bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, dfs["df_left_wrong"]   , X_left_wrong   , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, dfs["df_left_correct"] , X_left_correct , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, dfs["df_high_wrong"]   , X_high_wrong   , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, dfs["df_high_correct"] , X_high_correct , bdt_bins_double, binned_weights_double)

plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "q2_coll", 20, 0 ,12        , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "q2_lhcb_alt", 20, 0 ,12    , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pi_pt", 20, 0 ,15          , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "kk_deltaR", 20, 0 ,0.5     , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "phiPi_deltaR", 20, 0 ,0.5  , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "dsMu_deltaR", 20, 0 ,1     , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "cosPiK1", 20, -1 ,1        , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "cosMuW_lhcb_alt", 20, -1 ,1, flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pt_miss_coll", 20, 0 ,30   , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "phiPi_m", 20, 1.968, 2.028 , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "mu_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "mu_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pi_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "e_star_coll", 25,0, 3      , flag = "_double", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "dsMu_m", 25,0, 8           , flag = "_double", region = "Right SB")

plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "q2_coll", 20, 0 ,12          , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "q2_lhcb_alt", 20, 0 ,12      , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "bs_pt_lhcb_alt", 20, 0 ,60   , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pi_pt", 20, 0 ,15            , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "kk_deltaR", 20, 0 ,0.5       , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "phiPi_deltaR", 20, 0 ,0.5    , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "dsMu_deltaR", 20, 0 ,1       , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "cosPiK1", 20, -1 ,1          , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "cosMuW_lhcb_alt", 20, -1 ,1           , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pt_miss_coll", 20, 0 ,30     , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "phiPi_m", 20,  1.91, 1.968   , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "mu_pt", 20,  0, 15           , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "mu_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pi_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k1_pt", 20,  0, 15           , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k2_pt", 20,  0, 15           , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k1_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k2_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "e_star_coll", 25,0, 3        , flag = "_double", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "dsMu_m", 25,0, 8             , flag = "_double", region = "Left SB")

plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "q2_coll", 20, -12 ,12        , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "q2_lhcb_alt", 20, 0 ,12      , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "bs_pt_lhcb_alt", 20, 0 ,60   , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pi_pt", 20, 0 ,15            , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "kk_deltaR", 20, 0 ,0.5       , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "phiPi_deltaR", 20, 0 ,0.5    , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "dsMu_deltaR", 20, 0 ,1       , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "cosPiK1", 20, -1 ,1          , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "cosMuW_lhcb_alt", 20, -1 ,1           , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pt_miss_coll", 20, 0 ,30     , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "phiPi_m", 20,  1.91, 2.028   , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "mu_pt", 20,  0, 15           , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "mu_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pi_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k1_pt", 20,  0, 15           , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k2_pt", 20,  0, 15           , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k1_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k2_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "e_star_coll", 25,0, 3        , flag = "_double", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "dsMu_m", 25,0, 8             , flag = "_double", region = "High mass")


predictAndGetWeight(model_high, data_high, dfs["df_right_wrong"]  , X_right_wrong  , bdt_bins_high, binned_weights_high)
predictAndGetWeight(model_high, data_high, dfs["df_right_correct"], X_right_correct, bdt_bins_high, binned_weights_high)
predictAndGetWeight(model_high, data_high, dfs["df_left_wrong"]   , X_left_wrong   , bdt_bins_high, binned_weights_high)
predictAndGetWeight(model_high, data_high, dfs["df_left_correct"] , X_left_correct , bdt_bins_high, binned_weights_high)
predictAndGetWeight(model_high, data_high, dfs["df_high_wrong"]   , X_high_wrong   , bdt_bins_high, binned_weights_high)
predictAndGetWeight(model_high, data_high, dfs["df_high_correct"] , X_high_correct , bdt_bins_high, binned_weights_high)

plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "q2_coll", 20, 0 ,12        , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "q2_lhcb_alt", 20, 0 ,12    , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pi_pt", 20, 0 ,15          , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "kk_deltaR", 20, 0 ,0.5     , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "phiPi_deltaR", 20, 0 ,0.5  , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "dsMu_deltaR", 20, 0 ,1     , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "cosPiK1", 20, -1 ,1        , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "cosMuW_lhcb_alt", 20, -1 ,1         , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pt_miss_coll", 20, 0 ,30   , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "phiPi_m", 20, 1.968, 2.028 , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "mu_pt", 20,  0, 15         , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "mu_eta", 25,  -2.4, 2.4    , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "pi_eta", 25,  -2.4, 2.4    , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_pt", 20,  0, 15         , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_pt", 20,  0, 15         , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_eta", 25,  -2.4, 2.4    , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_eta", 25,  -2.4, 2.4    , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_pt", 20,  0, 15         , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_pt", 20,  0, 15         , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k1_eta", 25,  -2.4, 2.4    , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "k2_eta", 25,  -2.4, 2.4    , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "e_star_coll", 25,0, 3      , flag = "_high", region = "Right SB")
plotHist(dfs["df_right_wrong"],dfs["df_right_correct"], "dsMu_m", 25,5.3, 8         , flag = "_high", region = "Right SB")

plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "q2_coll", 20, 0 ,12          , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "q2_lhcb_alt", 20, 0 ,12      , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "bs_pt_lhcb_alt", 20, 0 ,60   , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pi_pt", 20, 0 ,15            , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "kk_deltaR", 20, 0 ,0.5       , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "phiPi_deltaR", 20, 0 ,0.5    , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "dsMu_deltaR", 20, 0 ,1       , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "cosPiK1", 20, -1 ,1          , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "cosMuW_lhcb_alt", 20, -1 ,1           , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pt_miss_coll", 20, 0 ,30     , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "phiPi_m", 20,  1.91, 1.968   , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "mu_pt", 20,  0, 15           , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "mu_eta", 25,  -2.4, 2.4      , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "pi_eta", 25,  -2.4, 2.4      , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k1_pt", 20,  0, 15           , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k2_pt", 20,  0, 15           , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k1_eta", 25,  -2.4, 2.4      , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "k2_eta", 25,  -2.4, 2.4      , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "e_star_coll", 25,0, 3        , flag = "_high", region = "Left SB")
plotHist(dfs["df_left_wrong"],dfs["df_left_correct"], "dsMu_m", 25,5.3, 8           , flag = "_high", region = "Left SB")

plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "q2_coll", 20, -12 ,12        , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "q2_lhcb_alt", 20, 0 ,12      , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "bs_pt_lhcb_alt", 20, 0 ,60   , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pi_pt", 20, 0 ,15            , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "kk_deltaR", 20, 0 ,0.5       , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "phiPi_deltaR", 20, 0 ,0.5    , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "dsMu_deltaR", 20, 0 ,1       , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "cosPiK1", 20, -1 ,1          , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "cosMuW_lhcb_alt", 20, -1 ,1           , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pt_miss_coll", 20, 0 ,30     , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "phiPi_m", 20,  1.91, 2.028   , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "mu_pt", 20,  0, 15           , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "mu_eta", 25,  -2.4, 2.4      , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "pi_eta", 25,  -2.4, 2.4      , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k1_pt", 20,  0, 15           , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k2_pt", 20,  0, 15           , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k1_eta", 25,  -2.4, 2.4      , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "k2_eta", 25,  -2.4, 2.4      , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "e_star_coll", 25,0, 3        , flag = "_high", region = "High mass")
plotHist(dfs["df_high_wrong"],dfs["df_high_correct"], "dsMu_m", 25,5.3, 8           , flag = "_high", region = "High mass")



#X_train['bdt_prob']          = model.predict      (X_train)
#X_train['bdt_prob'  ] = model.predict_proba(X_train[features])[:, 1]

#outfile = uproot.recreate("output_corrected.root")
#outfile['tree_left_wrong'         ] = dfs["df_left_wrong"]
#outfile['tree_left_correct'       ] = df_left_correct
#outfile['tree_right_wrong'        ] = df_right_wrong
#outfile['tree_right_correct'      ] = df_right_correct
#outfile['X_train'                 ] = X_train 

