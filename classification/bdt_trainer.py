import xgboost as xgb
import pdb
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
import random

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


# Set: export PYTHONNOUSERSITE=1 when running on CPU with cpu_bdt environment!

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


  random.shuffle(files)
  #train only on subset
  #files = files[:500]

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
pimu_wrong   = "&& ((mu_charge*pi_charge > 0) && (k1_charge*k2_charge < 0))" # pimu wrong
kk_wrong     = "&& ((k1_charge*k2_charge > 0))"                              # kk wrong 
correct_sign = "&& ((k1_charge*k2_charge < 0) && (mu_charge*pi_charge < 0))"

#                    # high mass region      # low mass region + sidebands
#train_region = f"(  (dsMu_m > {bsMass_})  || ((dsMu_m < {bsMass_}) && ((({mlow3} < phiPi_m) && (phiPi_m < {mlow2})) || (({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3})))) ) && (bs_pt_coll>10) && (cosMuW_coll > -0.95) "
train_region = f"(  ((dsMu_m < {bsMass_}) && ((({mlow3} < phiPi_m) && (phiPi_m < {mlow2})) || (({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3})))) )  "

#train_region = f"(  (dsMu_m > {bsMass_}) ) && (bs_pt_coll>10) && (cosMuW_coll > -0.95) "
low_mass     = f"&& (dsMu_m < {bsMass_})"
high_mass    = f"(  (dsMu_m > {bsMass_})  )  "

start_time = time.time()

#######################################
# Defining training features          #
#######################################


features = [
    "phiPi_deltaR", 
    "kk_deltaR", 
    "dsMu_deltaR", 
    #"bs_eta_lhcb_alt", 
    #"bs_phi_lhcb_alt", 
    "bs_pt_coll", 
    "bs_pt_lhcb_alt", 
    "bs_pt_reco_1", 
    "bs_pt_reco_2", 
    #"bs_phi_lhcb_alt", 
    #"q2_lhcb_alt", 
    "dsMu_m", 
    "phiPi_m", 
    "q2_coll", 
    "q2_lhcb_alt", 
    "q2_reco_1", 
    "q2_reco_2", 
    "pi_pt",
    #"pi_eta",
    "mu_pt",
    #"mu_eta",
    "k1_pt",
    "k2_pt",
    #"pv_prob",
    #"sv_prob",
    #"tv_prob",
    #"pt_miss_coll", #affects q2 coll 
    "cosPiK1", 
    "cosMuW_coll",
    "cosMuW_lhcb_alt",
    "cosMuW_reco_1",
    "cosMuW_reco_2",
    #"e_gamma"
    "e_star_coll",
    "e_star_lhcb_alt",
    "e_star_reco_1",
    "e_star_reco_2",
    #"mu_pt",
    "fv_chi2",
    "tv_chi2",
    "sv_chi2",
    "lxy_ds", 
    "mu_id_medium", 
    "rel_iso_03_pv", 
    "mu_is_global", 
    "ds_vtx_cosine_xyz_pv"

]

branches = [
    "q2_coll",
    "q2_lhcb_alt",
    "q2_reco_1",
    "q2_reco_2",
    "bs_pt_lhcb_alt",
    "pi_pt",
    "kk_deltaR",
    "kk_m",
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
    "bs_pt_reco_1",
    "bs_pt_reco_2",
    "dsMu_m",
    "fv_chi2",
    "sv_chi2",
    "tv_chi2",
    "cosMuW_coll",
    "cosMuW_lhcb_alt",
    "cosMuW_reco_1",
    "cosMuW_reco_2",
    "e_gamma",
    "event",
    "mu7_ip4",
    "mu9_ip6",
    "mu_charge",
    "pi_charge",
    "k1_charge",
    "k2_charge",
    "e_star_coll",
    "e_star_lhcb_alt",
    "e_star_reco_1",
    "e_star_reco_2",
    "lxy_ds", 
    "mu_id_medium", 
    "rel_iso_03_pv", 
    "mu_is_global", 
    "ds_vtx_cosine_xyz_pv"

]

#############
# Load data #
#############

dfs = {}

# WRONG SIGN EVENTS
df_wrong   = getDf(data_path, branches = branches, selection = train_region + pimu_wrong   + trigger, debug = debug  ) 

#split wrong sign events here to get control region
df_wrong_active  = df_wrong[df_wrong["event"].astype(int) % 2 == 0]
df_wrong_control = df_wrong[df_wrong["event"].astype(int) % 2 == 1]

# CORRECT SIGN EVENTS
df_correct  = getDf(data_path, branches = branches, selection = train_region + correct_sign + trigger, debug = debug  ) 

# LOADING TIME
end_time = time.time()
elapsed_time = end_time - start_time
print(f"File loading time: {elapsed_time:.4f} seconds")

#append target column
df_wrong_active["target"] = 0
df_correct     ["target"] = 1

#concatenate wrong and correct data into one for the training :D
df_train = pd.concat( [df_wrong_active, df_correct], ignore_index=True)

#keep only features for training
X = df_train[features]
#keep only target for training
y = df_train["target"]  #  column (0 or 1)

# this is to do the plots later :)
X_control = df_wrong_control[features]
X_correct = df_correct[features]

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

events_train_wrong   = len(X_train[y_train == 0])
print(f"=====> We have {events_train_wrong} pimu wrong sign events ")

events_train_correct = len(X_train[y_train == 1])
print(f"=====> We have {events_train_correct} right sign events for training ")

#balance classes: (negative is 0 and pos is 1)
weight_wrong   = events_train_correct / events_train_wrong 
#assign non-zero weight (weight_wrong) to wrong class (target 0), otherwise 1.0 for correct class (target 1)
X_train["weights"] = np.where(y_train == 0, weight_wrong, 1.0)

# define binning for weights
bdt_bins = [0.00, 0.4] + list(np.linspace(0.45, 0.6, 20).tolist()) + [0.65, 1.0]
 
#now only keep going with features and drop weights
#X_pimu_train = X_pimu_train[features]
#X_pimu_test  = X_pimu_test [features]
#X_pimu       = X_pimu      [features]

# use .train() rather than .fit() -> allows more complex handling

# prepare input as Dmatrix
dtrain  = xgb.DMatrix(X_train[features], label=y_train, weight=X_train["weights"])
dtest   = xgb.DMatrix(X_test [features], label=y_test)


# define the model here!
params = {
  "objective"  : "binary:logistic",
  "eval_metric": "logloss",
  "max_depth"  : 3, 
  "eta"        : 0.0005,
  "tree_method": "gpu_hist",
}

evals    = [(dtrain, "train"),   (dtest, "eval")]

#to save history
history = {}

rounds = 10000 
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


if not os.path.exists(dt):
  os.makedirs(dt)
  os.makedirs(dt + "/plots_pimu")
  os.makedirs(dt + "/plots_pimu_control")

with open( dt + f"/info_pimu.txt", "w") as f:
  f.write( f" These plots use the following params: {params}, with {rounds} rounds and early stopping after {es}\n")
  f.write( f" These plots use trigger: {trig}\n")

model.save_model( dt + '/bdt_model_pimu.json') 

print("=====> training finished")

#######################################
# Plot ROC curve and prob. histo      #
#######################################

def plotRoc(model, X_tt, y_tt, bdt_bins, flag = "", roc_type = "train"):


  # convert to DMatrix format
  X_tt_dmat = xgb.DMatrix(X_tt[features])  
  #X    = xgb.DMatrix(X)     

  # Make predictions on train/test X
  y_prob  = model.predict(X_tt_dmat) # this already outputs probabilites
  y_pred  = (y_prob > 0.5).astype(int) #convert to 0 or 1 

  #append column 
  X_tt["bdt_prob"] = y_prob 

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
  aa = axes[1].hist(X_tt.bdt_prob[y_tt==0], bins=bdt_bins, alpha=0.5, label='Wrong sign',    color='red'  , density=True)
  bb = axes[1].hist(X_tt.bdt_prob[y_tt==1], bins=bdt_bins, alpha=0.5, label='Correct sign' , color='green', density=True)

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
  return X_tt, binned_weights


#######################################
# Create weight signflip weight histo #
#######################################
  
def predictAndGetWeight(model, X_df, bdt_bins, binned_weights):

  X_df = X_df.copy()

  #first, lets predict the prob for this df
  X_df_dmat  = xgb.DMatrix(X_df[features])
  X_df["bdt_prob"] = model.predict(X_df_dmat)
  
  # apply the weights to our df 
  # bdt_bins are the edges, so we have to subtract 1
  masks = [(X_df['bdt_prob'] > bdt_bins[i]) & (X_df['bdt_prob'] <= bdt_bins[i+1]) for i in range(len(bdt_bins)-1)]
  
  sf_weights = []
  for i in range(len(masks)):
    # now scale all masks by the weight
    sf_weights.append(binned_weights[i] * masks[i])
  
  sf_weights = sum(sf_weights)

  #append the column
  #X_df["sf_weights"] = sf_weights

  #NEW
  s = np.clip(X_df["bdt_prob"], 1e-6, 1 - 1e-6)
  #X_df["sf_weights"] = (s / (1 - s)) / weight_wrong
  X_df["sf_weights"] = (s / (1 - s)) 

  pdb.set_trace()

  return X_df 

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

  y_pred_train = model.predict(xgb.DMatrix(X_train[features]))#.head(50000)))
  y_pred_test  = model.predict(xgb.DMatrix(X_test [features]))#.head(50000)))


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

X_train, binned_weights      = plotRoc(model,  X_train, y_train, bdt_bins, flag = "_pimu" , roc_type = "train")
X_test , _                   = plotRoc(model,  X_test , y_test , bdt_bins, flag = "_pimu" , roc_type = "test")

with open( dt + "/bdt_tools_pimu.json", "w") as f: 
  json.dump({"binned_weights": binned_weights.tolist(), "bdt_bins": bdt_bins, "features": features, "weight_wrong":weight_wrong}, f)

#plot loss
plotLoss(history, flag = "_pimu")

#plot KS test between train and test
plotKS(model, X_train, y_train, X_test, y_test, flag = "_pimu")

# predict for all the partial df and apply weight column
df_wrong_active  = predictAndGetWeight(model, df_wrong_active  ,  bdt_bins, binned_weights)
df_wrong_control = predictAndGetWeight(model, df_wrong_control ,  bdt_bins, binned_weights)
df_correct       = predictAndGetWeight(model, df_correct       ,  bdt_bins, binned_weights)

plotHist(df_wrong_active ,df_correct, "q2_coll", 20, 0 ,12        , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "q2_lhcb_alt", 20, 0 ,12    , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "q2_reco_1", 20, 0 ,12      , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "q2_reco_2", 20, 0 ,12      , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "bs_pt_coll", 20, 0 ,30     , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "bs_pt_lhcb_alt", 20, 0 ,30 , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "bs_pt_reco_2", 20, 0 ,30   , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "bs_pt_reco_1", 20, 0 ,30   , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "pi_pt", 20, 0 ,15          , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "kk_deltaR", 20, 0 ,0.3     , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "kk_m", 20, 1.0 , 1.035     , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "cosPiK1", 20, -1 ,1        , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "cosMuW_lhcb_alt", 20, -1 ,1, flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "phiPi_m", 20, 1.91, 2.028  , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "mu_pt", 20,  7, 15         , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "pi_pt", 20,  0, 6          , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "k1_pt", 20,  0, 6          , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "k2_pt", 20,  0, 6          , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "mu_eta", 25,  -2.4, 2.4    , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "pi_eta", 25,  -2.4, 2.4    , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "k1_eta", 25,  -2.4, 2.4    , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "k2_eta", 25,  -2.4, 2.4    , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "e_star_coll", 25,0, 3      , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "e_star_lhcb_alt", 25,0, 3  , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "e_star_reco_1", 25,0, 3    , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "e_star_reco_2", 25,0, 3    , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "dsMu_m", 25,0, 8           , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "dsMu_m", 25,0, 8           , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "sv_chi2", 25,0, 10         , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "tv_chi2", 25,0,  7         , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "fv_chi2", 25,0,  5         , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "lxy_ds", 25,0,  1          , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "mu_id_medium", 3,0,  1     , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "rel_iso_03_pv", 25,0,0.3   , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "mu_is_global", 3,0,1.0     , flag = "_pimu", region = "")
plotHist(df_wrong_active ,df_correct, "ds_vtx_cosine_xyz_pv", 25,0.8,1.0 , flag = "_pimu", region = "")


plotHist(df_wrong_control ,df_correct, "q2_coll", 20, 0 ,12        , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "q2_lhcb_alt", 20, 0 ,12    , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "q2_reco_1", 20, 0 ,12      , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "q2_reco_2", 20, 0 ,12      , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "bs_pt_coll", 20, 0 ,30     , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "bs_pt_lhcb_alt", 20, 0 ,30 , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "bs_pt_reco_2", 20, 0 ,30   , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "bs_pt_reco_1", 20, 0 ,30   , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "pi_pt", 20, 0 ,15          , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "kk_deltaR", 20, 0 ,0.3     , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "kk_m", 20, 1.0 , 1.035     , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "cosPiK1", 20, -1 ,1        , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "cosMuW_lhcb_alt", 20, -1 ,1, flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "phiPi_m", 20, 1.91, 2.028  , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "mu_pt", 20,  7, 15         , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "pi_pt", 20,  0, 6          , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "k1_pt", 20,  0, 6          , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "k2_pt", 20,  0, 6          , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "mu_eta", 25,  -2.4, 2.4    , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "pi_eta", 25,  -2.4, 2.4    , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "k1_eta", 25,  -2.4, 2.4    , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "k2_eta", 25,  -2.4, 2.4    , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "e_star_coll", 25,0, 3      , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "e_star_lhcb_alt", 25,0, 3  , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "e_star_reco_1", 25,0, 3    , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "e_star_reco_2", 25,0, 3    , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "dsMu_m", 25,0, 8           , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "dsMu_m", 25,0, 8           , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "sv_chi2", 25,0, 10         , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "tv_chi2", 25,0,  7         , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "fv_chi2", 25,0,  5         , flag = "_pimu_control", region = "")
plotHist(df_wrong_control ,df_correct, "lxy_ds", 25,0,  1          , flag = "_pimu", region = "")
plotHist(df_wrong_control ,df_correct, "mu_id_medium", 3,0,  1     , flag = "_pimu", region = "")
plotHist(df_wrong_control ,df_correct, "rel_iso_03_pv", 25,0,0.3   , flag = "_pimu", region = "")
plotHist(df_wrong_control ,df_correct, "mu_is_global", 3,0,1.0     , flag = "_pimu", region = "")
plotHist(df_wrong_control ,df_correct, "ds_vtx_cosine_xyz_pv", 25,0.8,1.0 , flag = "_pimu", region = "")

#predictAndGetWeight(model_kk, data_kk, dfs["df_kk_wrong"]  , X_kk_wrong      , bdt_bins_kk, binned_weights_kk)
#predictAndGetWeight(model_kk, data_kk, dfs["df_correct"]     , X_correct     , bdt_bins_kk, binned_weights_kk)
#
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "q2_coll", 20, 0 ,12        , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "q2_lhcb_alt", 20, 0 ,12    , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "q2_reco_1", 20, 0 ,12      , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "q2_reco_2", 20, 0 ,12      , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "bs_pt_coll", 20, 0 ,30     , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "bs_pt_lhcb_alt", 20, 0 ,30 , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "bs_pt_reco_1", 20, 0 ,30   , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "bs_pt_reco_2", 20, 0 ,30   , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "pi_pt", 20, 0 ,15          , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "kk_deltaR", 20, 0 ,0.5     , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "kk_m", 20, 1.0 , 1.035     , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "phiPi_deltaR", 20, 0 ,0.5  , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "dsMu_deltaR", 20, 0 ,1     , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "cosPiK1", 20, -1 ,1        , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "cosMuW_coll"    , 20, -1 ,1, flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "cosMuW_lhcb_alt", 20, -1 ,1, flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "cosMuW_reco_1"  , 20, -1 ,1, flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "cosMuW_reco_2"  , 20, -1 ,1, flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "phiPi_m", 20, 1.91, 2.028  , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "mu_pt", 20,  7, 15         , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "pi_pt", 20,  0, 6          , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "k1_pt", 20,  0, 6          , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "k2_pt", 20,  0, 6          , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "mu_eta", 25,  -2.4, 2.4    , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "pi_eta", 25,  -2.4, 2.4    , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "k1_eta", 25,  -2.4, 2.4    , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "k2_eta", 25,  -2.4, 2.4    , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "e_star_coll", 25,0, 3      , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "e_star_lhcb_alt", 25,0, 3  , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "e_star_reco_1", 25,0, 3    , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "e_star_reco_2", 25,0, 3    , flag = "_kk", region = "")
#plotHist(dfs["df_kk_wrong"],dfs["df_correct"], "dsMu_m", 25,0, 8           , flag = "_kk", region = "")

#X_train['bdt_prob']          = model.predict      (X_train)
#X_train['bdt_prob'  ] = model.predict_proba(X_train[features])[:, 1]

#outfile = uproot.recreate("output_corrected.root")
#outfile['tree_left_wrong'         ] = dfs["df_left_wrong"]
#outfile['tree_left_correct'       ] = df_left_correct
#outfile['tree_right_wrong'        ] = df_right_wrong
#outfile['tree_right_correct'      ] = df_right_correct
#outfile['X_train'                 ] = X_train 

