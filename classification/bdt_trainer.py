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
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *



################################################
# Logic: We train on wrong  vs right sign data #
# in the left data SB, interpret the score as  #
# a weight for the wrong sign data.            #
# Evaluation happens on the right SB as a test #
################################################


# Load dataset
files_data = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/20240724_170443/*.root"
chain = ROOT.TChain("tree")
chain.Add(files_data)
rdf_wrong   = ROOT.RDataFrame(chain)
rdf_correct = ROOT.RDataFrame(chain)

file_list = glob.glob(files_data)
print(f"=====> loading {len(file_list)} files")
#tree_name = "tree"

#define sidebands for now
sigma = 0.009 #GeV, eyeballing

maxevents = 500

#######################################
# Plot ROC curve and prob. histo      #
#######################################


def plotRoc(model, data, X, X_tt, y_tt, bdt_bins, flag = "train"):

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
  axes[1].set_title('Prediction Distribution by Target')
  axes[1].legend()
  axes[1].grid()
  
  plt.tight_layout()
  plt.savefig(f'roc_{flag}.pdf')
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


def plotLoss(history):
  '''
    Plot the loss for training and validation sets
  '''

  loss_train = history["train"]["logloss"]
  loss_val   = history["eval"]["logloss"]
  epochs     = range(1, len(loss_train)+1)

  plt.figure()
  plt.plot(epochs, loss_train, 'g', label='Training loss')
  plt.plot(epochs, loss_val, 'b', label='Validation loss')

  plt.title('Training and Validation Loss')
  plt.xlabel('Epochs')
  plt.ylabel('Loss')
  plt.legend()
  plt.savefig(f'loss.pdf')
  plt.clf()


def plotHist(df_wrong,df_correct, var, bins, start,stop, flag = ""):

  # Define binning
  bins = np.linspace(start, stop, bins + 1)
  
  # Normalize histograms
  norm_factor_wrong    = len(df_wrong)   if len(df_wrong)   > 0 else 1
  norm_factor_correct  = len(df_correct) if len(df_correct) > 0 else 1
  norm_factor_weighted = sum(df_wrong['sf_weights']) if sum(df_wrong['sf_weights']) > 0 else 1
  
  # Create the figure and axis
  plt.figure(figsize=(8, 6))
  
  # Plot the first normalized histogram (df_wrong)
  plt.hist(df_wrong[var], bins=bins, label=f'{flag} Wrong', histtype='step', linewidth=2, density=True, color = "red")
  
  # Plot the second normalized histogram (df_correct)
  plt.hist(df_correct[var], bins=bins, label=f'{flag}  Correct', histtype='step', linewidth=2, density=True, color = "green")
  
  # Plot the third normalized histogram (df_wrong with per event weights)
  plt.hist(df_wrong[var], bins=bins, weights=df_wrong['sf_weights']/norm_factor_weighted, label=f'{flag} Wrong (Weighted)', histtype='step', linewidth=2, density=True, color = "blue")
  
  # Labels and legend
  plt.xlabel(var)
  plt.ylabel('a.u.')
  plt.legend(loc='upper right', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
  plt.savefig(f"{var}_{flag}.pdf")


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
wrong_sign   = "&& ((mu_charge*pi_charge > 0) || (k1_charge*k2_charge > 0))"
correct_sign = "&& ((k1_charge*k2_charge < 0) && (mu_charge*pi_charge < 0))"

#high mass
high_mass    = f"(dsMu_m > {bsMass_})"
low_mass     = f"&& (dsMu_m < {bsMass_})"

start_time = time.time()

######################################
# wrong sign data for both sidebands #
######################################

rdf_left_wrong  = rdf_wrong.Filter(left_sb_cut  + wrong_sign + low_mass)
rdf_right_wrong = rdf_wrong.Filter(right_sb_cut + wrong_sign + low_mass)
rdf_high_wrong  = rdf_wrong.Filter(high_mass    + wrong_sign)

df_left_wrong   = rdf_left_wrong .AsNumpy()
df_right_wrong  = rdf_right_wrong.AsNumpy()
df_high_wrong   = rdf_high_wrong .AsNumpy()

df_left_wrong   = pd.DataFrame(df_left_wrong )
df_right_wrong  = pd.DataFrame(df_right_wrong)
df_high_wrong   = pd.DataFrame(df_high_wrong )

# get training events
neg_events_train =  len(df_left_wrong)
print(f"=====> We have {neg_events_train} wrong sign events for training in left SB")

########################################
# correct sign data for both sidebands #
########################################

rdf_left_correct  = rdf_correct.Filter(left_sb_cut  + correct_sign + low_mass)
rdf_right_correct = rdf_correct.Filter(right_sb_cut + correct_sign + low_mass)
rdf_high_correct  = rdf_correct.Filter(high_mass    + correct_sign)

df_left_correct   = rdf_left_correct.AsNumpy()
df_right_correct  = rdf_right_correct.AsNumpy()
df_high_correct   = rdf_high_correct.AsNumpy()

df_left_correct   = pd.DataFrame(df_left_correct)
df_right_correct  = pd.DataFrame(df_right_correct)
df_high_correct   = pd.DataFrame(df_high_correct)

# get training events
pos_events_train = len(df_left_correct)
print(f"=====> We have {pos_events_train} right sign events for training in left SB")

#balance classes: (negative is 0 and pos is 1)
weight        = len(df_left_correct) / len(df_left_wrong)
weight_double = (len(df_left_correct) + len(df_right_correct)) / (len(df_left_wrong) + len(df_right_wrong)) 

df_left_wrong  ["weights"] = weight 
df_left_correct["weights"] = 1.0

df_left_wrong   ["weights_double"] = weight_double
df_left_correct ["weights_double"] = 1.0
df_right_wrong  ["weights_double"] = weight_double 
df_right_correct["weights_double"] = 1.0

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Execution Time: {elapsed_time:.4f} seconds")
print(f"Loading Time per file: {elapsed_time/len(file_list)} seconds")


# Define features and target
features = [
    "phiPi_deltaR", 
    "kk_deltaR", 
    "dsMu_deltaR", 
    "bs_pt_lhcb_alt", 
    "q2_coll", 
    #"dsMu_m", 
    "phiPi_m", 
    #"q2_lhcb_alt", 
    #"pi_pt",
    #"mu_pt",
    #"k1_pt",
    #"k2_pt",
    #"pv_prob",
    #"sv_prob",
    #"tv_prob",
    "pt_miss_coll", #affects q2 coll 
    #"cosPiK1", 
    #"cosMuW_lhcb_alt",
    #"e_gamma"
    #"mu_pt",
    #"tv_prob",
]

print("=====> files loaded")


# train it on the left sideband only
df_left_wrong   ['target'] = 0
df_left_correct ['target'] = 1
df_right_wrong  ['target'] = 0
df_right_correct['target'] = 1

data = pd.concat([df_left_wrong, df_left_correct], ignore_index=True)
data_double = pd.concat([df_left_wrong, df_left_correct, df_right_correct, df_right_wrong], ignore_index=True)

X = data[features + ["weights"]]
X_double = data_double[features + ["weights_double"]]

X_right_wrong   = df_right_wrong  [features]
X_right_correct = df_right_correct[features]
X_left_wrong    = df_left_wrong   [features]
X_left_correct  = df_left_correct [features]
X_high_wrong    = df_high_wrong   [features]
X_high_correct  = df_high_correct [features]

y        = data['target']         #  column (0 or 1)
y_double = data_double['target']  #  column (0 or 1)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_double_train, X_double_test, y_double_train, y_double_test = train_test_split(X_double, y_double, test_size=0.2, random_state=42)

weights = X_train["weights"]
weights_double = X_double_train["weights_double"]

bdt_bins = [0.00, 0.25] + list(np.linspace(0.30, 0.75, 20).tolist()) + [0.8, 1.0]


#now only keep going with features
X_train = X_train[features]
X_test  = X_test[features]
X       = X[features]

X_double_train = X_double_train[features]
X_double_test  = X_double_test[features]
X_double       = X_double[features]

# use .train() rather than .fit() -> allows more complex handling
# prepare input as Dmatrix
dtrain = xgb.DMatrix(X_train[features], label=y_train, weight=weights)
dtest  = xgb.DMatrix(X_test[features], label=y_test)

dtrain_double = xgb.DMatrix(X_double_train[features], label=y_double_train, weight=weights_double)
dtest_double  = xgb.DMatrix(X_double_test[features], label=y_double_test)

# define the model here!
params = {
  "objective"  : "binary:logistic",
  "eval_metric": "logloss",
  "max_depth"  : 4, 
  "eta"        : 0.005
}

evals        = [(dtrain, "train"), (dtest, "eval")]
evals_double = [(dtrain_double, "train"), (dtest_double, "eval")]

#to save history
history= {}
history_double= {}

model = xgb.train(
    params,
    dtrain,
    num_boost_round=2000,
    evals=evals,
    evals_result=history,
    early_stopping_rounds=30,
    verbose_eval=True
)
model_double = xgb.train(
    params,
    dtrain_double,
    num_boost_round=2000,
    evals=evals_double,
    evals_result=history_double,
    early_stopping_rounds=30,
    verbose_eval=True
)
print(model)
model.save_model('bdt_model.json') 
#model_double.save_model('bdt_model.json') 

print("=====> training finished")

# plot roc and get the binned weights from hist ratio 
data,binned_weights = plotRoc(model, data,X, X_train, y_train, bdt_bins, flag = "train")
data,_              = plotRoc(model, data,X, X_test,  y_test,  bdt_bins, flag = "test")

data_double,binned_weights_double = plotRoc(model_double, data_double,X_double, X_double_train, y_double_train, bdt_bins, flag = "train_double")
data_double,_                     = plotRoc(model_double, data_double,X_double, X_double_test,  y_double_test,  bdt_bins, flag = "test_double")

with open("bdt_tools.json", "w") as f: 
  json.dump({"binned_weights": binned_weights.tolist(), "bdt_bins": bdt_bins, "features": features}, f)
#with open("bdt_tools.json", "w") as f: 
#  json.dump({"binned_weights": binned_weights_double.tolist(), "bdt_bins": bdt_bins, "features": features}, f)


#plot loss
plotLoss(history)
plotLoss(history_double)

# predict for all the partial df and apply weight column
predictAndGetWeight(model, data, df_right_wrong  , X_right_wrong  , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_right_correct, X_right_correct, bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_left_wrong   , X_left_wrong   , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_left_correct , X_left_correct , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_high_wrong   , X_high_wrong   , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_high_correct , X_high_correct , bdt_bins, binned_weights)

plotHist(df_right_wrong,df_right_correct, "q2_coll", 20, 0 ,12        , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "pi_pt", 20, 0 ,15          , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "kk_deltaR", 20, 0 ,0.5     , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "cosPiK1", 20, -1 ,1        , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "pt_miss_coll", 20, 0 ,30   , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "phiPi_m", 20, 1.968, 2.028 , flag = "rightSB")

plotHist(df_left_wrong,df_left_correct, "q2_coll", 20, 0 ,12        , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "pi_pt", 20, 0 ,15          , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "kk_deltaR", 20, 0 ,0.5     , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "cosPiK1", 20, -1 ,1        , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "pt_miss_coll", 20, 0 ,30   , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "phiPi_m", 20,  1.91, 1.968 , flag = "leftSB")

plotHist(df_high_wrong,df_high_correct, "q2_coll", 20, -12 ,12      , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "pi_pt", 20, 0 ,15          , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "kk_deltaR", 20, 0 ,0.5     , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "cosPiK1", 20, -1 ,1        , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "pt_miss_coll", 20, 0 ,30   , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_m", 20,  1.91, 2.028 , flag = "high mass")

predictAndGetWeight(model_double, data_double, df_right_wrong  , X_right_wrong  , bdt_bins, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_right_correct, X_right_correct, bdt_bins, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_left_wrong   , X_left_wrong   , bdt_bins, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_left_correct , X_left_correct , bdt_bins, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_high_wrong   , X_high_wrong   , bdt_bins, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_high_correct , X_high_correct , bdt_bins, binned_weights_double)

plotHist(df_right_wrong,df_right_correct, "q2_coll", 20, 0 ,12        , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "pi_pt", 20, 0 ,15          , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "kk_deltaR", 20, 0 ,0.5     , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "cosPiK1", 20, -1 ,1        , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "pt_miss_coll", 20, 0 ,30   , flag = "rightSB")
plotHist(df_right_wrong,df_right_correct, "phiPi_m", 20, 1.968, 2.028 , flag = "rightSB")

plotHist(df_left_wrong,df_left_correct, "q2_coll", 20, 0 ,12        , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "pi_pt", 20, 0 ,15          , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "kk_deltaR", 20, 0 ,0.5     , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "cosPiK1", 20, -1 ,1        , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "pt_miss_coll", 20, 0 ,30   , flag = "leftSB")
plotHist(df_left_wrong,df_left_correct, "phiPi_m", 20,  1.91, 1.968 , flag = "leftSB")

plotHist(df_high_wrong,df_high_correct, "q2_coll", 20, -12 ,12      , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "pi_pt", 20, 0 ,15          , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "kk_deltaR", 20, 0 ,0.5     , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "cosPiK1", 20, -1 ,1        , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "pt_miss_coll", 20, 0 ,30   , flag = "high mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_m", 20,  1.91, 2.028 , flag = "high mass")

#X_train['bdt_prob']          = model.predict      (X_train)
#X_train['bdt_prob'  ] = model.predict_proba(X_train[features])[:, 1]

outfile = uproot.recreate("output_corrected.root")
outfile['tree_left_wrong'         ] = df_left_wrong
outfile['tree_left_correct'       ] = df_left_correct
outfile['tree_right_wrong'        ] = df_right_wrong
outfile['tree_right_correct'      ] = df_right_correct
#outfile['X_train'                 ] = X_train 

