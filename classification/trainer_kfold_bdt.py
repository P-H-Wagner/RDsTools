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
from root_numpy import tree2array

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

###########################
# Parse arguments         #
###########################

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--prod",    required = True, help = "Specify '24' or '25' to specify the data production")
parser.add_argument("-t", "--trigger", required = True, help = "Specify '7' or '9' to specify trigger menu")
parser.add_argument('-d', '--debug'  , action='store_true' )
args = parser.parse_args()

if args.prod not in ["24", "25"]:
  raise ValueError ("Error: Not a valid key for --prod, please use '24' or '25'")
else: prod = args.prod

if args.trigger not in ["mu7", "mu9"]:
  raise ValueError ("Error: Not a valid key for --trigger, please use 'mu7' or 'mu9'")
else: trigger = args.trigger

debug = False
if args.debug: debug = args.debug

with open("/work/pahwagne/RDsTools/hammercpp/development_branch/weights/average_weights.yaml","r") as f:
  averages = yaml.safe_load(f)

#date time
now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")
path = f"./bdt_outputs/bdt_{dt}/"
os.makedirs(path)

###########################
# Define features         #
###########################

# Define features and target
features = [
'cosPhiDs_lhcb_alt',
'cosPhiDs_reco_1',
'cosPhiDs_reco_2',
#
'abs(cosPiK1)',

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

"lxy_bs_sig",
"ds_vtx_cosine_xy_pv",
"ds_vtx_cosine_xy",
"signed_decay_ip3d_mu_ds_sv",


"e_gamma / photon_pt",
"bs_mass_corr",
"bs_mass_corr_photon",
"ds_perp",
"ds_perp_photon",
"ds_mu_perp",
"ds_mu_perp_photon",

"rel_iso_03_pv"

]

#needed for filtering, selection, etc.
more_vars = ["k1_pt","k2_pt","lxy_ds","mu_id_medium","fv_prob","mu_is_global",
"ds_vtx_cosine_xyz_pv","mu7_ip4","mu9_ip6",
"k1_charge", "k2_charge", "pi_charge", "mu_charge"
]


###########################
# Define helper functions #
###########################


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
  for name in files:
    f = ROOT.TFile(name)
    tree_obj = f.Get("tree")
    arr = tree2array(tree_obj,selection = selection, branches = branches) #event will be removed later!
    pd_list.append(pd.DataFrame(arr))

  df = pd.concat(pd_list) 
 
  return df

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
  loss_val   = history["eval"]["mlogloss"]
  epochs     = range(1, len(loss_train)+1)

  plt.figure()
  plt.plot(epochs, loss_train, 'g', label='Training loss')
  plt.plot(epochs, loss_val, 'b', label='Validation loss')

  plt.xlabel('Rounds')
  plt.ylabel('Loss')
  plt.legend()
  plt.savefig(path + f'/loss.pdf')
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


def plotScore(model, X_test, y_test):

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
    plt.savefig(f"{path}/scores_for_class{i}.pdf")


    plt.figure()
    #also get argmax for these events
    pred_class = np.argmax(scores, axis = 1)
    counts     = np.bincount(pred_class)

    #plot distro into histogram
    plt.bar(range(len(counts)), counts, color = col[i], alpha = 0.5, label = f"class {i}")
    plt.legend()

    plt.savefig(f"{path}/class_pred_for_class{i}.pdf")



def doKsTest(X_train,X_test,y_train,y_test):

  #convert into dmatrix and predict
  y_pred_train = model.predict(xgb.DMatrix(X_train.head(50000)))
  y_pred_test  = model.predict(xgb.DMatrix(X_test .head(50000)))

  y_train = y_train.head(50000)
  y_test = y_test.head(50000)

  classes = 6
  for i in range(classes):
    #select only events where the true class is class {i}
    scores_train = y_pred_train[y_train == i]
    scores_test  = y_pred_test [y_test  == i]

    #only keep class i and save it
    scores_train = scores_train[:,i]
    scores_test  = scores_test [:,i]

    np.savetxt(f"{path}/scoretrain_{i}.csv", scores_train, delimiter = ",")
    np.savetxt(f"{path}/scoretest_{i}.csv" , scores_test , delimiter = ",")

    #fill root histos and do KS test
    h1 = ROOT.TH1F(f'train_{i}', f'train_{i}', 30, 0, 1)
    h2 = ROOT.TH1F(f'test_{i}' , f'test_{i}', 30, 0, 1)

    for s_train, s_test in zip(scores_train, scores_test):
      h1.Fill(s_train)
      h2.Fill(s_test )

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
    ks_value.AddText(f'Average KS score of class {i} = {round(ks_score,3)}')
    ks_value.SetFillColor(0)
    ks_value.Draw('EP SAME')
    leg = ROOT.TLegend(.55,.82,.83,.88)

    leg.AddEntry(h1 ,'Training' ,'F' )
    leg.AddEntry(h2 ,'Validation' ,'EP' )
    leg.Draw("SAME")
    c1.SaveAs(f"{path}/KS_test_{i}.pdf")
     

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

maxevents = 500

mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, signalRegion, anti_signalRegion, left_sb_cut, right_sb_cut = getRegions(sigma)

# sign conditions
wrong_sign   = "&& ((mu_charge*pi_charge > 0) && (k1_charge*k2_charge < 0))" # only flip pi-mu charge!
correct_sign = "&& ((k1_charge*k2_charge < 0) && (mu_charge*pi_charge < 0))"

#high mass
highMass    = f"&& (dsMu_m > {bsMass_})"
lowMass     = f"&& (dsMu_m < {bsMass_})"

#trigger
if trigger == "mu7":
  trigger_data = " && (mu7_ip4 == 1)" #on data we select on the trigger
  trigger_mc   = "" #" && (event % 2 == 0) " #on mc we use event nr (all triggers are on for mc!)
else:
  trigger_data = " && ((mu9_ip6 == 1) && ( mu7_ip4 == 0)) "
  trigger_mc   = "" # "&& (event % 2 == 1) "

# Lets collect everything which is not signal and not combinatorial into hb (for now)
hbSelec = " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) && (gen_match_success)"

#Sign flip selection for data
signFlip   = " && ((k1_charge*k2_charge < 0) && (pi_charge*mu_charge>0))"

start_time = time.time()

###########################
# define classes          #
###########################
classes = 6

#class selections
sel = {}
sel[0] =  base_wout_tv_25 + trigger_mc   + "&& (gen_sig == 1)"  + signalRegion + lowMass 
sel[1] =  base_wout_tv_25 + trigger_mc   + "&& (gen_sig == 11)" + signalRegion + lowMass 
sel[2] =  base_wout_tv_25 + trigger_mc   + "&& (gen_sig == 0)"  + signalRegion + lowMass 
sel[3] =  base_wout_tv_25 + trigger_mc   + "&& (gen_sig == 10)" + signalRegion + lowMass 
sel[4] =  base_wout_tv_25 + trigger_mc   + hbSelec              + signalRegion + lowMass 
sel[5] =  base_wout_tv_25 + trigger_data + signFlip             + signalRegion + lowMass 


###########################
# Load data and mc sets   #
###########################

sig_path   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{sig_cons_hammer_25}/*"
hb_path    = [
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_25[0]}/*",
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bs_cons_25[0]}/*",
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{b0_cons_25[0]}/*",
f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bplus_cons_25[0]}/*",
]
data_path  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_25}/*"

print("bonjour")

#filter df and define the seperate classes
dfs = {}
for i in range(4):
  dfs[i]  = getDf(sig_path , branches = features + ["central_w"] , selection = sel[i], debug = debug)
  dfs[i]["is_signal"] = i
dfs[4]    = getDf(hb_path  , branches = features , selection = sel[4], debug = debug)
dfs[4]["is_signal"] = 4
dfs[5]    = getDf(data_path, branches = features + ["sf_weights"], selection = sel[5], debug = debug)
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
  eff_pop[i] = (dfs[i]["central_w"] / averages[f"central_w_{names[i]}"]).sum()

#the efffective population of hb class is just the population (no weight)
eff_pop[4] = len(dfs[4]) 
#the efffective population of data class is just the bdt signflip weight
eff_pop[5] = (dfs[5]["sf_weights"]).sum() 

total_pop = sum(eff_pop.values())

#define total weights

# the total weight of class i is given by:
# hammer weights * class weight (signal)
# class weight (hb)
# resp. sf weight * class weight (data)

for i in range(4): 
  dfs[i]["total_w"] = (total_pop / eff_pop[i]) * dfs[i]["central_w"] / averages[f"central_w_{names[i]}"] 
dfs[4]["total_w"]   = (total_pop / eff_pop[4]) 
dfs[5]["total_w"]   = (total_pop / eff_pop[5]) * dfs[5]["sf_weights"]

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Loading Time: {elapsed_time:.4f} seconds")

###########################
# Prepare train/test df   #
###########################

#stack all the different pandas for training/testing
#replace all nans with 1.0 (f.e. if not sharing the same columns)
df = pd.concat(dfs, ignore_index = True).fillna(1.0)
#shuffle rows
df = df.sample(frac=1, random_state=1968).reset_index(drop=True)

#get x and y 
X  = df[features   ]
y  = df["is_signal"]
w  = df["total_w"  ]

# Split data into training and testing sets
X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(X, y, w, test_size=0.2, random_state=1968)

# use .train() rather than .fit() -> allows more complex handling
# prepare input as Dmatrix
dtrain = xgb.DMatrix(data=X_train, label=y_train, weight=w_train) 
dtest  = xgb.DMatrix(data=X_test , label=y_test , weight=w_test )

# define the model here!

def lr_schedule(step):
    initial_lr = 0.1
    decay_rate = 0.99
    return initial_lr * (decay_rate ** step)

params = {
  #"n_jobs"     : 8,
  "num_class"  : 6,
  "objective"  : "multi:softprob",
  "eval_metric": "mlogloss",
  "max_depth"  : 6, 
  "eta"        : 0.1,
  "lambda"     : 10.0
}

evals        = [(dtrain, "train"), (dtest, "eval")]

#to save history
history= {}

rounds = 10000 
es = 30

callbacks = [xgb.callback.LearningRateScheduler(lr_schedule)]
model = xgb.train(
    params,
    dtrain,
    num_boost_round=rounds,
    evals=evals,
    evals_result=history,
    early_stopping_rounds=es,
    verbose_eval=True,
#    callbacks = callbacks
)

print("=====> training finished")

sys.exit()

with open( path + f"/info.txt", "w") as f:
  f.write( f" These plots use the following params: {params}, with {rounds} rounds and early stopping after {es}\n")

model.save_model       ( path + '/bdt_model.json') 

#only do this once since it takes forever
#y_pred_train = model.predict(xgb.DMatrix(X_train))
#y_pred_test  = model.predict(xgb.DMatrix(X_test ))

# plot roc and get the binned weights from hist ratio 
#data,binned_weights = plotRoc(model, data,X, X_train, y_train, bdt_bins, roc_type  = "train")
#data,_              = plotRoc(model, data,X, X_test,  y_test,  bdt_bins, roc_type  = "test")

#data_double,binned_weights_double = plotRoc(model_double, data_double,X_double, X_double_train, y_double_train, bdt_bins_double, flag = "_double" , roc_type = "train")
#data_double,_                     = plotRoc(model_double, data_double,X_double, X_double_test,  y_double_test,  bdt_bins_double, flag = "_double" , roc_type = "test")

with open( path + "/bdt_tools.json", "w") as f: 
  json.dump({"features": features}, f)


#plot loss
plotLoss(history)

# predict for all the partial df and apply weight column
predictAndGetWeight(model, data, df_right_wrong  , X_right_wrong  , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_right_correct, X_right_correct, bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_left_wrong   , X_left_wrong   , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_left_correct , X_left_correct , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_high_wrong   , X_high_wrong   , bdt_bins, binned_weights)
predictAndGetWeight(model, data, df_high_correct , X_high_correct , bdt_bins, binned_weights)

plotHist(df_right_wrong,df_right_correct, "q2_coll", 20, 0 ,12        , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "pi_pt", 20, 0 ,15          , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "kk_deltaR", 20, 0 ,0.5     , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "phiPi_deltaR", 20, 0 ,0.5  , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "dsMu_deltaR", 20, 0 ,1     , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "cosPiK1", 20, -1 ,1        , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "pt_miss_coll", 20, 0 ,30   , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "phiPi_m", 20, 1.968, 2.028 , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "mu_pt", 20,  0, 15         , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "mu_eta", 25,  -2.4, 2.4    , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "pi_eta", 25,  -2.4, 2.4    , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k1_pt", 20,  0, 15         , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k2_pt", 20,  0, 15         , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k1_eta", 25,  -2.4, 2.4    , region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k2_eta", 25,  -2.4, 2.4    , region = "Right SB")

plotHist(df_left_wrong,df_left_correct, "q2_coll", 20, 0 ,12          , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "bs_pt_lhcb_alt", 20, 0 ,60   , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "pi_pt", 20, 0 ,15            , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "kk_deltaR", 20, 0 ,0.5       , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "phiPi_deltaR", 20, 0 ,0.5    , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "dsMu_deltaR", 20, 0 ,1       , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "cosPiK1", 20, -1 ,1          , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "pt_miss_coll", 20, 0 ,30     , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "phiPi_m", 20,  1.91, 1.968   , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "mu_pt", 20,  0, 15           , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "mu_eta", 25,  -2.4, 2.4      , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "pi_eta", 25,  -2.4, 2.4      , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k1_pt", 20,  0, 15           , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k2_pt", 20,  0, 15           , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k1_eta", 25,  -2.4, 2.4      , region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k2_eta", 25,  -2.4, 2.4      , region = "Left SB")


plotHist(df_high_wrong,df_high_correct, "q2_coll", 20, -12 ,12        , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "bs_pt_lhcb_alt", 20, 0 ,60   , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "pi_pt", 20, 0 ,15            , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "kk_deltaR", 20, 0 ,0.5       , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_deltaR", 20, 0 ,0.5    , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "dsMu_deltaR", 20, 0 ,1       , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "cosPiK1", 20, -1 ,1          , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "pt_miss_coll", 20, 0 ,30     , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_m", 20,  1.91, 2.028   , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "mu_pt", 20,  0, 15           , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "mu_eta", 25,  -2.4, 2.4      , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "pi_eta", 25,  -2.4, 2.4      , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k1_pt", 20,  0, 15           , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k2_pt", 20,  0, 15           , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k1_eta", 25,  -2.4, 2.4      , region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k2_eta", 25,  -2.4, 2.4      , region = "High mass")

predictAndGetWeight(model_double, data_double, df_right_wrong  , X_right_wrong  , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_right_correct, X_right_correct, bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_left_wrong   , X_left_wrong   , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_left_correct , X_left_correct , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_high_wrong   , X_high_wrong   , bdt_bins_double, binned_weights_double)
predictAndGetWeight(model_double, data_double, df_high_correct , X_high_correct , bdt_bins_double, binned_weights_double)

plotHist(df_right_wrong,df_right_correct, "q2_coll", 20, 0 ,12        , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "q2_lhcb_alt", 20, 0 ,12    , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "bs_pt_lhcb_alt", 20, 0 ,60 , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "pi_pt", 20, 0 ,15          , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "kk_deltaR", 20, 0 ,0.5     , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "phiPi_deltaR", 20, 0 ,0.5  , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "dsMu_deltaR", 20, 0 ,1     , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "cosPiK1", 20, -1 ,1        , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "pt_miss_coll", 20, 0 ,30   , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "phiPi_m", 20, 1.968, 2.028 , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "mu_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "mu_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "pi_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k1_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k2_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k1_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k2_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k1_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k2_pt", 20,  0, 15         , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k1_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")
plotHist(df_right_wrong,df_right_correct, "k2_eta", 25,  -2.4, 2.4    , flag = "_double", region = "Right SB")

plotHist(df_left_wrong,df_left_correct, "q2_coll", 20, 0 ,12          , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "q2_lhcb_alt", 20, 0 ,12      , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "bs_pt_lhcb_alt", 20, 0 ,60   , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "pi_pt", 20, 0 ,15            , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "kk_deltaR", 20, 0 ,0.5       , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "phiPi_deltaR", 20, 0 ,0.5    , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "dsMu_deltaR", 20, 0 ,1       , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "cosPiK1", 20, -1 ,1          , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "pt_miss_coll", 20, 0 ,30     , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "phiPi_m", 20,  1.91, 1.968   , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "mu_pt", 20,  0, 15           , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "mu_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "pi_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k1_pt", 20,  0, 15           , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k2_pt", 20,  0, 15           , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k1_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")
plotHist(df_left_wrong,df_left_correct, "k2_eta", 25,  -2.4, 2.4      , flag = "_double", region = "Left SB")

plotHist(df_high_wrong,df_high_correct, "q2_coll", 20, -12 ,12        , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "q2_lhcb_alt", 20, 0 ,12      , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "bs_pt_lhcb_alt", 20, 0 ,60   , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "pi_pt", 20, 0 ,15            , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "kk_deltaR", 20, 0 ,0.5       , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_deltaR", 20, 0 ,0.5    , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "dsMu_deltaR", 20, 0 ,1       , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "cosPiK1", 20, -1 ,1          , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "pt_miss_coll", 20, 0 ,30     , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "phiPi_m", 20,  1.91, 2.028   , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "mu_pt", 20,  0, 15           , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "mu_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "pi_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k1_pt", 20,  0, 15           , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k2_pt", 20,  0, 15           , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k1_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")
plotHist(df_high_wrong,df_high_correct, "k2_eta", 25,  -2.4, 2.4      , flag = "_double", region = "High mass")

#X_train['bdt_prob']          = model.predict      (X_train)
#X_train['bdt_prob'  ] = model.predict_proba(X_train[features])[:, 1]

#outfile = uproot.recreate("output_corrected.root")
#outfile['tree_left_wrong'         ] = df_left_wrong
#outfile['tree_left_correct'       ] = df_left_correct
#outfile['tree_right_wrong'        ] = df_right_wrong
#outfile['tree_right_correct'      ] = df_right_correct
#outfile['X_train'                 ] = X_train 

