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
from datetime import datetime
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

#import the model
import xgboost as xgb

now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")

print("Bonjour!")

#load model
model     = xgb.Booster()
model_mu9 = xgb.Booster()
model     .load_model(f'/work/pahwagne/RDsTools/classification/{bdt_model_25_mu7}/HOOK_MODEL.json')
model_mu9 .load_model(f'/work/pahwagne/RDsTools/classification/{bdt_model_25_mu9}/HOOK_MODEL.json')
#model     .load_model(f'/work/pahwagne/RDsTools/classification/{bdt_data_25_mu7}/bdt_model_double.json')
#model_mu9 .load_model(f'/work/pahwagne/RDsTools/classification/{bdt_data_25_mu9}/bdt_model_double.json')

# Load datasets into chain

chain = ROOT.TChain("tree")
chain.Add("HOOK_FILE_IN")
#chain.Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/20250227_155416/skimmed_base_wout_tv_25_20250227_155416_chunk_51.root")
#chain.Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/20250227_155416/skimmed_base_wout_tv_25_20250227_155416_chunk_45.root")
#chain.Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/20250227_155416/skimmed_base_wout_tv_25_20250227_155416_chunk_46.root")
#chain.Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/20250227_155416/skimmed_base_wout_tv_25_20250227_155416_chunk_44.root")
#chain.Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/20250227_155416/skimmed_base_wout_tv_25_20250227_155416_chunk_38.root")
#chain.Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/20250227_155416/skimmed_base_wout_tv_25_20250227_155416_chunk_49.root")

rdf = ROOT.RDataFrame(chain)
rdf = rdf.AsNumpy()
df  = pd.DataFrame(rdf)

# open also binned weights
with open(f"/work/pahwagne/RDsTools/classification/{bdt_model_25_mu7}/HOOK_TOOLS.json", "r") as f:
    data = json.load(f)

    binned_weights = np.array(data["binned_weights"])
    bdt_bins       = np.array(data["bdt_bins"      ])
    features       = np.array(data["features"      ])

binned_weights = np.array(binned_weights)
bdt_bins       = np.array(bdt_bins      )
features       = np.array(features      )

# open also binned weights
with open(f"/work/pahwagne/RDsTools/classification/{bdt_model_25_mu9}/HOOK_TOOLS.json", "r") as f:
    data_mu9 = json.load(f)

    binned_weights_mu9  = np.array(data_mu9["binned_weights"])
    bdt_bins_mu9        = np.array(data_mu9["bdt_bins"      ])
    features_mu9        = np.array(data_mu9["features"      ])

binned_weights_mu9      = np.array(binned_weights_mu9)
bdt_bin_mu9             = np.array(bdt_bins_mu9      )
feature_mu9             = np.array(features_mu9      )


# turn into dmatrix and predict score (on features column only, then append to full df)
X_df = df[features]
X_df = xgb.DMatrix(X_df)
df['bdt_prob']     = model.predict(X_df)
df['bdt_prob_mu9'] = model_mu9.predict(X_df)

#apply weights column

# for each event, this is a list of length bins with exactly one true [[0,0,0,1,0,0,0,0], ....]
trig_mu7 = (df["mu7_ip4"] == 1.0)                             #mask for mu7 (they are orthogonal!)
trig_mu9 = ((df["mu9_ip6"] == 1.0) & (df["mu7_ip4"] != 1.0) ) #mask for mu9

bins     = [(df['bdt_prob']     > bdt_bins[i]) & (df['bdt_prob']     <= bdt_bins[i+1]) for i in range(len(bdt_bins)-1)]
bins_mu9 = [(df['bdt_prob_mu9'] > bdt_bins[i]) & (df['bdt_prob_mu9'] <= bdt_bins[i+1]) for i in range(len(bdt_bins)-1)]

weights     = [bins[i]     * binned_weights[i]     for i in range(len(bins))]
weights_mu9 = [bins_mu9[i] * binned_weights_mu9[i] for i in range(len(bins))]

# now we want to pick mask_mu9 for events which are triggered by mu9, otherwise mu7
sf_weights  = [ w.where(trig_mu7, w_mu9) for w,w_mu9 in zip(weights, weights_mu9)]

sf_weights = sum(sf_weights)
                                                                                                             
#append the column
df["sf_weights"] = sf_weights


# save into tree:
outfile = uproot.recreate("HOOK_FILE_OUT")
outfile["tree"] = df
