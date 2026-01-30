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
#model_mu7_kk   = xgb.Booster()
model_mu7_pimu = xgb.Booster()
#model_mu9_kk   = xgb.Booster()
model_mu9_pimu = xgb.Booster()

#model_mu7_kk  .load_model(f'/work/pahwagne/RDsTools/classification/{bdt_model_25_mu7}/bdt_model_kk.json')
model_mu7_pimu.load_model(f'/work/pahwagne/RDsTools/classification/{bdt_model_25_mu7}/bdt_model_pimu.json')

#model_mu9_kk  .load_model(f'/work/pahwagne/RDsTools/classification/{bdt_model_25_mu9}/bdt_model_kk.json')
model_mu9_pimu.load_model(f'/work/pahwagne/RDsTools/classification/{bdt_model_25_mu9}/bdt_model_pimu.json')

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
with open(f"/work/pahwagne/RDsTools/classification/{bdt_model_25_mu7}/bdt_tools_pimu.json", "r") as f:
    data_mu7 = json.load(f)

    binned_weights_mu7_pimu = np.array(data_mu7["binned_weights"])
    bdt_bins_mu7_pimu       = np.array(data_mu7["bdt_bins"      ])
    features_mu7_pimu       = np.array(data_mu7["features"      ])

binned_weights_mu7_pimu     = np.array(binned_weights_mu7_pimu)
bdt_bins_mu7_pimu           = np.array(bdt_bins_mu7_pimu      )
features_mu7_pimu           = np.array(features_mu7_pimu      )

#with open(f"/work/pahwagne/RDsTools/classification/{bdt_model_25_mu7}/bdt_tools_kk.json", "r") as f:
#    data_mu7 = json.load(f)
#
#    binned_weights_mu7_kk = np.array(data_mu7["binned_weights"])
#    bdt_bins_mu7_kk       = np.array(data_mu7["bdt_bins"      ])
#    features_mu7_kk       = np.array(data_mu7["features"      ])
#
#binned_weights_mu7_kk     = np.array(binned_weights_mu7_kk)
#bdt_bins_mu7_kk           = np.array(bdt_bins_mu7_kk      )
#features_mu7_kk           = np.array(features_mu7_kk      )


# open also binned weights
with open(f"/work/pahwagne/RDsTools/classification/{bdt_model_25_mu9}/bdt_tools_pimu.json", "r") as f:
    data_mu9 = json.load(f)

    binned_weights_mu9_pimu  = np.array(data_mu9["binned_weights"])
    bdt_bins_mu9_pimu        = np.array(data_mu9["bdt_bins"      ])
    features_mu9_pimu        = np.array(data_mu9["features"      ])

binned_weights_mu9_pimu      = np.array(binned_weights_mu9_pimu)
bdt_bin_mu9_pimu             = np.array(bdt_bins_mu9_pimu      )
feature_mu9_pimu             = np.array(features_mu9_pimu      )

#with open(f"/work/pahwagne/RDsTools/classification/{bdt_model_25_mu9}/bdt_tools_kk.json", "r") as f:
#    data_mu9 = json.load(f)
#
#    binned_weights_mu9_kk  = np.array(data_mu9["binned_weights"])
#    bdt_bins_mu9_kk        = np.array(data_mu9["bdt_bins"      ])
#    features_mu9_kk        = np.array(data_mu9["features"      ])
#
#binned_weights_mu9_kk      = np.array(binned_weights_mu9_kk)
#bdt_bin_mu9_kk             = np.array(bdt_bins_mu9_kk      )
#feature_mu9_kk             = np.array(features_mu9_kk      )

features = features_mu7_pimu #the same for all

# turn into dmatrix and predict score (on features column only, then append to full df)
X_df = df[features]
X_df = xgb.DMatrix(X_df)
df['bdt_prob_mu7_pimu'] = model_mu7_pimu.predict(X_df)
#df['bdt_prob_mu7_kk']   = model_mu7_kk  .predict(X_df)
df['bdt_prob_mu9_pimu'] = model_mu9_pimu.predict(X_df)
#df['bdt_prob_mu9_kk']   = model_mu9_kk  .predict(X_df)

print("used moels are:")
print(bdt_model_25_mu7)
print(bdt_model_25_mu9)

#apply weights column

#mu7 trigger (careful bc it is a float! use > 0.5)
trig_mu7_pimu = ((df["mu7_ip4"] > 0.5)                          & (df["pi_charge"] * df["mu_charge"] > 0) & (df["k1_charge"] * df["k2_charge"] < 0)) # mask for mu7 pimu 
#trig_mu7_kk   = ((df["mu7_ip4"] > 0.5)                          & (df["k1_charge"] * df["k2_charge"] > 0)) # mask for mu7 kk

#mu9 trigger (orthogonal)
trig_mu9_pimu = ((df["mu9_ip6"] > 0.5) & (df["mu7_ip4"] < 0.5) & (df["pi_charge"] * df["mu_charge"] > 0) & (df["k1_charge"] * df["k2_charge"] < 0) ) # mask for mu9 pimu
#trig_mu9_kk   = ((df["mu9_ip6"] > 0.5) & (df["mu7_ip4"] < 0.5) & (df["k1_charge"] * df["k2_charge"] > 0) ) # mask for mu9 kk

# for each event, this is a list of length bins with exactly one true [[0,0,0,1,0,0,0,0], ....]
#bins_mu7_kk    = [(df['bdt_prob_mu7_kk']    > bdt_bins_mu7_kk[i])   & (df['bdt_prob_mu7_kk']     <= bdt_bins_mu7_kk[i+1])       for i in range(len(bdt_bins_mu7_kk)-1)  ]
bins_mu7_pimu  = [(df['bdt_prob_mu7_pimu']  > bdt_bins_mu7_pimu[i]) & (df['bdt_prob_mu7_pimu']   <= bdt_bins_mu7_pimu[i+1])     for i in range(len(bdt_bins_mu7_pimu)-1)]
#bins_mu9_kk    = [(df['bdt_prob_mu9_kk']    > bdt_bins_mu9_kk[i])   & (df['bdt_prob_mu9_kk']     <= bdt_bins_mu9_kk[i+1])       for i in range(len(bdt_bins_mu9_kk)-1)  ]
bins_mu9_pimu  = [(df['bdt_prob_mu9_pimu']  > bdt_bins_mu9_pimu[i]) & (df['bdt_prob_mu9_pimu']   <= bdt_bins_mu9_pimu[i+1])     for i in range(len(bdt_bins_mu9_pimu)-1)]

weights_mu7_pimu = [bins_mu7_pimu[i] * binned_weights_mu7_pimu[i] for i in range(len(bins_mu7_pimu))]
#weights_mu7_kk   = [bins_mu7_kk[i]   * binned_weights_mu7_kk[i]   for i in range(len(bins_mu7_kk))  ]
weights_mu9_pimu = [bins_mu9_pimu[i] * binned_weights_mu9_pimu[i] for i in range(len(bins_mu9_pimu))]
#weights_mu9_kk   = [bins_mu9_kk[i]   * binned_weights_mu9_kk[i]   for i in range(len(bins_mu9_kk))  ]

# now we want to pick mask_mu9 for events which are triggered by mu9, otherwise mu7
#sf_weights = [
#     np.select([trig_mu7_pimu, trig_mu7_kk   , trig_mu9_pimu   , trig_mu9_kk], [w1, w2, w3, w4], default=1.0)
#     for w1, w2, w3, w4 in zip( weights_mu7_pimu, weights_mu7_kk, weights_mu9_pimu, weights_mu9_kk
#    )
#]

sf_weights = [
     np.select([trig_mu7_pimu, trig_mu9_pimu], [w1, w2], default=1.0)
     for w1, w2 in zip( weights_mu7_pimu, weights_mu9_pimu)
]

sf_weights = sum(sf_weights)
                                                                                                             
#append the column
df["sf_weights"] = sf_weights


# save into tree:
outfile = uproot.recreate("HOOK_FILE_OUT")
outfile["tree"] = df
