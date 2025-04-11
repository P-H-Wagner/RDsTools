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

model = xgb.Booster()
model.load_model('/work/pahwagne/RDsTools/classification/bdt_model.json')


# Load datasets into chain

chain = ROOT.TChain("tree")

chain.Add("HOOK_FILE_IN")

rdf = ROOT.RDataFrame(chain)
rdf = rdf.AsNumpy()
df  = pd.DataFrame(rdf)

# open also binned weights
with open("/work/pahwagne/RDsTools/classification/bdt_tools.json", "r") as f:
    data = json.load(f)

    binned_weights = np.array(data["binned_weights"])
    bdt_bins       = np.array(data["bdt_bins"      ])
    features       = np.array(data["features"      ])

binned_weights = np.array(binned_weights)
bdt_bins       = np.array(bdt_bins      )
features       = np.array(features      )


# turn into dmatrix and predict score (on features column only, then append to full df)
X_df = df[features]
X_df = xgb.DMatrix(X_df)
df['bdt_prob'] = model.predict(X_df)

#apply weights column

masks = [(df['bdt_prob'] > bdt_bins[i]) & (df['bdt_prob'] <= bdt_bins[i+1]) for i in range(len(bdt_bins)-1)]

sf_weights = []
for i in range(len(masks)):
  # now scale all masks by the weight
  sf_weights.append(binned_weights[i] * masks[i])

sf_weights = sum(sf_weights)
                                                                                                             
#append the column
df["sf_weights"] = sf_weights


# save into tree:
outfile = uproot.recreate("HOOK_FILE_OUT")
outfile["tree"] = df
