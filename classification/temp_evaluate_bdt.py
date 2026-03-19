import pdb
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
model = xgb.Booster()
model.load_model(f'/work/pahwagne/RDsTools/classification/{bdt_model}/bdt_model_pimu.json')

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
with open(f"/work/pahwagne/RDsTools/classification/{bdt_model}/bdt_tools_pimu.json", "r") as f:
    data_ = json.load(f)

    binned_weights = np.array(data_["binned_weights"])
    bdt_bins       = np.array(data_["bdt_bins"      ])
    features       = np.array(data_["features"      ])
    weight_wrong   = np.array(data_["weight_wrong"  ])

binned_weights     = np.array(binned_weights)
bdt_bins           = np.array(bdt_bins      )
features           = np.array(features      )
weight_wrong       = float(weight_wrong  )


# turn into dmatrix and predict score (on features column only, then append to full df)
X_df = df[features]
X_df = xgb.DMatrix(X_df)
df['bdt_prob'] = model.predict(X_df)

print("used moels is:")
print(bdt_model)

#pdb.set_trace()
#apply weights column

# NEW!!
s       = np.clip(df['bdt_prob'], 1e-6, 1 - 1e-6)
weights = s / (1 - s) 
df["sf_weights"] = weights

# save into tree:
outfile = uproot.recreate("HOOK_FILE_OUT")
outfile["tree"] = df
