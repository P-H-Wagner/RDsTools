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
model = {}
nfolds = 10

for n in range(nfolds):
  model[n]     = xgb.Booster()
  model[n]     .load_model(f'/work/pahwagne/RDsTools/classification/bdt_outputs/bdt_HOOK_DATE_TIME/bdt_model_{n}.json')

# get features (can choose fold 0 wlog) 
with open(f"/work/pahwagne/RDsTools/classification/bdt_outputs/bdt_HOOK_DATE_TIME/bdt_tools_0.json", "r") as f:
    data = json.load(f)
    features       = np.array(data["features"])

# Load datasets into chain

#chain = ROOT.TChain("tree")
#chain.Add("HOOK_FILE_IN")

files =  [
"HOOK_FILE_IN",
]

df_list = []

chain = ROOT.TChain("tree")
for f in files:
    chain.Add(f)

rdf = ROOT.RDataFrame(chain)

cols = rdf.GetColumnNames()

arr = rdf.AsNumpy(cols)
df = pd.DataFrame(arr)

#for i,name in enumerate(files):
#  f = ROOT.TFile(name)
#  tree_obj = f.Get("tree")
#
#  branches = [branch.GetName() for branch in tree_obj.GetListOfBranches()]
#  branches += [var for var in features if var not in branches]
#
#  arr = tree2array(tree_obj, branches = branches) #event will be removed later!
#  df_list.append(pd.DataFrame(arr))
#
#df = pd.concat(df_list)

#rdf = ROOT.RDataFrame(chain)
#rdf = rdf.AsNumpy()
#df  = pd.DataFrame(rdf)

# turn into dmatrix and predict score (on features column only, then append to full df)

to_save = []
for n in range(nfolds):

  df_slice = df [df["event"] % nfolds == n].reset_index(drop=True)

  #predict for every fold, and then append

  X_df = df_slice[features] 
  X_df = xgb.DMatrix(X_df)
  y_pred = model[n].predict(X_df)

  classes = 6
  for i in range(classes):
    df_slice[f"score{i}"] = y_pred[:,i]
    df_slice["class"    ] = np.argmax(y_pred, axis = 1)

  to_save.append(df_slice)


# save into tree:
#outfile = uproot.recreate("HOOK_FILE_OUT")
#outfile["tree"] = pd.concat(to_save)

with uproot.recreate("HOOK_FILE_OUT") as f:
  f["tree"] = pd.concat(to_save)
    




