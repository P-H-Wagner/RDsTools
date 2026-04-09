import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument("-dt","--datetime")
args   = parser.parse_args()


#load xx pandas and features from pck

df   = pd.read_pickle(f"./outputs/test_{args.datetime}/xx_train.pck"      )
feat = pd.read_pickle(f"./outputs/test_{args.datetime}/input_features.pck")


os.system(f"mkdir -p ./outputs/test_{args.datetime}/logx/")
os.system(f"mkdir -p ./outputs/test_{args.datetime}/logy/")

#load settings and get nr of folds

with open(f"./outputs/test_{args.datetime}/settings.txt", "r") as f:
  for line in f:
    #\s* looks for spaces
    #(\d+) for one or more digits
    #() what we want to capture
    match = re.search(r"Nfolds:\s*(\d+)", line)

    if match:
      nfolds = int(match.group(1))
      break

 
#patch them together into a pandas :)
df_pandas = {}
for n in range(nfolds): df_pandas[n] =  pd.DataFrame(df[n],  columns=list(set(feat)))


#plot into hist

#loop over variables
for i,var in enumerate(feat):

  print(f"====> Plotting for variable {var}") 

  #loop over folds and concatenate
  data = pd.concat([df_pandas[n][var] for n in range(1)], ignore_index=True)
 
 
  #get min and max elements
  xx_min = min(data)  
  xx_max = max(data)  

  print(f"====> plotting ...") 
  plt.figure()
  plt.xlabel(f"{var} (scaled)")
  plt.ylabel("#events")
  plt.hist(data, bins = 50, range= (xx_min, xx_max), histtype = "stepfilled")
  plt.savefig(f"./outputs/test_{args.datetime}/scaled_{var}.pdf")

  plt.yscale("log")
  plt.savefig(f"./outputs/test_{args.datetime}/logy/scaled_{var}.pdf")
  plt.clf()
  plt.close()
