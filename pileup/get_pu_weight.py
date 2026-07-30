import numpy as np
import uproot
import matplotlib.pyplot as plt


#import mc signal hist, directly as list with 100 entries
#https://cmsweb.cern.ch/couchdb/reqmgr_config_cache/9317e4ada7a7ec55c5242871e82ab2d1/configFile
from mc_simu_2024_pileup import mc_pu as nmc   


#import mu7 data (get total histo)
#fdata = uproot.open("pileup_2018D_HLT_Mu7_IP4_GOLDEN_combined.root")
fdata = uproot.open("pileup_GOLDEN_combined.root")
hdata      = fdata["total"]
hdata_up   = fdata["total_up"]
hdata_down = fdata["total_down"]

ndata      = list(hdata.values()) #list with 200 entries
ndata_up   = list(hdata_up.values()) #list with 200 entries
ndata_down = list(hdata_down.values()) #list with 200 entries

#get the weight as a per bin ratio

nbins = min([len(ndata),len(nmc)])

weights      = []
weights_up   = []
weights_down = []

for i in range(nbins):

  weights     .append(ndata     [i] / nmc[i])
  weights_up  .append(ndata_up  [i] / nmc[i])
  weights_down.append(ndata_down[i] / nmc[i])

#now save the weights into a csv:
np.savetxt("pu_weights.txt"     , weights     , fmt="%.6f")
np.savetxt("pu_weights_up.txt"  , weights_up  , fmt="%.6f")
np.savetxt("pu_weights_down.txt", weights_down, fmt="%.6f")

#plot to check
plt.figure()
x = np.arange(1,101)
plt.bar(x, ndata[:nbins], color = "k",    alpha = 0.5, label = "Data - BParking 2018")
plt.bar(x, nmc  [:nbins], color = "blue", alpha = 0.5, label = "MC")
plt.xlabel("Pileup")
plt.ylabel("a.u.")
#plt.bar(x, [nmc  [:nbins][i] * weights[i] for i in range(nbins)], color = "pink", alpha = 0.5)
plt.legend()
plt.savefig("pu_distros.pdf")


plt.figure()
x = np.arange(1,101)
plt.scatter(x, weights     , color = "k", marker = "+", label = "Central")
plt.scatter(x, weights_up  , color = "r", marker = "+", label = "Up")
plt.scatter(x, weights_down, color = "b", marker = "+", label = "Down")
plt.xlabel("Pileup")
plt.ylabel("Weight")
plt.legend()
plt.savefig("pu_weights.pdf")
