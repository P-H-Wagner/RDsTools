from uncertainties import ufloat
import matplotlib.pyplot as plt
import sys
import numpy as np

from fs_over_fu   import bins   as bins_cms
from fs_over_fu   import values as values_cms 
from fs_over_fu   import unc    as unc_cms 

from fl_over_fufd import bins   as bins_lhcb 
from fl_over_fufd import values as values_lhcb
from fl_over_fufd import unc    as unc_lhcb


#we take the same binning for fs

#sanity check first

if len(bins_cms) != len(bins_lhcb)  : 
  print(f"Matching error! Bins of CMS and LHCb not the same size"); sys.exit();

if len(values_cms) != len(values_lhcb)  : 
  print(f"Matching error! Values of CMS and LHCb not the same size"); sys.exit();

if len(unc_cms) != len(unc_lhcb)  : 
  print(f"Matching error! Uncertainties of CMS and LHCb not the same size"); sys.exit();


for (i,j) in zip(bins_cms,bins_lhcb): 
  if (i != j): 
    print(f"Matching error! Bin edges of CMS and LHCb are different !", i, " vs. ", j );sys.exit();

#put them in a u float  and rename
# alpha = fs/fu (CMS)
# beta  = f_lambda / (fu + fd) (LHCb)

alpha = []
beta = []

for val,unc in zip (values_cms, unc_cms):
  alpha.append(ufloat(val,unc))

for val,unc in zip (values_lhcb, unc_lhcb):
  beta.append(ufloat(val,unc))

print(alpha)
print(beta)

def fs(alpha, beta):

  #alpha = fs/fu
  #beta  = f_lambda / (fu + fd)

  fs = 1.0 / (2/alpha + 1 + 2*beta / alpha)
  return fs

#the binning is the same
bins_fs = bins_cms
values_fs = []
unc_fs    = []
                      #fs/fu       #f_lambda / (fu + fd)
for a,b in zip (alpha, beta):

  values_fs.append(fs(a,b).nominal_value)
  unc_fs   .append(fs(a,b).std_dev      )


print("f_s values: \n", values_fs)

# PLOTTING

# get the bin centers
bin_center = []

for i in range(len(values_fs)):
  bin_center.append((bins_fs[i] + bins_fs[i+1])/2)

plt.figure()
plt.errorbar(bin_center, values_fs, yerr = unc_fs, marker = "v", color = "r", linestyle="none", ecolor="k",capsize=3)
plt.xlabel(r"$p_{T}$")
plt.ylabel(r"$f_{s}$")
plt.grid()
plt.savefig("fs.pdf")
