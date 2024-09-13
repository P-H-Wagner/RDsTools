import ROOT
import numpy as np
import uproot
import glob
import sys

#use uproot to load file
samples = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/13_09_2024_15_49_45/*"
sample_list = glob.glob(samples)

#branches to load
branches = [
"gen_sig",
"q2_coll",

"mu_id_loose", 
"mu_id_medium",              
"mu_id_tight",               

"mu_id_medium_promt"       , 
"mu_id_global_high_pt"     , 
"mu_id_trk_high_pt"        , 
"mu_pf_iso_very_loose"     , 
"mu_pf_iso_loose"          , 
"mu_pf_iso_medium"         , 
"mu_pf_iso_tight"          , 
"mu_pf_iso_very_tight"     , 

"mu_tk_iso_loose"          , 
"mu_tk_iso_tight"          , 
"mu_id_soft"               , 
"mu_id_soft_mva"           , 
"mu_mva_loose"             , 
"mu_mva_medium"            , 
"mu_mvs_tight"             , 
"mu_mini_iso_loose"        , 
"mu_mini_iso_medium"       , 
"mu_mini_iso_tight"        , 

"mu_mini_iso_very_tight"   , 
"mu_trigger_id_loose"      , 
"mu_in_time_muon"          , 
"mu_pf_iso_very_very_tight",
"mu_multi_iso_loose"       , 
"mu_multi_iso_medium"      , 
"mu_puppi_iso_loose"       , 
"mu_puppi_iso_medium"      , 
"mu_puppi_iso_tight"       , 
"mu_mva_v_tight"           , 

"mu_mva_vv_tight"          , 
"mu_low_pt_mva_loose"      , 
"mu_low_pt_mva_medium"     , 
"mu_mv_id_wp_medium"       , 
"mu_mv_id_wp_tight"        , 

]

ids = branches[2:]

channels  = [
-10,
0,
1,
10,
11, 
-1
]

labels    = [
"fakes", 
"ds mu", 
"ds tau", 
"ds star mu", 
"ds star tau", 
"hb"
]

labels = {chan: labels[i] for i,chan in enumerate(channels)}

selection = [
"(gen_sig == -10)", 
"(gen_sig == 0)", 
"(gen_sig == 1)", 
"(gen_sig == 10)",
"(gen_sig == 11)",
"(gen_sig != -9999) & (gen_sig != 0) & (gen_sig != 1) & (gen_sig != 10) & (gen_sig != 11) & (gen_sig != -10) & (gen_match_success ==1)"
]

selection = {chan: selection[i] for i,chan in enumerate(channels)}

df = {}


for i,sample in enumerate(sample_list):

  with uproot.open(sample) as f:
    
    for chan in channels:

      tree = f["tree"]

      if i == 0:
        #starting value
        df[chan]   = tree.arrays(branches, library = "pd", entry_start=None, entry_stop=None, cut=selection[chan]) 

      else:
        #append
        df[chan] = df[chan].append(tree.arrays(branches, library = "pd", entry_start=None, entry_stop=None, cut=selection[chan]) ) 


#best is to compare dsmu against fakes
 
import matplotlib.pyplot as plt

merits = []
for sel in ids:

  S = len( df[10][ df[10][sel]  == 1 ])
  B = len( df[-10][df[-10][sel] == 1 ])

  print("Selection: ", sel)
  print("S: ", S)
  print("B: ", B)

  if ((S+B) != 0):
    merit = S / np.sqrt(S + B)
  else:
    merit = np.nan
  
  merits.append(merit)

#plotting
x = list(range(len(merits)))

plt.figure(figsize=(15, 15))
plt.barh(x,merits, align = 'center')
plt.xlabel("S / sqrt(S+B)")
plt.yticks(x,ids)
plt.gcf().subplots_adjust(left=0.2) #space on the left side
plt.show()
plt.savefig("merit.pdf")
plt.savefig("merit.png")
  






