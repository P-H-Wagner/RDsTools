import json


################ 
# 2025 samples #
################ 


#### FULLY CONSTRAINED FIT     

#sig_cons_25    = ["16_03_2025_20_56_34"]
#hb_cons_25     = ["29_08_2025_18_35_52"] #["17_03_2025_08_36_11"]
#bs_cons_25     = ["17_03_2025_08_37_23"]
#b0_cons_25     = ["17_03_2025_08_37_03"]
#bplus_cons_25  = ["17_03_2025_08_36_31"]
                  #part 1            #part 2            #part 3            #part 4            #part 5
#data_cons_25     = ["20250227_155416", "20250227_161007", "20250227_161505", "20250227_161842", "20250227_161914"]

#### PARTIALLY UNCONSTRAINED FIT     

# kk constrained
data_cons_25     = ["20250915_094746", "20250915_131727", "20250915_131849", "20250915_131952", "20250915_132020"]
sig_cons_25      = ["23_09_2025_13_17_17"]
hb_cons_25       = ["26_01_2026_17_50_33"] #["23_09_2025_13_12_55"] #(wout pileup info): with pileup: 30_09_2025_09_15_30

#kk unconstrained
#data_cons_25 = ["20251027_161544", "20251027_161733", "20251027_161828", "20251027_161959", "20251027_162140"]
#sig_cons_25  = ["03_11_2025_16_01_07"]
#hb_cons_25   = ["03_11_2025_10_48_49"]

#### BDT 

# kk constrained
bdt_model_25_mu7= "28_01_2026_09_11_15" #"16_01_2026_12_31_25"  
bdt_model_25_mu9= "28_01_2026_09_49_54" #"16_01_2026_11_30_13" 

# kk unconstrained
#bdt_model_25_mu7= "26_11_2025_13_40_02" #"21_11_2025_15_49_43" #less features: "11_11_2025_13_41_34"
#bdt_model_25_mu9= "26_11_2025_16_17_08" #less features: "11_11_2025_13_41_37"

# kk constrained 
#bdt_data_25     = "26_11_2025_16_36_16" #"21_11_2025_16_51_16" #less features: "12_11_2025_13_51_29" #"16_09_2025_21_41_00"
bdt_data_25     = "16_01_2026_11_49_11" #"16_01_2026_12_54_24" #trained in highmass and low mass: "16_01_2026_11_49_11" 

# kk unconstrained
#bdt_data_25     = "26_11_2025_16_33_07" #"21_11_2025_16_32_00" #less features: missing?  #didnt cut at track pt 1 here: "12_11_2025_13_57_45"

#### signflip fit for relative normalization used in analysis note     : 08_12_2025_13_34_12 

# stacked plots before bdt weights used in analysis note               : 22_07_2025_10_11_46 
# stacked plots after bdt weights used in analysis note                : 22_07_2025_08_54_02 (same as below) 

# closure plots after bdt weights for high mass region in analysis note: 09_12_2025_00_09_43 #mu7  
# closure plots after bdt weights for high mass region in analysis note: 08_12_2025_18_41_16 #mu9 
# closure plots after bdt weights for sidebands in analysis note       : 09_12_2025_00_10_19 #mu7
# closure plots after bdt weights for sidebands in analysis note       : 08_12_2025_18_19_03 #mu9

# ds peak on dsmu signal mc plot used in analysis note                 : 08_12_2025_10_05_28 #on mu7 
# shapes plot of different reconstructions: used in analysis note      : 08_12_2025_10_01_51 #on mu7
# shapes plot of different reconstructions: not yet in analysis note   : #on mu9
# more shapes to put in appendix if needed: used in analysis note      : 
# hb inclusive cocktail list is taken from this folder                 : inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1 in the /gen repo

# kk constrained 
nn_25_mu7 = "05Dec2025_11h14m48s" #"26Nov2025_16h40m23s" #"25Nov2025_16h14m45s" #"12Nov2025_14h02m26s" # find the norm score plot: test_12Nov2025_16h56m34s resp to 0.4 only: test_12Nov2025_17h14m00s
nn_25_mu9 = "05Dec2025_11h15m22s"

# kk unconstrained 
#nn_25_mu7 = "01Dec2025_15h16m00s" #"26Nov2025_19h21m09s" #"21Nov2025_16h56m09s"#"19Nov2025_09h20m20s" # this is with no track pt cut: "12Nov2025_09h57m07s" # 500 epochs but wrong hammer: "11Nov2025_16h20m49s" find the norm score plot: test_12Nov2025_16h50m09s resp to 0.4 only: test_12Nov2025_17h15m06s 
#nn_25_mu9 = "11Nov2025_15h57m14s"

cons_pastNN_25_mu7    = {}
cons_pastNN_25_mu7["sig"  ] = "sig_"    + nn_25_mu7
cons_pastNN_25_mu7["hb"   ] = "hb_"     + nn_25_mu7
cons_pastNN_25_mu7["bs"   ] = "bs_"     + nn_25_mu7
cons_pastNN_25_mu7["b0"   ] = "b0_"     + nn_25_mu7
cons_pastNN_25_mu7["bplus"] = "bplus_"  + nn_25_mu7
cons_pastNN_25_mu7["data" ] = "data_"   + nn_25_mu7

cons_pastNN_25_mu9    = {}
cons_pastNN_25_mu9["sig"  ] = "sig_"    + nn_25_mu9
cons_pastNN_25_mu9["hb"   ] = "hb_"     + nn_25_mu9
cons_pastNN_25_mu9["bs"   ] = "bs_"     + nn_25_mu9
cons_pastNN_25_mu9["b0"   ] = "b0_"     + nn_25_mu9
cons_pastNN_25_mu9["bplus"] = "bplus_"  + nn_25_mu9
cons_pastNN_25_mu9["data" ] = "data_"   + nn_25_mu9

#baseline selection
base_wout_tv_25 = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 8.0)', 
'(k1_pt > 1.0)',
'(k2_pt > 1.0)',
'(pi_pt > 1.0)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(rel_iso_03_pv < 0.3)',
#'(kk_m < 1.025) ',
#'(kk_m > 1.015) ',
'(fv_prob > 0.1)',
'(mu_is_global == 1)',
'(ds_vtx_cosine_xyz_pv > 0.8)',
])

# kk constrained 
sig_cons_hammer_25 = "signal_default_17_10_2025_16_16_23" 
# kk unconstrained
#sig_cons_hammer_25 = "signal_default_04_11_2025_17_41_46"

isoflip = ' && '.join([ #remove the charge and ds+mu mass cuts! 
'(rel_iso_03_pv > 0.3)',
'(ds_vtx_cosine_xyz_pv < 0.8)',
])

################
# Hammer tools #
################

# unfiltered gen-level MC samples used for hammer circular test
# and calculation of average weights

dsmu_gen       = "15_11_2024_06_57_58"
dsmu_isgw2_gen = "17_11_2024_17_26_37"
dsstarmu_gen   = "15_11_2024_09_43_21"
#dsstarmu_isgw2_gen = "19_11_2024_08_23_27"
dsstarmu_isgw2_gen = "07_10_2025_08_03_59" #with explicit neutrino tagging
dstau_gen      = "15_11_2024_13_50_26"
dsstartau_gen  = "16_11_2024_09_45_34"


#systematic unc. list
scalar_model = "Bcl"
systematics_scalar = [
"e1",
"e2",
"e3",
"e4",
"e5",
"e6"
]
vector_model = "Bgl"
systematics_vector = [
"e1",
"e2",
"e3",
"e4",
"e5",
"e6",
"e7",
"e8",
"e9",
"e10"
]

# unfiltered gen-level weighted files of circular test

# dsmu with isgw2 re-weighted to HQET2 (CLN)
dsmu_isgw2_to_cln = "dsmu_isgw2_CLN_04_06_2025_11_29_41"
#dsmu with HQET2 (CLN) to ISGW2
dsmu_to_isgw2     = "dsmu_ISGW2_04_06_2025_11_34_00"
# plots for circular test used in the analysis note: /work/pahwagne/RDsTools/hammercpp/tests/22_07_2025_08_25_22
circularTestPlots = "22_07_2025_08_25_22"

# unfiltered gen-level weighted files used to calculate the average weight for BCL/BGL
dsmu_to_bcl       = "dsmu_BCLVar_04_06_2025_14_40_48"
dsmu_isgw2_to_bcl = "dsmu_isgw2_BCLVar_26_06_2025_12_00_43"
dsstarmu_to_bgl   = "dsstarmu_BGLVar_20_10_2025_09_48_08" # hammer with factor 2: "dsstarmu_BGLVar_30_09_2025_22_37_38" #the old (january 25 harrison weights): "dsstarmu_BGLVar_04_06_2025_15_03_24"
dstau_to_bcl      = "dstau_BCLVar_04_06_2025_14_58_02"
dsstartau_to_bgl  = "dsstartau_BGLVar_20_10_2025_09_37_44" # hammer with factor 2: "dsstartau_BGLVar_30_09_2025_22_37_46" # the old (january 25 harrison weights): "dsstartau_BGLVar_04_06_2025_15_25_47"
# yaml file used in the analysis note: /work/pahwagne/RDsTools/hammercpp/development_branch/weights/20_10_2025_09_57_12/
averageWeightsYaml = "20_10_2025_09_57_12" #the old (wrong factor 2 weights) "01_10_2025_09_01_19" # the old (nauary 25 harrison weights): "04_06_2025_16_55_31"
# plots which show central and variational weight effect used in analysis note: 
# /work/pahwagne/RDsTools/hammercpp/development_branch/weights/plots/09_12_2025_09_01_30/ #mu7
# /work/pahwagne/RDsTools/hammercpp/development_branch/weights/plots/09_12_2025_09_16_52/ #mu9


# prefit plots
# 09_12_2025_11_37_43 #mu7
# 09_12_2025_11_37_40 #mu9

#shapes to fit:
# 09_12_2025_18_31_30 #mu7
# 09_12_2025_22_53_09 #mu9

# Mass constants

dsMass_       = 1.96834
bsMass_       = 5.36688
phiMass_      = 1.019461

# Define sideband region

nSignalRegion = 3 # signal region is 3 sigma
nSidebands    = 4 # sideband region starts after 5 sigma
sbWidth       = 2 # sideband region is 1 sigma broad


#fill all relevant into dictionary
baselines = {
"base_wout_tv_25": base_wout_tv_25
}

# need this when we want to import (some)
# variables via .sh files (print and then read from shell)

if __name__ == "__main__":

  config = {
    "sig_cons_25": sig_cons_25,
    "hb_cons_25": hb_cons_25,
    "data_cons_25": data_cons_25,
    "dsmu_gen"      : dsmu_gen,           
    "dsmu_isgw2_gen": dsmu_isgw2_gen, 
    "dsstarmu_gen"  : dsstarmu_gen,   
    "dsstarmu_isgw2_gen"  : dsstarmu_isgw2_gen,   
    "dstau_gen"     : dstau_gen,       
    "dsstartau_gen" : dsstartau_gen,  
    "dsmu_isgw2_to_cln": dsmu_isgw2_to_cln,
    "dsmu_to_isgw2": dsmu_to_isgw2,
    "dsmu_to_bcl": dsmu_to_bcl, 
    "dsmu_isgw2_to_bcl": dsmu_isgw2_to_bcl, 
    "dsstarmu_to_bgl": dsstarmu_to_bgl,
    "dstau_to_bcl": dstau_to_bcl,
    "dsstartau_to_bgl": dsstartau_to_bgl,
    "averageWeightsYaml": averageWeightsYaml,
    "sig_cons_hammer_25":sig_cons_hammer_25,
    "base_wout_tv_25": base_wout_tv_25
    }
  print(json.dumps(config))

