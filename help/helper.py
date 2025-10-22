import json

################ 
# 2024 samples #
################ 

## Flat ntuples
# unconstrained fit
sig_unc        = ["26_07_2024_14_44_54"]
hb_unc         = ["25_07_2024_14_23_01"]
bs_unc         = ["30_07_2024_11_40_54"]
b0_unc         = ["30_07_2024_11_40_34"]
bplus_unc      = ["30_07_2024_11_40_45"]

                  #part 1            #part 2            #part 3           #part 4             #part5
data_unc       = ["20240802_111807", "20240730_223445", "20240728_200528","20240729_090628" , "20240806_090127"]

# constrained fit
sig_cons_24    = ["26_07_2024_14_46_03"]
hb_cons_24     = ["25_07_2024_14_23_42"]
bs_cons_24     = ["31_07_2024_10_07_17"]
b0_cons_24     = ["31_07_2024_10_07_43"]
bplus_cons_24  = ["31_07_2024_10_08_00"]

                  #part 1            #part 2            #part 3            #part 4            #part 5
data_cons_24   = ["20240724_170443", "20240804_220622", "20240809_082548", "20240809_235436", "20240811_203518"]
bdt_data_24    = "25_04_2025_16_43_51" # with skimmed samples for all parts

fakes          = ["19_09_2024_17_38_06","19_09_2024_19_05_01","19_09_2024_21_24_27","20240919_184219"]

# NN model
code_24        = "26Sep2024_07h46m21s_cons" # used for inaugural talk

sig_cons_pastNN_24     = "sig_"    +code_24 
hb_cons_pastNN_24      = "hb_"     +code_24
bs_cons_pastNN_24      = "bs_"     +code_24
b0_cons_pastNN_24      = "b0_"     +code_24
bplus_cons_pastNN_24   = "bplus_"  +code_24
data_cons_pastNN_24    = "data_"   +code_24

#baseline selection
base_wout_tv_24 = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 8)', 
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'(fv_prob > 0.1)'
])
#hammered signals
sig_cons_hammer_24  = "signal_default_29_04_2025_13_51_27" #all events and all branches

################ 
# 2025 samples #
################ 

## Flat ntuples
#sig_cons_25    = ["16_03_2025_20_56_34"]
#hb_cons_25     = ["29_08_2025_18_35_52"] #["17_03_2025_08_36_11"]
#bs_cons_25     = ["17_03_2025_08_37_23"]
#b0_cons_25     = ["17_03_2025_08_37_03"]
#bplus_cons_25  = ["17_03_2025_08_36_31"]

                  #part 1            #part 2            #part 3            #part 4            #part 5
#data_cons_25     = ["20250227_155416", "20250227_161007", "20250227_161505", "20250227_161842", "20250227_161914"]

#UNCONSTRAINED FIT     
data_cons_25     = ["20250915_094746", "20250915_131727", "20250915_131849", "20250915_131952", "20250915_132020"]
sig_cons_25      = ["23_09_2025_13_17_17"]
hb_cons_25       = ["23_09_2025_13_12_55"] #(wout pileup info): with pileup: 30_09_2025_09_15_30

# BDT (not even needed for unconstrained fit)

#bdt models, one separate for each trigger
#bdt_model_25_mu7 = "21_07_2025_16_09_20" # this is in the Analysis Note 
#bdt_model_25_mu9 = "21_07_2025_17_14_06" # trained on mu9

#this is with both OR sign data trained
bdt_model_25_mu7="15_09_2025_20_27_30"
bdt_model_25_mu9="15_09_2025_20_58_00"
#bdt_model_25_mu7="16_09_2025_16_10_56" 
#bdt_model_25_mu9="16_09_2025_22_13_41"

#data file evaluated on the different trigger bdts
#bdt_data_25     = "21_07_2025_17_56_45"
#rhis is with both OR sign data trained
bdt_data_25     = "16_09_2025_21_41_00"



#signflip fit for relative normalization used in analysis note: 22_07_2025_07_52_30

#stacked plots before bdt weights used in analysis note: 22_07_2025_10_11_46 
#stacked plots after bdt weights used in analysis note: 22_07_2025_08_54_02 (same as below) 

#closure plots after bdt weights for high mass region in analysis note: 23_07_2025_09_16_32 
#closure plots after bdt weights for sidebands in analysis note: 23_07_2025_08_47_54 


#ds peak on dsmu signal mc plot used in analysis note: 06_06_2025_12_05_54
#shapes plot of different reconstructions: used in analysis note: 22_07_2025_07_52_30  (including BDT weights)
#more shapes to put in appendix if needed: used in analysis note: 22_07_2025_08_54_02  (including BDT weights)
#hb inclusive cocktail list is taken from this folder: inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1 in the /gen repo

# NN model used in analysis note
#nn_25_mu7 = "31Aug2025_20h37m02s"
#this is with both OR sign data trained
#nn_25_mu7 = "16Sep2025_08h40m52s"
#this is with unconstrained data but including phipi mass
#nn_25_mu7 = "23Sep2025_00h05m05s" 
#this is with unconstrained data excluding phipi mass from training
#nn_25_mu7 = "23Sep2025_22h35m11s"

#this is with new hammer weights
nn_25_mu7 = "20Oct2025_14h53m39s" #unconstrained data
nn_25_mu9 = "31Aug2025_22h05m54s"


#nn_25_mu7 = "22Oct2025_06h22m43s" # constrained data

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
f'(mu_pt > 8)', 
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(rel_iso_03_pv < 0.3)',
'(fv_prob > 0.1)',
'(mu_is_global == 1)',
'(ds_vtx_cosine_xyz_pv > 0.8)',
])

#hammered signals
#sig_cons_hammer_25 = "signal_default_10_06_2025_08_05_53"
#hammered signals with new harrison weights
#sig_cons_hammer_25 = "signal_default_01_10_2025_06_55_14"
# new hammered signals with unconstrained fit
#sig_cons_hammer_25 ="signal_default_23_09_2025_17_04_03"
# new hammered signals with unconstrained fit and new harrison weights
#sig_cons_hammer_25= "signal_default_30_09_2025_22_38_54"
#if i remove the factor of 2 everywhere..
#sig_cons_hammer_25 = "signal_default_17_10_2025_16_16_23" #on unconstraiend signal
sig_cons_hammer_25 = "signal_default_21_10_2025_21_39_09" #on constrained signal


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
# yaml file used in the analysis note: /work/pahwagne/RDsTools/hammercpp/development_branch/weights/04_06_2025_16_55_31
averageWeightsYaml = "01_10_2025_09_01_19" # the old (nauary 25 harrison weights): "04_06_2025_16_55_31"
# plots which show central and variational weight effect used in analysis note: 
# /work/pahwagne/RDsTools/hammercpp/development_branch/weights/plots/10_06_2025_13_04_57/



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
"base_wout_tv_24": base_wout_tv_24,
"base_wout_tv_25": base_wout_tv_25
}

# need this when we want to import (some)
# variables via .sh files (print and then read from shell)

if __name__ == "__main__":

  config = {
    "sig_cons_24": sig_cons_24,
    "hb_cons_24": hb_cons_24,
    "bs_cons_24": bs_cons_24,
    "b0_cons_24": b0_cons_24,
    "bplus_cons_24": bplus_cons_24,
    "data_cons_24": data_cons_24,
    "fakes": fakes,
    "sig_cons_25": sig_cons_25,
    "hb_cons_25": hb_cons_25,
    "bs_cons_25": bs_cons_25,
    "b0_cons_25": b0_cons_25,
    "bplus_cons_25": bplus_cons_25,
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

