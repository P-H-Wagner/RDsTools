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
sig_cons_25    = ["16_03_2025_20_56_34"]
hb_cons_25     = ["17_03_2025_08_36_11"]
bs_cons_25     = ["17_03_2025_08_37_23"]
b0_cons_25     = ["17_03_2025_08_37_03"]
bplus_cons_25  = ["17_03_2025_08_36_31"]

                  #part 1            #part 2            #part 3            #part 4            #part 5
data_cons_25   = ["20250227_155416", "20250227_161007", "20250227_161505", "20250227_161842", "20250227_161914"]
bdt_data_25    = "29_05_2025_16_27_10" # this is in the Analysis Note 

# NN model
code_25 = "28May2025_16h33m29s"

sig_cons_pastNN_25     = "sig_"    +code_25
hb_cons_pastNN_25      = "hb_"     +code_25
bs_cons_pastNN_25      = "bs_"     +code_25
b0_cons_pastNN_25      = "b0_"     +code_25
bplus_cons_pastNN_25   = "bplus_"  +code_25
data_cons_pastNN_25    = "data_"   +code_25

#baseline selection
base_wout_tv_25 = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 8)', 
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(rel_iso_03_pv < 0.3)',
'(fv_prob > 0.1)'
])

#hammered signals
sig_cons_hammer_25 = "signal_default_03_06_2025_17_06_15"

################
# Hammer tools #
################

# unfiltered gen-level MC samples used for hammer circular test
# and calculation of average weights

dsmu_gen       = "15_11_2024_06_57_58"
dsmu_isgw2_gen = "17_11_2024_17_26_37"
dsstarmu_gen   = "15_11_2024_09_43_21"
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
# plots for circular test used in the analysis note: /work/pahwagne/RDsTools/hammercpp/tests/04_06_2025_16_25_50
circularTestPlots = "04_06_2025_16_25_50"

# unfiltered gen-level weighted files used to calculate the average weight for BCL/BGL
dsmu_to_bcl       = "dsmu_BCLVar_04_06_2025_14_40_48"
dsstarmu_to_bgl   = "dsstarmu_BGLVar_04_06_2025_15_03_24"
dstau_to_bcl      = "dstau_BCLVar_04_06_2025_14_58_02"
dsstartau_to_bgl  = "dsstartau_BGLVar_04_06_2025_15_25_47"
# yaml file used in the analysis note: /work/pahwagne/RDsTools/hammercpp/development_branch/weights/04_06_2025_16_55_31
averageWeightsYaml = "04_06_2025_16_55_31"

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
    "dstau_gen"     : dstau_gen,       
    "dsstartau_gen" : dsstartau_gen,  
    "dsmu_isgw2_to_cln": dsmu_isgw2_to_cln,
    "dsmu_to_isgw2": dsmu_to_isgw2,
    "dsmu_to_bcl": dsmu_to_bcl, 
    "dsstarmu_to_bgl": dsstarmu_to_bgl,
    "dstau_to_bcl": dstau_to_bcl,
    "dsstartau_to_bgl": dsstartau_to_bgl,
    "averageWeightsYaml": averageWeightsYaml,
    "sig_cons_hammer_25":sig_cons_hammer_25,
    "base_wout_tv_25": base_wout_tv_25
    }
  print(json.dumps(config))

