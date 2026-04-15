import json


#data_cons_25     = ["20250915_094746", "20250915_131727", "20250915_131849", "20250915_131952", "20250915_132020"]


#production with all triggers and pt cut 0.7


data_cons_25     = [
# BPH 1,2,3 of part D
"20260227_105644", 
"20260227_104615", 
"20260227_102607"
] 

sig_cons_25      = ["17_03_2026_11_31_56"] #["23_09_2025_13_17_17"]
hb_cons_25       = ["25_02_2026_16_36_40"]  

#### BDT 
bdt_model_25_mu7        = "27_02_2026_08_25_42" #"19_02_2026_13_23_30" #<-- cont. weights #"28_01_2026_09_11_15"   
bdt_model_afternn_25_mu7= "09_03_2026_10_14_51"   
#bdt_model_afternn_25_mu7= "12_03_2026_16_46_06"   

bdt_model_25_mu9        = "27_02_2026_08_40_28" #"19_02_2026_13_28_24" #<-- cont. weights # "28_01_2026_09_49_54"   
bdt_model_afternn_25_mu9= "24_02_2026_09_31_30"   
#bdt_model_afternn_25_mu9= "12_03_2026_16_50_12"   

bdt_model_afternn = "15_04_2026_09_30_45" #for minimal model: 2026_04_13_12_15_26


bdt_model = "23_03_2026_15_05_14" #all triggers, minimal selection
bdt_data = "23_03_2026_15_30_57" #all triggers, minimal selection

bdt_data_25     = "27_02_2026_08_43_51"#"26_02_2026_11_52_09"#"19_02_2026_14_01_51" # <--- cont. weight "30_01_2026_09_38_39"    

bdt_data_afternn_mu7 ="09_03_2026_10_22_26"
#bdt_data_afternn_mu7 ="12_03_2026_16_53_36"
bdt_data_afternn_mu9 ="24_02_2026_09_46_20" 
#bdt_data_afternn_mu9 ="12_03_2026_16_55_15" 

bdt_data_afternn = "15_04_2026_10_42_07"


#### signflip fit for relative normalization used in analysis note     : 08_12_2025_13_34_12 

# stacked plots before bdt weights used in analysis note               : 22_07_2025_10_11_46 
# stacked plots after bdt weights used in analysis note in SR          : 29_01_2026_11_52_28, 29_01_2026_11_50_43, 29_01_2026_11_51_55, 29_01_2026_11_52_54, 29_01_2026_11_51_14 #mu7 
# stacked plots after bdt weights used in analysis note in SR          : 30_01_2026_13_45_24,30_01_2026_13_44_35, 30_01_2026_13_42_15,30_01_2026_13_43_21, 30_01_2026_13_44_57 #mu9 

# closure plots after bdt weights for high mass region in analysis note: 28_01_2026_10_42_04 #mu7  
# closure plots after bdt weights for high mass region in analysis note: 30_01_2026_09_55_20, 30_01_2026_09_53_31, 30_01_2026_09_54_30, 30_01_2026_09_56_08, 30_01_2026_10_24_26 #mu9 

# closure plots after bdt weights for left sideband in analysis note   : 28_01_2026_10_43_30 #mu7
# closure plots after bdt weights for left sideband in analysis note   : 30_01_2026_11_35_05,30_01_2026_11_35_49,30_01_2026_11_35_21, 30_01_2026_11_36_23, 30_01_2026_11_24_51 #mu9
# closure plots after bdt weights for right sideband in analysis note  : 29_01_2026_09_59_08, 29_01_2026_10_03_21, 29_01_2026_10_00_25, 29_01_2026_10_02_26, 29_01_2026_10_01_23 #mu7 (ch0)

# closure plots after bdt weights for right sideband in analysis note  : 30_01_2026_12_31_27,30_01_2026_12_31_02, 30_01_2026_12_33_24,30_01_2026_12_30_11, 30_01_2026_12_32_55 #mu9

# ds peak on dsmu signal mc plot used in analysis note                 : 08_12_2025_10_05_28 #on mu7 
# shapes plot of different reconstructions: used in analysis note      : 08_12_2025_10_01_51 #on mu7
# shapes plot of different reconstructions: not yet in analysis note   : #on mu9
# more shapes to put in appendix if needed: used in analysis note      : 
# hb inclusive cocktail list is taken from this folder                 : inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1 in the /gen repo

# kk constrained 
nn_25_mu7 = "02Feb2026_13h29m45s" #"05Dec2025_11h14m48s" #with new hb: 02Feb2026_13h29m45s
nn_25_mu9 = "02Feb2026_09h24m08s" #"05Dec2025_11h15m22s" #with new hb: 02Feb2026_09h24m08s

#nn_25_mu7 = "2026_04_10_15_16_27" #with the same trainer as for the all triggers but only on mu7
#nn_25_mu7 = "2026_04_10_15_22_39" #with the same trainer as for the all triggers but only on mu9

#all triggers!
#nn_model = "2026_04_08_15_32_42" #with sf weights on minimal
#nn_model = "2026_04_09_09_40_49" #without sf weights on minimal
#nn_model = "2026_04_10_10_18_23" #with sf weights on minimal2
#nn_model = "2026_04_10_09_05_51" #with sf weights on base_wout_tv
#nn_model  = "2026_04_10_13_31_41" #with sf weights on base_wout_tv
#nn_model = "test_10Apr2026_16h48m13s" #mu7 only with old script (before splitting)
#nn_model = "test_02Feb2026_13h29m45s"
#nn_model= "2026_04_13_09_46_22" #mu7 only and same arch as 02Feb model but with new script

##nn_model = "2026_04_10_15_16_27" #with the same trainer as for the all triggers but only on mu7
#nn_model = "2026_04_10_15_22_39" #with the same trainer as for the all triggers but only on mu9

#nn_model= "2026_04_13_12_09_32" #all triggers, odl arch, base_wout_tv_25
nn_model= "2026_04_13_12_15_26" #all triggers, odl arch, minimal 

#nn_25_mu7 = "09Mar2026_10h05m52s" #"05Dec2025_11h14m48s" #with new hb: 02Feb2026_13h29m45s
#nn_25_mu9 = "09Mar2026_15h07m00s" #"05Dec2025_11h15m22s" #with new hb: 02Feb2026_09h24m08s

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


pastNN_sig  = "sig_"    + nn_model
pastNN_hb   = "hb_"     + nn_model
pastNN_data = "data_"   + nn_model


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

offline = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 7.0)', 
'(k1_pt > 0.7)',
'(k2_pt > 0.7)',
'(pi_pt > 0.7)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(rel_iso_03_pv < 0.3)',
'(fv_prob > 0.1)',
'(mu_is_global == 1)',
'(ds_vtx_cosine_xyz_pv > 0.8)',
])

minimal = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 7.0)', 
'(k1_pt > 0.7)',
'(k2_pt > 0.7)',
'(pi_pt > 0.7)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
#'(rel_iso_03_pv < 0.3)',
#'(fv_prob > 0.1)',
'(mu_is_global == 1)',
'(ds_vtx_cosine_xyz_pv > 0.8)',
])

minimal2 = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 7.0)', 
'(k1_pt > 0.7)',
'(k2_pt > 0.7)',
'(pi_pt > 0.7)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(rel_iso_03_pv < 0.3)',
#'(fv_prob > 0.1)',
'(mu_is_global == 1)',
'(ds_vtx_cosine_xyz_pv > 0.8)',
])

# kk constrained 
#sig_cons_hammer_25 = "signal_default_17_10_2025_16_16_23" 
#sig_cons_hammer_25 = "signal_default_16_03_2026_13_34_54" #same as line above but with scale factors :D! 
#sig_cons_hammer_25 = "signal_default_17_03_2026_18_28_42" #same as line above but with scale factors :D! 
sig_cons_hammer_25 = "signal_default_20_03_2026_09_37_22" #same as line above but with minimal selection :D! 

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


#unfiltered gen lvl sample including pv,sv, beta etc. used to calc. Bs decay time
dsmu_for_time = "19_03_2026_12_19_50" #inside pnfs/.../flatNano
#unfiltered gen lvl sample including time 
dsmu_with_time = "19_03_2026_16_32_11"

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
"base_wout_tv_25": base_wout_tv_25,
"offline"        : offline, 
"minimal"        : minimal,
"minimal2"        : minimal2,
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

