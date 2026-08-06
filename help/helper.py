import json


#data_cons_25     = ["20250915_094746", "20250915_131727", "20250915_131849", "20250915_131952", "20250915_132020"]


#production with all triggers and pt cut 0.7
data_flatNanos=[
# BPH 1 total skimmed: 22087190 (with iso)
"22_05_2026_11_24_52", # part A skimmed: 1174665 (with iso)
"24_05_2026_17_46_53", # part B skimmed: 2388069 (with iso) 
"25_05_2026_11_46_22", # part C skimmed: 2344270 (with iso) 
"20260227_105644"    , # part D skimmed: 16180186 (with iso)

# BPH 2 total skimmed: 22370503 (with iso)
"27_05_2026_17_18_28", # part A skimmed: 1875601 (with iso)
"27_05_2026_22_31_02", # part B skimmed: 2389219 (with iso)  
"29_05_2026_07_05_38", # part C skimmed: 2324820 (with iso)  
"20260227_104615"    , # part D skimmed: 15780863 (with iso) 

# BPH 3 total skimmed: 22946220 (with iso)
"30_05_2026_09_11_47", # part A skimmed: 1887908 (with iso)
"31_05_2026_12_52_37", # part B skimmed: 2485389 (with iso)
"01_06_2026_02_11_46", # part C skimmed: 2334498 (with iso)
"20260227_102607"    ,  # part D skimmed: 16238425 (with iso)

# BPH 4 total skimmed 22929931 (with_iso):
"20260618_145647", # part A skimmed: 1933081 (with iso) 
"20260618_145823", # part B skimmed: 2512963 (with iso) 
"20260618_145947", # part C skimmed: 2426333 (with iso) 
"20260603_144926", # part D skimmed: 16057554 (with iso) 


# BPH 5 total skimmed: 21842060 (with_iso)
"02_06_2026_22_55_03", # part A skimmed: 1903193 (with_iso)
"03_06_2026_17_34_42", # part B skimmed: 2475207 (with_iso)
"04_06_2026_22_26_54", # part C skimmed: 1931850 (with_iso)
"20260409_141324"    , # part D skimmed: 15531810 (with_iso)

#BPH 6
"20260702_112403", # part A
"20260702_112800", # part B

] 

split_by_parts = {
"data_bph1" : data_flatNanos[ :4  ],
"data_bph2" : data_flatNanos[4:8  ],
"data_bph3" : data_flatNanos[8:12 ],
"data_bph4" : data_flatNanos[12:16],
"data_bph5" : data_flatNanos[16:20],
"data_bph6" : data_flatNanos[20:22],
}

#split_by_parts = {
#"data_bph1" :[ data_flatNanos[ 3  ]],
#"data_bph2" :[ data_flatNanos[ 7  ]],
#"data_bph3" :[ data_flatNanos[ 11 ]],
#"data_bph4" :[ data_flatNanos[ 15 ]],
#"data_bph5" :[ data_flatNanos[ 19 ]],
#}


#sig_flatNanos     = ["17_03_2026_11_31_56"] 
#sig_flatNanos     = ["check_sf_full/26_05_2026_17_56_10_fullBPark"] #including 8p5 and 10p5 
sig_flatNanos     = ["26_05_2026_17_56_10"] #including 8p5 and 10p5 
#hb_flatNanos      = ["25_02_2026_16_36_40"]  
#hb_flatNanos      = ["minimal/25_02_2026_16_36_40_D_only"]  
hb_flatNanos      = ["26_05_2026_17_59_33"] #including 8p5 and 10p5   

#### BDT 

#bdt_model = "27_02_2026_08_25_42" #model for 02Feb2026_13h29m45s
#bdt_model="27_02_2026_08_43_51" #evaluated data for model 02Feb2026_13h29m45s (back then different filename)

##bdt_data_afternn = "09_03_2026_10_22_26" #for model 02Feb2026_13h29m45s
#bdt_data_afternn = "24_02_2026_09_46_20" #for model 02Feb2026_09h24m08s 

#bdt_model="27_02_2026_08_43_51" #for model 02Feb2026_09h24m08s  


#bdt_model = "04_06_2026_09_34_56" #correction for 30 epochs NNmodel (2026_06_04_08_08_42), big bdt, bkg closes!
#bdt_model = "04_06_2026_15_02_06" #correction for 100 epochs NNmodel (2026_06_04_10_48_38)
#bdt_model = "09_06_2026_16_34_03" #correction for 100 epochs NNmodel (2026_06_04_10_48_38), but only trained in region score5 < 0.1
#bdt_model = "04_06_2026_16_39_36" #bdt trained with score5<0.9 for NNmodel (2026_06_04_10_48_38)
#bdt_model = "05_06_2026_07_50_44" #for bdt classifier
#bdt_model = "10_06_2026_13_33_07" #pre NN correction
#bdt_model = "11_06_2026_07_13_52" #pre NN correction, trained on high mass region only!
#bdt_model = "11_06_2026_10_58_47" #pre NN, trained in abs(cosPiK1) < 0.2, without training on mass
#bdt_model = "11_06_2026_16_28_21" #this is for the nn wout mass training
#bdt_model = "16_06_2026_08_43_47" # for nn 2026_06_15_10_03_36
#bdt_model = "18_06_2026_09_59_55" # for bdt 2026_06_15_10_23_15
#bdt_model = "23_06_2026_09_19_47" #for bdt 2026_06_19_09_32_52
#bdt_model="01_07_2026_09_22_13" #for bdt 2026_07_01_00_48_47
#bdt_model = "02_07_2026_07_53_24" #for bdt 2026_07_01_16_29_30

bdt_model="20_07_2026_22_07_48" ##for bdt 2026_07_01_00_48_47

####### REDO OLD LOGIC

#bdt_model = "16_07_2026_08_28_03" # this is a bdt model BEFORE NN training on base_wout and mu7 only


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
#nn_model = "02Feb2026_13h29m45s"

#all triggers!
#nn_model = "2026_06_04_08_08_42" #small NN, 30 epochs only, closes together with BDT above
#nn_model = "2026_06_04_10_48_38" #big NN model
#nn_model = "2026_05_29_16_52_19" #this is actually a bdt!
#nn_model = "2026_06_10_14_37_21" #trained on 30 epochs with sf weights
#nn_model = "2026_06_11_13_49_23" #trained on 30 epochs, without sf_weights, and without mass!!
#nn_model = "2026_06_12_07_12_31" #trained on 30 epochs, without sf_weights, even less kin. information!!
#nn_model = "2026_06_15_10_03_36" #big NN model, wout mass in training
#nn_model="2026_06_15_10_23_15" #big BDT model wout mass in training
#nn_model = "2026_06_19_09_32_52" #even bigger bdt :))

# HAMMER ISGW2  everywhere
nn_model = "2026_07_01_00_48_47" # bdt with photon info :))
#nn_model = "2026_07_01_16_29_30" # bdt wout photon info

#nn_model = "2026_07_14_08_42_58" #nn wirth photon info and mass
#nn_model = "2026_07_14_13_44_52" #same as for cms week presi (same features, mu7 only, no SF, but no bdt corr)

#nn_model = "2026_07_16_09_02_09"

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


sig_pastNN  = "sig_"    + nn_model
hb_pastNN   = "hb_"     + nn_model
#data_pastNN = "data_"   + nn_model


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

check_sf = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 7.0)', 
'(mu_pt < 20.0)',
'(abs(dxy_mu_sig_pv) > 3.0)',
'(abs(dxy_mu_sig_pv) < 15.0)',
'(k1_pt > 1.0)',
'(k2_pt > 1.0)',
'(pi_pt > 1.0)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(rel_iso_03_pv < 0.3)',
'(fv_prob > 0.1)',
'(mu_is_global == 1)',
'(ds_vtx_cosine_xyz_pv > 0.8)',
])

check_sf_D_only= ' && '.join([ #remove the charge and ds+mu mass cuts!
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

check_sf_full = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 7.0)', 
'(mu_pt < 50)', #match with SF binning 
'(abs(dxy_mu_sig_pv) > 3)',
'(abs(dxy_mu_sig_pv) < 100)',
'(k1_pt > 1.0)',
'(k2_pt > 1.0)',
'(pi_pt > 1.0)',
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

with_iso = ' && '.join([ #remove the charge and ds+mu mass cuts!
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
#sig_hammer_flatNano = "signal_default_20_03_2026_09_37_22" #same as line above but with minimal selection :D! 
#sig_hammer_flatNano = "signal_default_22_05_2026_13_20_18" #same as line above but with minimal selection  + lifetime uncertainty:D! 

#sig_hammer_flatNano = "signal_default_28_05_2026_12_38_44" #same as line above but with minimal selection  + lifetime uncertainty + riccardos SF - D only!!!! 
#sig_hammer_flatNano= "signal_default_28_05_2026_15_04_16" # UNFILTERED AND HAMMERED WITH CLN
#sig_hammer_flatNano= "signal_default_29_05_2026_07_47_53" # FILTERED: with_iso AND HAMMERED 
#sig_hammer_flatNano= "signal_default_30_06_2026_12_06_43" # UNFILTERED AND HAMMERED WITH ISGW2 
#sig_hammer_flatNano= "signal_default_01_07_2026_00_43_46" # FILTERED: with_iso AND HAMMERED 
sig_hammer_flatNano= "signal_default_30_07_2026_11_01_16" # FILTERED: with_iso AND HAMMERED AND PU WEIGHTS

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
dsstarmu_isgw2_to_bgl = "dsstarmu_isgw2_BGLVar_02_07_2026_15_11_07" # hammer with factor 2: "dsstarmu_BGLVar_30_09_2025_22_37_38" #the old (january 25 harrison weights): "dsstarmu_BGLVar_04_06_2025_15_03_24"
dstau_to_bcl      = "dstau_BCLVar_04_06_2025_14_58_02"
dsstartau_to_bgl  = "dsstartau_BGLVar_20_10_2025_09_37_44" # hammer with factor 2: "dsstartau_BGLVar_30_09_2025_22_37_46" # the old (january 25 harrison weights): "dsstartau_BGLVar_04_06_2025_15_25_47"
# yaml file used in the analysis note: /work/pahwagne/RDsTools/hammercpp/development_branch/weights/20_10_2025_09_57_12/
#averageWeightsYaml = "20_10_2025_09_57_12" #with the CLN MC #the old (wrong factor 2 weights) "01_10_2025_09_01_19" # the old (nauary 25 harrison weights): "04_06_2025_16_55_31"
averageWeightsYaml = "02_07_2026_15_35_35" #with the new all-ISGW2 MC
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
"with_iso"         : with_iso,
"check_sf"       : check_sf,
"check_sf_full"       : check_sf_full,
"check_sf_D_only"       : check_sf_D_only,
}

# need this when we want to import (some)
# variables via .sh files (print and then read from shell)

if __name__ == "__main__":

  config = {
    "sig_flatNanos"   : sig_flatNanos,
    "hb_flatNanos"    : hb_flatNanos,
    "data_flatNanos"  : data_flatNanos,
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
    "dsstarmu_isgw2_to_bgl": dsstarmu_isgw2_to_bgl,
    "dstau_to_bcl": dstau_to_bcl,
    "dsstartau_to_bgl": dsstartau_to_bgl,
    "averageWeightsYaml": averageWeightsYaml,
    "sig_hammer_flatNano":sig_hammer_flatNano,
    "base_wout_tv_25": base_wout_tv_25,
    "with_iso":with_iso 
    }
  print(json.dumps(config))

