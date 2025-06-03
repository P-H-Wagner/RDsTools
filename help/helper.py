import json

# samples

sig_unc     = ["26_07_2024_14_44_54"]
hb_unc      = ["25_07_2024_14_23_01"]
bs_unc      = ["30_07_2024_11_40_54"]
b0_unc      = ["30_07_2024_11_40_34"]
bplus_unc   = ["30_07_2024_11_40_45"]

               #part 1            #part 2            #part 3           #part 4             #part5
data_unc    = ["20240802_111807", "20240730_223445", "20240728_200528","20240729_090628" , "20240806_090127"]

sig_cons_24    = ["26_07_2024_14_46_03"]
hb_cons_24     = ["25_07_2024_14_23_42"]
bs_cons_24     = ["31_07_2024_10_07_17"]
b0_cons_24     = ["31_07_2024_10_07_43"]
bplus_cons_24  = ["31_07_2024_10_08_00"]
               #part 1            #part 2            #part 3            #part 4            #part 5
data_cons_24   = ["20240724_170443", "20240804_220622", "20240809_082548", "20240809_235436", "20240811_203518"]

fakes = ["19_09_2024_17_38_06","19_09_2024_19_05_01","19_09_2024_21_24_27","20240919_184219"]


# new production round february 25

sig_cons_25    = ["16_03_2025_20_56_34"]
hb_cons_25     = ["17_03_2025_08_36_11"]
bs_cons_25     = ["17_03_2025_08_37_23"]
b0_cons_25     = ["17_03_2025_08_37_03"]
bplus_cons_25  = ["17_03_2025_08_36_31"]
               #part 1            #part 2            #part 3            #part 4            #part 5
data_cons_25   = ["20250227_155416", "20250227_161007", "20250227_161505", "20250227_161842", "20250227_161914"]

bdt_data_24 = "25_04_2025_16_43_51" # with skimmed samples for all parts
#bdt_data_25 = "25_04_2025_10_29_02" # with skimmed samples for all parts
bdt_data_25 = "22_05_2025_18_01_26" # with skimmed samples for all parts

# aftter NN samples
code = "05Sep2024_15h09m02s"
code = "08Sep2024_10h21m10s" 
code = "07Sep2024_13h22m00s_unc" # lets try to evaluate it on an constrained NN model!
code = "12Sep2024_13h28m24s_unc"
code = "22Sep2024_15h26m18s_unc" #unconstraiend data with signflip method trained
code = "23Sep2024_18h11m55s_unc" #unconstrained data with 3-5 sigma sideband trained
code = "24Sep2024_10h41m22s_unc" #signflip method without mass - still bad modeling!
code = "26Sep2024_07h46m21s_unc" #signflip + sideband
#code = "02Oct2024_09h45m36s_unc" #only signflip
#code = "19Aug2024_19h52m02s"
sig_unc_pastNN      = "sig_"    +code
hb_unc_pastNN       = "hb_"     +code
bs_unc_pastNN       = "bs_"     +code
b0_unc_pastNN       = "b0_"     +code
bplus_unc_pastNN    = "bplus_"  +code
data_unc_pastNN     = "data_"   +code

code = "05Sep2024_13h36m33s"
code = "07Sep2024_13h22m00s" #corresponding model: test_07Sep2024_13h22m00s
code = "10Sep2024_10h57m30s_cons" #for vfs presi 
code = "12Sep2024_11h21m46s_cons" #50 epochs and swish
code = "12Sep2024_20h19m27s_cons" #100 epochs
code = "16Sep2024_08h01m39s_cons" #200 epochs
code = "22Sep2024_11h49m06s_cons" #corresponds to 22Sep2024_15h26m18s_unc (same model but constrained) 
code = "24Sep2024_10h41m22s_cons" #wout mass - still bad modeling!
code = "26Sep2024_07h46m21s_cons" #sb and sf #for inaugural --> what we usually use!!
#code = "21Mar2025_14h25m24s_cons" #using hammer weights for training
#code = "21Mar2025_18h04m03s_cons" # using again sb and sf
#code = "23Aug2024_19h41m42s" #6 classes from optuna
code = "23Apr2025_18h57m36s_cons" #new NN with pimu wrong only + weights 
code = "25Apr2025_17h05m56s_cons" #new NN with pimu wrong only + weights hammer
code_24 = "29Apr2025_22h58m14s_cons" #new NN with pimu wrong + weights + new hammer + old data 24
#code_25 = "20May2025_11h13m16s_cons" # 25 prod :)
#code_25 = "22May2025_09h31m46s_cons" #25 prod with all features but only 50 epochs
code_25 = "23May2025_01h04m20s_cons" # with new isolation cut
#code_25 = "23May2025_14h33m00s_cons" #no hammer

code_25 = "28May2025_16h33m29s"

sig_cons_pastNN_24     = "sig_"    +code_24 
hb_cons_pastNN_24      = "hb_"     +code_24
bs_cons_pastNN_24      = "bs_"     +code_24
b0_cons_pastNN_24      = "b0_"     +code_24
bplus_cons_pastNN_24   = "bplus_"  +code_24
data_cons_pastNN_24    = "data_"   +code_24


sig_cons_pastNN_25     = "sig_"    +code_25
hb_cons_pastNN_25      = "hb_"     +code_25
bs_cons_pastNN_25      = "bs_"     +code_25
b0_cons_pastNN_25      = "b0_"     +code_25
bplus_cons_pastNN_25   = "bplus_"  +code_25
data_cons_pastNN_25    = "data_"   +code_25


#after hammer samples
#2024 prdouction
sig_cons_hammer     = "signal_BGLVar_13_01_2025_20_38_04"
sig_cons_hammer     = "signal_BGLVar_16_01_2025_13_40_56"
sig_cons_hammer     = "signal_default_10_02_2025_10_07_42" #includes weights for all signals, only part of signal
sig_cons_hammer     = "signal_default_07_03_2025_13_33_38" #all events
sig_cons_hammer     = "signal_default_29_04_2025_13_51_27" #all events and all branches
sig_unc_hammer = ""

dsStarTau_w = 3.304854088595039
dsStarMu_w  = 3.198165968764498

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


#gen production samples (without filter, for hammer average)
dstau_gen     = "14_10_2024_16_31_39" 
dsstartau_gen = "14_10_2024_16_33_54"

# Mass constants

dsMass_       = 1.96834
bsMass_       = 5.36688
phiMass_      = 1.019461

# Define sideband region

nSignalRegion = 3 # signal region is 3 sigma
nSidebands    = 4 # sideband region starts after 5 sigma
sbWidth       = 2 # sideband region is 1 sigma broad

# alpha = efficiency * f_i * BR 
# where f_i is the fragmentation for quark type i = u,s,d, ...
#               #epsilon    #BR                                             #the old ones without filter
alpha_b0      = 0.0002836 * 0.0177 * 597/1412 *0.05603#0.00012625 * 0.0177 /0.05603   #0.0002836 * 0.0177
alpha_bplus   = 0.000258  * 0.0171 * 433/1285 *0.05943#0.00012125 * 0.0171 /0.05943  #0.000258  * 0.0171
alpha_bs      = 0.0000992 * 0.0144 * 173/494  *0.02642#0.00003    * 0.0144 /0.02642  #0.0000992 * 0.0144
alpha_lambdab = 0.00001125 * 0.04400  #0.0000544 * 0.04400 

# taken from Bs4mu analysis note br_fi = (1-(1-f_i)^2) 
br_fs = 0.189
br_fd = 0.676
br_fu = 0.676


#B0 meson (d - anti b)
# counted on 1 signal sample: 597 #(total events are 1412, total genmatched are 682, selected process is 213) 
eff_b0     = 0.0002836 * 597.0 / 1412.0
b0_to_cc   = 0.0177 #* 0.05603
# counted on whole exclusive B0 sample (all decays)
n_b0_to_cc = 1081753.0

b0_lumi = n_b0_to_cc / (eff_b0 * b0_to_cc)

#Bs meson (s - anti b)
# counted on 1 signal sample: 173 #(total events are 494, total genmatched are 187, selected process is 315) 
eff_bs     = 0.0000992 * 173.0 / 494.0
bs_to_cc   = 0.0144
# counted on whole exclusive Bs sample (all decays)
n_bs_to_cc = 549523

bs_lumi = n_bs_to_cc / (eff_bs * bs_to_cc)

#bs but now with another channel Bs -> Ds- D+
eff_bs     = 0.0001084 * 363.0/539.0 
bs_to_cc   = 0.00028  #* 0.02642
n_bs_to_cc = 23920.0

bs_lumi = n_bs_to_cc / (eff_bs * bs_to_cc)

#B+ meson (u - anti b)
#counted on 1 signal sample 433 #(total events are 1285, total genmatched are 504, selected process is 114) 
eff_bplus     = 0.000258 * 433 / 1285
bplus_to_cc   = 0.0171  #* 0.14045
# counted on whole exclusive B+ sample (all decays)
n_bplus_to_cc = 377432.0

bplus_lumi = n_bplus_to_cc / (eff_bplus * bplus_to_cc)



# constants from jira production (taken on 10.06.2024)
# no idea what we did here ???

#sigma_bb = 5.7e11
#
#n_b0   = 4585117 #9329044  
#eff_b0 = 3.33e-4
#br_b0  = 0.05603
#lumi_b0 = n_b0 / (eff_b0 * br_b0 * sigma_bb)
#
#n_bs   = 1096223 #2956307 
#eff_bs = 1.03e-4
#br_bs  = 0.02642
#lumi_bs = n_bs / (eff_bs * br_bs * sigma_bb)
#
#n_bplus   = 1792418 #4988382 
#eff_bplus = 1.67e-4
#br_bplus  = 0.05943
#lumi_bplus = n_bplus / (eff_bplus * br_bplus * sigma_bb)
#
## processed events and effective lumi
#n_processed_b0    = 4585432 #9280783
#n_processed_bs    = 1095772#2909717
#n_processed_bplus = 1876619#4866520
#
#lumi_b0_eff    = lumi_b0    * n_processed_b0    / n_b0
#lumi_bs_eff    = lumi_bs    * n_processed_bs    / n_bs
#lumi_bplus_eff = lumi_bplus * n_processed_bplus / n_bplus
#

# Selections



base = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 8)', 
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'(fv_prob > 0.1)',
'(tv_prob > 0.1)'
])

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

ma_cut_wout_tv = ' && '.join([
f'(dsMu_m < {bsMass_})',
'(k1_charge*k2_charge <0)',
'(mu_charge*pi_charge < 0)',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'(fv_prob > 0.1)'
])



ma_cut = ' && '.join([
f'(dsMu_m < {bsMass_})',
'(k1_charge*k2_charge <0)',
'(mu_charge*pi_charge < 0)',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'(fv_prob > 0.1)',
'(tv_prob > 0.1)'
])

ma_cut_high_mass = ' && '.join([
f'(dsMu_m > {bsMass_})',
'(k1_charge*k2_charge <0)',
'(mu_charge*pi_charge < 0)',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'(fv_prob > 0.1)',
'(tv_prob > 0.1)'
])

ma_cut_sign_flip = ' && '.join([
f'(dsMu_m < {bsMass_})',
'((k1_charge*k2_charge > 0) || (mu_charge*pi_charge > 0))',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'(fv_prob > 0.1)',
'(tv_prob > 0.1)'
])

ma_cut_sign_flip_high_mass = ' && '.join([
f'(dsMu_m > {bsMass_})',
'((k1_charge*k2_charge > 0) || (mu_charge*pi_charge > 0))',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'(fv_prob > 0.1)',
'(tv_prob > 0.1)'
])

###########################

ma_cut_wout_fv = ' & '.join([
f'(dsMu_m < {bsMass_})',
'(k1_charge*k2_charge <0)',
'(mu_charge*pi_charge < 0)',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
])

#baseline = ' & '.join([
#f'(dsMu_m < {bsMass_})',
#'(k1_charge*k2_charge <0)',
#'(mu_charge*pi_charge < 0)',
#'(mu_pt > 8)',
#'(k1_pt > 1)',
#'(k2_pt > 1)',
#'(pi_pt > 1)',
#'(lxy_ds < 1)',
#'(mu_id_medium == 1)',
#'(mu_rel_iso_03 < 0.3)',
#'(tv_prob > 0.1)',
#'((cosPiK1 < -0.3) || (cosPiK1 > 0.3))',
#'(fv_prob > 0.1)'
#])

addOn1 = ' & '.join([
'(bs_pt_reco_2 > 10)',
'(cosMuW_lhcb_alt > -0.5)',
'((cosPiK1 < -0.6) || (cosPiK1 > 0.6))',
'(e_star_reco_2 < 2.5)'
])

addOn2 = ' & '.join([
'(pt_miss_lhcb_alt < 25)',
'(m2_miss_lhcb_alt > -2)',
f'(abs(kk_m - {phiMass_}) < 0.010 )',
])

addOn3 = ' & '.join([
'(m2_miss_lhcb_alt > -2)',
'(e_star_reco_weighted < 1.7)',
'(bs_pt_reco_2 > 20)',
'((cosPiK1 < -0.4) || (cosPiK1 > 0.4))',
])

high_mass = ' & '.join([
f'dsMu_m > {bsMass_}',
#'mu_rel_iso_03 > 0.3'
])

flip_iso = ' & '.join([
#f'dsMu_m > {bsMass_}',
'mu_rel_iso_03 > 0.3'
])

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
    "data_cons_25": data_cons_25
    }
  print(json.dumps(config))

