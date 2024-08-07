
# samples

sig_unc     = ["26_07_2024_14_44_54"]
hb_unc      = ["25_07_2024_14_23_01"]
bs_unc      = ["30_07_2024_11_40_54"]
b0_unc      = ["30_07_2024_11_40_34"]
bplus_unc   = ["30_07_2024_11_40_45"]

               #part 1            #part 2            #part 3           #part 4
data_unc    = ["20240802_111807", "20240730_223445", "20240728_200528","20240729_090628" ]

sig_cons    = ["26_07_2024_14_46_03"]
hb_cons     = ["25_07_2024_14_23_42"]
bs_cons     = ["31_07_2024_10_07_17"]
b0_cons     = ["31_07_2024_10_07_43"]
bplus_cons  = ["31_07_2024_10_08_00"]
               #part 1
data_cons   = ["20240724_170443"]

#samples after NN (sgd)
sig_unc_pastNN     = "sig_pastNN_06Aug2024_11h58m45s"
hb_unc_pastNN      = "hb_pastNN_06Aug2024_11h56m50s"
bs_unc_pastNN      = "bs_pastNN_06Aug2024_11h56m53s"
b0_unc_pastNN      = "b0_pastNN_06Aug2024_11h56m57s"
bplus_unc_pastNN   = "bplus_pastNN_06Aug2024_11h56m53s"
data_unc_pastNN    = "data_pastNN_06Aug2024_11h59m43s"

#adam (mu and tau in class 0)
sig_unc_pastNN     = "sig_pastNN_06Aug2024_19h02m14s"
hb_unc_pastNN      = "hb_pastNN_06Aug2024_19h00m37s"
bs_unc_pastNN      = "bs_pastNN_06Aug2024_19h00m38s"
b0_unc_pastNN      = "b0_pastNN_06Aug2024_19h00m47s"
bplus_unc_pastNN   = "bplus_pastNN_06Aug2024_19h00m38s"
data_unc_pastNN    = "data_pastNN_06Aug2024_19h05m54s"

#adam (mu and hb in class 1)
sig_unc_pastNN     = "sig_pastNN_07Aug2024_14h04m03s"
hb_unc_pastNN      = "hb_pastNN_07Aug2024_14h03m05s"
bs_unc_pastNN      = "bs_pastNN_07Aug2024_14h03m09s"
b0_unc_pastNN      = "b0_pastNN_07Aug2024_14h03m18s"
bplus_unc_pastNN   = "bplus_pastNN_07Aug2024_14h03m08s"
data_unc_pastNN    = "data_pastNN_07Aug2024_14h09m37s"

#adam (mu and hb in class 1)
sig_cons_pastNN     = "sig_pastNN_07Aug2024_17h22m02s"
hb_cons_pastNN      = "hb_pastNN_07Aug2024_17h21m16s"
bs_cons_pastNN      = "bs_pastNN_07Aug2024_17h21m18s"
b0_cons_pastNN      = "b0_pastNN_07Aug2024_17h21m22s"
bplus_cons_pastNN   = "bplus_pastNN_07Aug2024_17h21m18s"
data_cons_pastNN    = "data_pastNN_07Aug2024_17h21m38s"


# Mass constants

dsMass_       = 1.96834
bsMass_       = 5.36688
phiMass_      = 1.019461

# Define sideband region

nSignalRegion = 3 # signal region is 3 sigma
nSidebands    = 5 # sideband region starts after 5 sigma
sbWidth       = 1 # sideband region is 1 sigma broad

# alpha = efficiency * f_i * BR 
# where f_i is the fragmentation for quark type i = u,s,d, ...
#               #epsilon    #BR                                             #the old ones without filter
alpha_b0      = 0.0002836 * 0.0177 * 597/1412 *0.05603#0.00012625 * 0.0177 /0.05603   #0.0002836 * 0.0177
alpha_bplus   = 0.000258  * 0.0171 * 433/1285 *0.05943#0.00012125 * 0.0171 /0.05943  #0.000258  * 0.0171
alpha_bs      = 0.0000992 * 0.0144 * 173/494  *0.02642#0.00003    * 0.0144 /0.02642  #0.0000992 * 0.0144
alpha_lambdab = 0.00001125 * 0.04400  #0.0000544 * 0.04400 


# constants from jira production (taken on 10.06.2024)
sigma_bb = 5.7e11

n_b0   = 9329044  
eff_b0 = 3.33e-4
br_b0  = 0.05603
lumi_b0 = n_b0 / (eff_b0 * br_b0 * sigma_bb)

n_bs   = 2956307 
eff_bs = 1.03e-4
br_bs  = 0.02642
lumi_bs = n_bs / (eff_bs * br_bs * sigma_bb)

n_bplus   = 4988382 
eff_bplus = 1.67e-4
br_bplus  = 0.05943
lumi_bplus = n_bplus / (eff_bplus * br_bplus * sigma_bb)

# processed events and effective lumi
n_processed_b0    = 9280783
n_processed_bs    = 2909717
n_processed_bplus = 4866520

lumi_b0_eff    = lumi_b0    * n_processed_b0    / n_b0
lumi_bs_eff    = lumi_bs    * n_processed_bs    / n_bs
lumi_bplus_eff = lumi_bplus * n_processed_bplus / n_bplus


# Selections

#########################
## work with this sel  ##
#########################

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

base_wout_tv = ' && '.join([ #remove the charge and ds+mu mass cuts!
f'(mu_pt > 8)', 
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
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

baseline = ' & '.join([
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
'(tv_prob > 0.1)',
'((cosPiK1 < -0.3) || (cosPiK1 > 0.3))',
'(fv_prob > 0.1)'
])

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

bkg = ' & '.join([
f'dsMu_m > {bsMass_}',
#'mu_rel_iso_03 > 0.3'
])

