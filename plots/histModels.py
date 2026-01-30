import ROOT
import yaml
import numpy as np
from array import array

##############################################################################
## Define the structure of your histograms which are plotted in plotNTuples ##
##############################################################################



models = {}
modelsSR = {}
pastNN_models = {}
pastNN_2Dmodels = {}

##masses 
models["phiPi_m" ]                   = (ROOT.RDF.TH1DModel("phiPi_m"                           , '',   30,   1.91, 2.028), r"m(KK#pi) (GeV)"                           ,    0)
#models["kk_m"    ]                   = (ROOT.RDF.TH1DModel("kk_m"                              , '',   30,   1.0,  1.040), r"m(KK) (GeV)"                           ,    0)
models["dsMu_m" ]                   = (ROOT.RDF.TH1DModel("dsMu_m"                             , '',   25,      2,   8), r" m(KK#pi#mu) (GeV)"                       ,    0)
#models["dsMu_m" ]                   = (ROOT.RDF.TH1DModel("dsMu_m"                             , '',   30,      2,   5.366), r"m(KK#pi#mu)(GeV)"                       ,    0)
#models["dsMu_m" ]                   = (ROOT.RDF.TH1DModel("dsMu_m"                             , '',   30,      5.366,   8), r"m(KK#pi#mu)(GeV)"                       ,    0)
#
#models["mu_pt" ]                     = (ROOT.RDF.TH1DModel("mu_pt"                             , '',   21,      8,    25), r" #mu(p_{T}) (GeV)"                         ,    0)
#models["pi_pt" ]                     = (ROOT.RDF.TH1DModel("pi_pt"                             , '',   21,      1,    12), r" #pi(p_{T}) (GeV)"                         ,    0)
#models["k1_pt" ]                     = (ROOT.RDF.TH1DModel("k1_pt"                             , '',   21,      1,    12), r" K1(p_{T}) (GeV)"                         ,    0)
#models["k2_pt" ]                     = (ROOT.RDF.TH1DModel("k2_pt"                             , '',   21,      1,    12), r" K2(p_{T}) (GeV)"                         ,    0)
#
#
#models["mu_eta" ]                     = (ROOT.RDF.TH1DModel("mu_eta"                             , '',   21,      0,   2.5), r" #eta(#mu)"                         ,    0)
#models["pi_eta" ]                     = (ROOT.RDF.TH1DModel("pi_eta"                             , '',   21,      0,   2.5), r" #eta(#pi)"                         ,    0)
#models["k1_eta" ]                     = (ROOT.RDF.TH1DModel("k1_eta"                             , '',   21,      0,   2.5), r" #eta(K1)"                         ,    0)
#models["k2_eta" ]                     = (ROOT.RDF.TH1DModel("k2_eta"                             , '',   21,      0,   2.5), r" #eta(K2)"                         ,    0)
#
#
#models["m2_miss_coll" ]              = (ROOT.RDF.TH1DModel("m2_miss_coll"                      , '',   21,      0,     6), r"m^{2}_{miss,coll} (GeV^{2})"                     ,    0)
#models["m2_miss_lhcb_alt" ]          = (ROOT.RDF.TH1DModel("m2_miss_lhcb_alt"                  , '',   21,      0,     6), r"m^{2}_{miss,xyz} (GeV^{2})"                     ,    0)
#models["m2_miss_reco_1" ]            = (ROOT.RDF.TH1DModel("m2_miss_reco_1"                    , '',   21,      0,     6), r"m^{2}_{miss,math,1} (GeV^{2})"                     ,    0)
#models["m2_miss_reco_2" ]            = (ROOT.RDF.TH1DModel("m2_miss_reco_2"                    , '',   21,      0,     6), r"m^{2}_{miss,math,2} (GeV^{2})"                     ,    0)
#
#
#models["pt_miss_coll" ]              = (ROOT.RDF.TH1DModel("pt_miss_coll"                      , '',   21,      0,    30), r"p_{T}^{miss}_{coll} (GeV)"                  ,    0)
#models["pt_miss_lhcb_alt" ]          = (ROOT.RDF.TH1DModel("pt_miss_lhcb_alt"                  , '',   21,      0,    30), r"p_{T}^{miss}_{xyz} (GeV)"                   ,    0)
#models["pt_miss_reco_1" ]            = (ROOT.RDF.TH1DModel("pt_miss_reco_1"                    , '',   21,      0,    30), r"p_{T}^{miss}_{math,1} (GeV)"                ,    0)
#models["pt_miss_reco_2" ]            = (ROOT.RDF.TH1DModel("pt_miss_reco_2"                    , '',   21,      0,    30), r"p_{T}^{miss}_{math,2} (GeV)"                ,    0)
#
#
#models["bs_pt_coll" ]                = (ROOT.RDF.TH1DModel("bs_coll_pt"                        , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#models["bs_pt_lhcb_alt" ]            = (ROOT.RDF.TH1DModel("bs_pt_lhcb_alt"                    , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#models["bs_pt_reco_1" ]              = (ROOT.RDF.TH1DModel("bs_pt_reco_1"                      , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#models["bs_pt_reco_2" ]              = (ROOT.RDF.TH1DModel("bs_pt_reco_2"                      , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
##
models["q2_coll" ]               = (ROOT.RDF.TH1DModel("q2_coll"                       , '',   31,      0,     12), r"q^{2}_{coll} (GeV)"                             ,    0)
#models["q2_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("q2_lhcb_alt"                   , '',   31,      0,     12), r"q^{2}_{xyz} (GeV)"                             ,    0)
#models["q2_reco_1" ]             = (ROOT.RDF.TH1DModel("q2_reco_1"                     , '',   31,      0,     12), r"q^{2}_{math,1} (GeV)"                             ,    0)
#models["q2_reco_2" ]             = (ROOT.RDF.TH1DModel("q2_reco_2"                     , '',   31,      0,     12), r"q^{2}_{math,2} (GeV)"                             ,    0)
#models["q2_coll" ]               = (ROOT.RDF.TH1DModel("q2_coll"                       , '',   31,      -12,     12), r"q^{2}_{coll} (GeV)"                             ,    0)
#models["q2_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("q2_lhcb_alt"                   , '',   31,      -12,     12), r"q^{2}_{xyz} (GeV)"                             ,    0)
#models["q2_reco_1" ]             = (ROOT.RDF.TH1DModel("q2_reco_1"                     , '',   31,      -12,     12), r"q^{2}_{math,1} (GeV)"                             ,    0)
#models["q2_reco_2" ]             = (ROOT.RDF.TH1DModel("q2_reco_2"                     , '',   31,      -12,     12), r"q^{2}_{math,2} (GeV)"                             ,    0)
#
#
#models["e_star_coll" ]               = (ROOT.RDF.TH1DModel("e_star_coll"                       , '',   31,      0,     3), r"E*_{coll} (GeV)"                             ,    0)
#models["e_star_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("e_star_lhcb_alt"                   , '',   31,      0,     3), r"E*_{xyz} (GeV)"                             ,    0)
#models["e_star_reco_1" ]             = (ROOT.RDF.TH1DModel("e_star_reco_1"                     , '',   31,      0,     3), r"E*_{math,1} (GeV)"                             ,    0)
#models["e_star_reco_2" ]             = (ROOT.RDF.TH1DModel("e_star_reco_2"                     , '',   31,      0,     3), r"E*_{math,2} (GeV)"                             ,    0)
#
#
#models["e_gamma" ]                   = (ROOT.RDF.TH1DModel("e_gamma"                           , '',   21,      0,     4), r"E^{#gamma}_{coll} (GeV)"                           ,    0)
#models["disc_negativity" ]           = (ROOT.RDF.TH1DModel("disc_negativity"                   , '',   21,      0,     6), r"Disc_{#le 0}"                                   ,    0)
#
#
#models["kk_deltaR" ]                 = (ROOT.RDF.TH1DModel("kk_deltaR"                         , '',   21,      0,   0.3), r"#Delta R(K,K)"                            ,    0)
#models["phiPi_deltaR" ]              = (ROOT.RDF.TH1DModel("phiPi_deltaR"                      , '',   21,      0,     1), r"#Delta R(KK,#pi)"                         ,    0)
#models["dsMu_deltaR" ]               = (ROOT.RDF.TH1DModel("dsMu_deltaR"                       , '',   21,      0,   0.8), r"#Delta R(KK#pi,#mu)"                        ,    0)
#
#
#models["bs_boost_coll" ]             = (ROOT.RDF.TH1DModel("bs_boost_coll"                     , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
#models["bs_boost_lhcb_alt" ]         = (ROOT.RDF.TH1DModel("bs_boost_lhcb_alt"                 , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
#models["bs_boost_reco_1" ]           = (ROOT.RDF.TH1DModel("bs_boost_reco_1"                   , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
#models["bs_boost_reco_2" ]           = (ROOT.RDF.TH1DModel("bs_boost_reco_2"                   , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
#
#models["cosMuW_coll" ]               = (ROOT.RDF.TH1DModel("cosMuWColl"                        , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{coll}"                           ,    0)
#models["cosMuW_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("cosMuWLhcbAlt"                     , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{xyz}"                            ,    0)
#models["cosMuW_reco_1" ]             = (ROOT.RDF.TH1DModel("cosMuWReco1"                       , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{math,1}"                         ,    0)
#models["cosMuW_reco_2" ]             = (ROOT.RDF.TH1DModel("cosMuWReco2"                       , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{math,2}"                         ,    0)
#
#models["cosPhiDs_coll" ]             = (ROOT.RDF.TH1DModel("cosPhiDsColl"                      , '',   21,     -1,     1), r"cos(#theta^{\Phi, D_{s}}_{hel})_{coll}"                            ,    0)
#models["cosPhiDs_lhcb_alt" ]         = (ROOT.RDF.TH1DModel("cosPhiDsLhcbAlt"                   , '',   21,     -1,     1), r"cos(#theta^{\Phi, D_{s}}_{hel})_{xyz}"                            ,    0)
#models["cosPhiDs_reco_1" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco1"                     , '',   21,     -1,     1), r"cos(#theta^{\Phi, D_{s}}_{hel})_{math,1}"                            ,    0)
#models["cosPhiDs_reco_2" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco2"                     , '',   21,     -1,     1), r"cos(#theta^{\Phi, D_{s}}_{hel})_{math,2}"                            ,    0)
##
##models["cosPlaneBs_coll" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsColl"                    , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
##models["cosPlaneBs_lhcb_alt" ]       = (ROOT.RDF.TH1DModel("cosPlaneBsLhcbAlt"                 , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
##models["cosPlaneBs_reco_1" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco1"                   , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
##models["cosPlaneBs_reco_2" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco2"                   , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
##
##
#models["lxy_ds_sig" ]                = (ROOT.RDF.TH1DModel("lxy_ds_sig"                        , '',   21,      0,    10), r"L^{sig}_{xy}(D_{s})"                        ,    0)
#models["cosPiK1" ]                   = (ROOT.RDF.TH1DModel("cosPiK1"                           , '',   21,     -1,     1), r"cos(#theta_{4})"                            ,    0)
#
#models["sv_chi2" ]                   = (ROOT.RDF.TH1DModel("sv_chi2"                           , '',   30,    0.0,  10.0), r"#chi^{2}_{Vtx} (V_{B^{0}_{s}}) "                           ,    0)
#models["tv_chi2" ]                   = (ROOT.RDF.TH1DModel("tv_chi2"                           , '',   30,    0.0,  10.0), r"#chi^{2}_{Vtx} (V_{D_{s}}) "                           ,    0)
#models["fv_chi2" ]                   = (ROOT.RDF.TH1DModel("fv_chi2"                           , '',   30,    0.0,  10.0), r"#chi^{2}_{Vtx} (V_{#Phi}) "                            ,    0)
#
#
#models["rel_iso_03_pv" ]                   = (ROOT.RDF.TH1DModel("rel_iso_03_pv"                           , '',   30,    0.0,  0.3), r"Iso^{rel}_{#mu}(PV)"                  ,    0)
#
#models["dsPhoton_m" ]                      = (ROOT.RDF.TH1DModel("dsPhoton_m"                              , '',   30,    2.0,  2.3), r"m(Ds + #gamma)"                     ,    0)
#models["photon_pt" ]                      = (ROOT.RDF.TH1DModel("photon_pt"                              , '',   30,    0.0,  2.0), r"p_{T}(#gamma) (GeV)"                     ,    0)
#
#models["signed_decay_ip3d_mu_ds_sig_sv" ]             = (ROOT.RDF.TH1DModel("signed_decay_ip3d_mu_ds_sv"                     , '',   30,    -0.5, 0.5), r"IP3D"                      ,    0)
#
#models["ds_perp" ]                      = (ROOT.RDF.TH1DModel("ds_perp"                              , '',   30,    0,  5), r"m_{#perp}(KK #pi) (GeV)"                     ,    0)
#models["ds_mu_perp" ]                   = (ROOT.RDF.TH1DModel("ds_mu_perp"                           , '',   30,    0,  8), r"m_{#perp}(KK #pi #mu) (GeV)"                 ,    0)
#models["ds_perp_photon" ]               = (ROOT.RDF.TH1DModel("ds_perp_photon"                       , '',   30,    0,  5), r"m_{#perp}(KK #pi #gamma) (GeV)"              ,    0)
#models["ds_mu_perp_photon" ]            = (ROOT.RDF.TH1DModel("ds_mu_perp_photon"                    , '',   30,    0,  8), r"m_{#perp}(KK #pi #mu #gamma) (GeV)"          ,    0)
#
#models["bs_mass_corr" ]            = (ROOT.RDF.TH1DModel("bs_mass_corr"                    , '',   30,    2.0,  12), r"m_{corr}(KK #pi #mu ) (GeV)"          ,    0)
#models["bs_mass_corr_photon" ]     = (ROOT.RDF.TH1DModel("bs_mass_corr_photon"             , '',   30,    2.0,  12), r"m_{corr}(KK #pi #mu #gamma) (GeV)"          ,    0)
#

"""
models["gen_pv_x" ]                  = (ROOT.RDF.TH1DModel("gen_pv_x"                        , '', 30,    0.0,  0.02), r"(PV)_{x} (cm)"                            ,    0)
models["gen_pv_y" ]                  = (ROOT.RDF.TH1DModel("gen_pv_y"                        , '', 30,   0.02,  0.06), r"(PV)_{y} (cm)"                            ,    0)
models["gen_pv_z" ]                  = (ROOT.RDF.TH1DModel("gen_pv_z"                        , '', 30,     -5,     5), r"(PV)_{z} (cm)"                            ,    0)
models["pv_x" ]                      = (ROOT.RDF.TH1DModel("pv_x"                            , '', 30,    0.0,  0.02), r"(PV)_{x} (cm)"                            ,    0)
models["pv_y" ]                      = (ROOT.RDF.TH1DModel("pv_y"                            , '', 30,   0.02,    12), r"(PV)_{y} (cm)"                            ,    0)
models["pv_z" ]                      = (ROOT.RDF.TH1DModel("pv_z"                            , '', 30,     -5,     5), r"(PV)_{z} (cm)"                            ,    0)

models["gen_sv_x" ]                  = (ROOT.RDF.TH1DModel("gen_sv_x"                        , '', 30,     -1,     1), r"(SV)_{x} (cm)"                            ,    0)
models["gen_sv_y" ]                  = (ROOT.RDF.TH1DModel("gen_sv_y"                        , '', 30,     -1,     1), r"(SV)_{y} (cm)"                            ,    0)
models["gen_sv_z" ]                  = (ROOT.RDF.TH1DModel("gen_sv_z"                        , '', 30,     -5,     5), r"(SV)_{z} (cm)"                            ,    0)
models["sv_x" ]                      = (ROOT.RDF.TH1DModel("sv_x"                            , '', 30,     -1,     1), r"(SV)_{x} (cm)"                            ,    0)
models["sv_y" ]                      = (ROOT.RDF.TH1DModel("sv_y"                            , '', 30,     -1,     1), r"(SV)_{y} (cm)"                            ,    0)
models["sv_z" ]                      = (ROOT.RDF.TH1DModel("sv_z"                            , '', 30,     -5,     5), r"(SV)_{z} (cm)"                            ,    0)

models["gen_tv_x" ]                  = (ROOT.RDF.TH1DModel("gen_tv_x"                        , '', 30,     -1,     1), r"(TV)_{x} (cm)"                            ,    0)
models["gen_tv_y" ]                  = (ROOT.RDF.TH1DModel("gen_tv_y"                        , '', 30,     -1,     1), r"(TV)_{y} (cm)"                            ,    0)
models["gen_tv_z" ]                  = (ROOT.RDF.TH1DModel("gen_tv_z"                        , '', 30,     -5,     5), r"(TV)_{z} (cm)"                            ,    0)
models["tv_x" ]                      = (ROOT.RDF.TH1DModel("tv_x"                            , '', 30,     -1,     1), r"(TV)_{x} (cm)"                            ,    0)
models["tv_y" ]                      = (ROOT.RDF.TH1DModel("tv_y"                            , '', 30,     -1,     1), r"(TV)_{y} (cm)"                            ,    0)
models["tv_z" ]                      = (ROOT.RDF.TH1DModel("tv_z"                            , '', 30,     -5,     5), r"(TV)_{z} (cm)"                            ,    0)
 
models["gen_fv_x" ]                  = (ROOT.RDF.TH1DModel("gen_fv_x"                        , '', 30,     -1,     1), r"(FV)_{x} (cm)"                            ,    0)
models["gen_fv_y" ]                  = (ROOT.RDF.TH1DModel("gen_fv_y"                        , '', 30,     -1,     1), r"(FV)_{y} (cm)"                            ,    0)
models["gen_fv_z" ]                  = (ROOT.RDF.TH1DModel("gen_fv_z"                        , '', 30,     -5,     5), r"(FV)_{z} (cm)"                            ,    0)
models["fv_x" ]                      = (ROOT.RDF.TH1DModel("fv_x"                            , '', 30,     -1,     1), r"(FV)_{x} (cm)"                            ,    0)
models["fv_y" ]                      = (ROOT.RDF.TH1DModel("fv_y"                            , '', 30,     -1,     1), r"(FV)_{y} (cm)"                            ,    0)
models["fv_z" ]                      = (ROOT.RDF.TH1DModel("fv_z"                            , '', 30,     -5,     5), r"(FV)_{z} (cm)"                            ,    0)


models["bs_fitted_pt" ]           = (ROOT.RDF.TH1DModel("bs_fitted_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["dsMu_pt" ]                = (ROOT.RDF.TH1DModel("dsMu_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)

models["gen_phi_pt" ]             = (ROOT.RDF.TH1DModel("gen_phi_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["phi_fitted_pt" ]          = (ROOT.RDF.TH1DModel("phi_fitted_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["phi_refitted_pt" ]        = (ROOT.RDF.TH1DModel("phi_refitted_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["kk_pt" ]                  = (ROOT.RDF.TH1DModel("kk_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)

models["gen_ds_pt" ]              = (ROOT.RDF.TH1DModel("gen_ds_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["ds_fitted_pt" ]           = (ROOT.RDF.TH1DModel("ds_fitted_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["ds_refitted_pt" ]         = (ROOT.RDF.TH1DModel("ds_refitted_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["phiPi_pt" ]               = (ROOT.RDF.TH1DModel("phiPi_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)

"""



#pastNN_models["class" ]                    = (ROOT.RDF.TH1DModel("class"                            , '',   6,       0,    5.9999), r"Class"                            ,    0)

pastNN_models["score0" ]                   = (ROOT.RDF.TH1DModel("score0"                           , '',   31,       0,    1), r"Score 0"                            ,    0)
pastNN_models["score1" ]                   = (ROOT.RDF.TH1DModel("score1"                           , '',   31,       0,    0.7), r"Score 1"                            ,    0)
pastNN_models["score2" ]                   = (ROOT.RDF.TH1DModel("score2"                           , '',   31,       0,    1), r"Score 2"                            ,    0)
pastNN_models["score3" ]                   = (ROOT.RDF.TH1DModel("score3"                           , '',   31,       0,    1), r"Score 3"                            ,    0)
pastNN_models["score4" ]                   = (ROOT.RDF.TH1DModel("score4"                           , '',   31,       0,    1), r"Score 4"                            ,    0)
pastNN_models["score5" ]                   = (ROOT.RDF.TH1DModel("score5"                           , '',   50,       0,    0.1), r"Score 5"                            ,    0)
#pastNN_models["score5" ]                = (ROOT.RDF.TH1DModel("score5_WP_lower"                           , '',   100,      0,    0.1), r"Score 5"                            ,    0)
#pastNN_models["score0" ]                = (ROOT.RDF.TH1DModel("score0_WP_upper"                           , '',   200,      0,    0.2), r"Score 0"                            ,    0)
#
#pastNN_2Dmodels["score0_score1" ]           = (ROOT.RDF.TH2DModel("score0_score1"                    , '',   2 ,       0,    0.6, 2, 0, 0.6), r"Score 0x1"                            ,    0)
#pastNN_2Dmodels["score2_score1" ]           = (ROOT.RDF.TH2DModel("score2_score1"                    , '',   4 ,       0,    1  , 4, 0, 1)  , r"Score 2x1"                            ,    0)
#pastNN_2Dmodels["score3_score1" ]           = (ROOT.RDF.TH2DModel("score3_score1"                    , '',   4 ,       0,    1  , 4, 0, 1)  , r"Score 3x1"                            ,    0)


#save this also into yaml file for hammer :)
to_export = {}
for name, model in models.items():
  to_export[name] = {"bins": model[0].fNbinsX, "xmin": model[0].fXLow, "xmax": model[0].fXUp }

for name, model in pastNN_models.items():
  to_export[name] = {"bins": model[0].fNbinsX, "xmin": model[0].fXLow, "xmax": model[0].fXUp }

with open("/work/pahwagne/RDsTools/hammercpp/development_branch/weights/plottingModels.yaml", "w") as f:
    yaml.dump(to_export, f)

print("Included models are:", models.keys())
print("Included past NN models are:", pastNN_models.keys())


#here we define special binnings for the final fit
#binning of score3 when fitted in bins of class

special_models = {}
special_models["score1_bin0" ]                   = (ROOT.RDF.TH1DModel("score1_bin0"                           , '',   31,       0  ,    0.45), r"Score 1"                            ,    0)
special_models["score1_bin1" ]                   = (ROOT.RDF.TH1DModel("score1_bin1"                           , '',   31,       0.2,    0.8 ), r"Score 1"                            ,    0)
special_models["score1_bin2" ]                   = (ROOT.RDF.TH1DModel("score1_bin2"                           , '',   31,       0  ,    0.2 ), r"Score 1"                            ,    0)
special_models["score1_bin3" ]                   = (ROOT.RDF.TH1DModel("score1_bin3"                           , '',   31,       0  ,    0.45), r"Score 1"                            ,    0)
special_models["score1_bin4" ]                   = (ROOT.RDF.TH1DModel("score1_bin4"                           , '',   31,       0  ,    0.45), r"Score 1"                            ,    0)
special_models["score1_bin5" ]                   = (ROOT.RDF.TH1DModel("score1_bin5"                           , '',   31,       0  ,    1   ), r"Score 1"                            ,    0)

special_models_q2_coll = {}

######### 2D fit ############

# bin 0,1
#edges = list(np.linspace(0,0.04,20)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin0" ]                   = (ROOT.RDF.TH1DModel("score1_bin0"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
##bin 2,3
#edges = list(np.linspace(0,0.11,20)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin1" ]                   = (ROOT.RDF.TH1DModel("score1_bin1"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
##bin 4,5
#edges = list(np.linspace(0,0.21,20)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin2" ]                   = (ROOT.RDF.TH1DModel("score1_bin2"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
##bin 6,7
#edges = list(np.linspace(0,0.35,20)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin3" ]                   = (ROOT.RDF.TH1DModel("score1_bin3"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
##bin 8,9
#edges = list(np.linspace(0,0.5,20)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin4" ]                   = (ROOT.RDF.TH1DModel("score1_bin4"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
##bin 10,11
#edges = list(np.linspace(0,0.6,20)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin5" ]                   = (ROOT.RDF.TH1DModel("score1_bin5"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
##bin 12,13
#edges = list(np.linspace(0,0.7,12)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin6" ]                   = (ROOT.RDF.TH1DModel("score1_bin6"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
##bin 14,15
#edges = list(np.linspace(0,0.7,12)) 
#edgescpp = array('d', edges)
#special_models_q2_coll["score1_bin7" ]                   = (ROOT.RDF.TH1DModel("score1_bin7"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)
#
#
#
#




############ 3D fit in the mass ###########

#bin 0,1
edges = list(np.linspace(0,0.04,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin0" ]                   = (ROOT.RDF.TH1DModel("score1_bin0"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.04,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin1" ]                   = (ROOT.RDF.TH1DModel("score1_bin1"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

#bin 2,3
edges = list(np.linspace(0,0.11,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin2" ]                   = (ROOT.RDF.TH1DModel("score1_bin2"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.11,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin3" ]                   = (ROOT.RDF.TH1DModel("score1_bin3"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

#bin 4,5
edges = list(np.linspace(0,0.21,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin4" ]                   = (ROOT.RDF.TH1DModel("score1_bin4"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.21,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin5" ]                   = (ROOT.RDF.TH1DModel("score1_bin5"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

#bin 6,7
edges = list(np.linspace(0,0.35,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin6" ]                   = (ROOT.RDF.TH1DModel("score1_bin6"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.35,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin7" ]                   = (ROOT.RDF.TH1DModel("score1_bin7"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

#bin 8,9
edges = list(np.linspace(0,0.5,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin8" ]                   = (ROOT.RDF.TH1DModel("score1_bin8"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.5,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin9" ]                   = (ROOT.RDF.TH1DModel("score1_bin9"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

#bin 10,11
edges = list(np.linspace(0,0.6,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin10" ]                   = (ROOT.RDF.TH1DModel("score1_bin10"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.6,20)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin11" ]                   = (ROOT.RDF.TH1DModel("score1_bin11"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

#bin 12,13
edges = list(np.linspace(0,0.7,12)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin12" ]                   = (ROOT.RDF.TH1DModel("score1_bin12"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.7,12)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin13" ]                   = (ROOT.RDF.TH1DModel("score1_bin13"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

#bin 14,15
edges = list(np.linspace(0,0.7,12)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin14" ]                   = (ROOT.RDF.TH1DModel("score1_bin14"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)

edges = list(np.linspace(0,0.7,12)) 
edgescpp = array('d', edges)
special_models_q2_coll["score1_bin15" ]                   = (ROOT.RDF.TH1DModel("score1_bin15"                           , '', len(edges)-1, edgescpp  ), r"Score 1"                            ,    0)





















