import ROOT
import yaml

##############################################################################
## Define the structure of your histograms which are plotted in plotNTuples ##
##############################################################################



models = {}
pastNN_models = {}
pastNN_2Dmodels = {}
## 
models["phiPi_m" ]                   = (ROOT.RDF.TH1DModel("phiPi_m"                           , '',   30,   1.91, 2.028), r"D_{s} mass (GeV)"                           ,    0)
 
#models["gen_mu_pt" ]                 = (ROOT.RDF.TH1DModel("gen_mu_pt"                       , '',   21,      0,    30), r"p_{T} (GeV)"                         ,    0)
#models["mu_pt" ]                     = (ROOT.RDF.TH1DModel("mu_pt"                             , '',   21,      0,    30), r" #mu(p_{T}) (GeV)"                         ,    0)
#models["pi_pt" ]                     = (ROOT.RDF.TH1DModel("pi_pt"                             , '',   21,      0,    30), r" #pi(p_{T}) (GeV)"                         ,    0)
#models["kk_pt" ]                     = (ROOT.RDF.TH1DModel("kk_pt"                             , '',   21,      0,    30), r" #kk(p_{T}) (GeV)"                         ,    0)
#
## Gens
#models["gen_q2" ]                    = (ROOT.RDF.TH1DModel("gen_q2"                            , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#models["gen_m2_miss" ]               = (ROOT.RDF.TH1DModel("gen_m2_miss"                       , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
#models["gen_pt_miss" ]               = (ROOT.RDF.TH1DModel("gen_pt_miss"                       , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
#models["gen_bs_pt" ]                 = (ROOT.RDF.TH1DModel("gen_bs_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#models["gen_e_star" ]                = (ROOT.RDF.TH1DModel("gen_e_star"                        , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
#models["gen_e_gamma" ]               = (ROOT.RDF.TH1DModel("gen_e_gamma"                       , '',   21,      0,     4), r"E_{#gamma} (GeV)"                           ,    0)
#models["gen_mu_rel_iso" ]            = (ROOT.RDF.TH1DModel("gen_mu_rel_iso"                    , '',   21,      0,    20), r"Iso^{rel}_{#mu}"                            ,    0)
#models["gen_cosMuW" ]                = (ROOT.RDF.TH1DModel("gen_cosMuWGen"                     , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
#models["gen_cosPhiDs" ]              = (ROOT.RDF.TH1DModel("gen_cosPhiDsGen"                   , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
#models["gen_cosPlaneBs" ]            = (ROOT.RDF.TH1DModel("gen_cosPlaneBsGen"                 , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
#models["gen_lxy_ds_sig" ]            = (ROOT.RDF.TH1DModel("gen_lxy_ds_sig"                    , '',   21,     -1,     2), r"cos(#chi)"                                  ,    0)
#models["gen_dxy_mu_sig" ]            = (ROOT.RDF.TH1DModel("gen_dxy_mu_sig"                    , '',   21,    -50,    50), r"cos(#chi)"                                  ,    0)
#models["gen_cosPiK1" ]               = (ROOT.RDF.TH1DModel("cosPiK1"                           , '',   21,     -1,     1), r"cos(#theta_{4})"                            ,    0)
#
## Basics
#
#models["dsMu_m" ]                   = (ROOT.RDF.TH1DModel("dsMu_m"                             , '',   30,      2,   8), r"D_{s}+#mu mass (GeV)"                       ,    0)
##models["dsMu_pt" ]                   = (ROOT.RDF.TH1DModel("dsMu_m"                             , '',   30,      0,  70), r"p_{T}(D_{s}+#mu) (GeV)"                       ,    0)
###
#models["q2_coll" ]                   = (ROOT.RDF.TH1DModel("q2_coll"                           , '',   21,       0,    12), r"q^{2}_{coll} (GeV^{2})"                          ,    0)
#models["q2_coll" ]                   = (ROOT.RDF.TH1DModel("q2_coll"                           , '',   21,       -12,    12), r"q^{2}_{coll} (GeV^{2})"                          ,    0)
####models["gen_q2" ]                   = (ROOT.RDF.TH1DModel("gen_q2"                           , '',   21,       0,    12), r"q^{2} (GeV^{2})"                            ,    0)
####models["q2_lhcb" ]                   = (ROOT.RDF.TH1DModel("q2_lhcb"                           , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#models["q2_lhcb_alt" ]               = (ROOT.RDF.TH1DModel("q2_lhcb_alt"                       , '',   21,      0,    12), r"q^{2}_{xyz} (GeV^{2})"                            ,    0)
#models["q2_lhcb_alt" ]               = (ROOT.RDF.TH1DModel("q2_lhcb_alt"                       , '',   21,      -12,    12), r"q^{2}_{xyz} (GeV^{2})"                            ,    0)
##models["q2_reco_1" ]                 = (ROOT.RDF.TH1DModel("q2_reco_1"                         , '',   21,      0,    12), r"q^{2}_{math,1} (GeV^{2})"                         ,    0)
##models["q2_reco_2" ]                 = (ROOT.RDF.TH1DModel("q2_reco_2"                         , '',   21,      0,    12), r"q^{2}_{math,2} (GeV^{2})"                         ,    0)
###models["q2_reco_weighted" ]          = (ROOT.RDF.TH1DModel("q2_reco_weighted"                  , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#######
##models["m2_miss_coll" ]              = (ROOT.RDF.TH1DModel("m2_miss_coll"                      , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
#####models["m2_miss_lhcb" ]              = (ROOT.RDF.TH1DModel("m2_miss_lhcb"                      , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
##models["m2_miss_lhcb_alt" ]          = (ROOT.RDF.TH1DModel("m2_miss_lhcb_alt"                  , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
##models["m2_miss_reco_1" ]          = (ROOT.RDF.TH1DModel("m2_miss_reco_1"                  , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
##models["m2_miss_reco_2" ]          = (ROOT.RDF.TH1DModel("m2_miss_reco_2"                  , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
######
##models["pt_miss_coll" ]              = (ROOT.RDF.TH1DModel("pt_miss_coll"                      , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
####models["pt_miss_lhcb" ]              = (ROOT.RDF.TH1DModel("pt_miss_lhcb"                      , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
##models["pt_miss_lhcb_alt" ]          = (ROOT.RDF.TH1DModel("pt_miss_lhcb_alt"                  , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
###models["pt_miss_reco_1" ]          = (ROOT.RDF.TH1DModel("pt_miss_reco_1"                  , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
###models["pt_miss_reco_2" ]          = (ROOT.RDF.TH1DModel("pt_miss_reco_2"                  , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
###models["pt_miss_reco_weighted" ]          = (ROOT.RDF.TH1DModel("pt_miss_lhcb_alt"                  , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
######
#models["bs_pt_coll" ]                = (ROOT.RDF.TH1DModel("bs_coll_pt"                        , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#######models["bs_pt_lhcb" ]                = (ROOT.RDF.TH1DModel("bs_pt_lhcb"                        , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#models["bs_pt_lhcb_alt" ]            = (ROOT.RDF.TH1DModel("bs_pt_lhcb_alt"                    , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
##models["bs_pt_reco_1" ]              = (ROOT.RDF.TH1DModel("bs_pt_reco_1"                      , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
##models["bs_pt_reco_2" ]              = (ROOT.RDF.TH1DModel("bs_pt_reco_2"                      , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#####models["bs_pt_reco_weighted" ]       = (ROOT.RDF.TH1DModel("bs_pt_reco_weighted"               , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
######
##models["e_star_coll" ]               = (ROOT.RDF.TH1DModel("e_star_coll"                       , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
#######models["e_star_lhcb" ]               = (ROOT.RDF.TH1DModel("e_star_lhcb"                       , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
##models["e_star_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("e_star_lhcb_alt"                   , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
##models["e_star_reco_1" ]             = (ROOT.RDF.TH1DModel("e_star_reco_1"                     , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
##models["e_star_reco_2" ]             = (ROOT.RDF.TH1DModel("e_star_reco_2"                     , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
##models["e_star_reco_weighted" ]      = (ROOT.RDF.TH1DModel("e_star_reco_weighted"              , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
######
##models["e_gamma" ]                   = (ROOT.RDF.TH1DModel("e_gamma"                           , '',   21,      0,     4), r"E_{#gamma} (GeV)"                           ,    0)
####models["mu_rel_iso_03" ]             = (ROOT.RDF.TH1DModel("mu_rel_iso_03"                     , '',   21,      0,   0.3), r"Iso^{rel}_{#mu}"                            ,    0)
######models["disc_is_negative" ]          = (ROOT.RDF.TH1DModel("disc_is_negative"                  , '',   21,     -1,     2), r"Disc < 0"                                   ,    0)
####models["disc_negativity" ]           = (ROOT.RDF.TH1DModel("disc_negativity"                   , '',   21,      0,     6), r"Disc < 0"                                   ,    0)
#####
#models["kk_deltaR" ]                 = (ROOT.RDF.TH1DModel("kk_deltaR"                         , '',   21,      0,   0.3), r"#Delta R(K,K)"                            ,    0)
#models["phiPi_deltaR" ]              = (ROOT.RDF.TH1DModel("phiPi_deltaR"                      , '',   21,      0,     1), r"#Delta R(#phi,#pi)"                         ,    0)
#models["dsMu_deltaR" ]               = (ROOT.RDF.TH1DModel("dsMu_deltaR"                       , '',   21,      0,   1.3), r"#Delta R(D_{s},#mu)"                        ,    0)
#####
###### kinematics
#####
##models["bs_boost_coll" ]             = (ROOT.RDF.TH1DModel("bs_boost_coll"                     , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
#####models["bs_boost_lhcb" ]             = (ROOT.RDF.TH1DModel("bs_boost_lhcb"                     , '',   21,      0,     1), r"Boost(B_{s})     "                         ,    0)
####models["bs_boost_lhcb_alt" ]         = (ROOT.RDF.TH1DModel("bs_boost_lhcb_alt"                 , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
#####models["bs_boost_reco_1" ]           = (ROOT.RDF.TH1DModel("bs_boost_reco_1"                   , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
#####models["bs_boost_reco_2" ]           = (ROOT.RDF.TH1DModel("bs_boost_reco_2"                   , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
####models["bs_boost_reco_weighted" ]    = (ROOT.RDF.TH1DModel("bs_boost_reco_weighted"            , '',   21,      0,     1), r"Boost(B_{s})      "                         ,    0)
####
####
##### Angles
##models["cosMuW_coll" ]               = (ROOT.RDF.TH1DModel("cosMuWColl"                        , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{coll}"                           ,    0)
#####models["cosMuW_lhcb" ]               = (ROOT.RDF.TH1DModel("cosMuWLhcb"                        , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
##models["cosMuW_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("cosMuWLhcbAlt"                     , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{xyz}"                            ,    0)
##models["cosMuW_reco_1" ]             = (ROOT.RDF.TH1DModel("cosMuWReco1"                       , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{math,1}"                         ,    0)
##models["cosMuW_reco_2" ]             = (ROOT.RDF.TH1DModel("cosMuWReco2"                       , '',   21,     -1,     1), r"cos(#theta^{W,\mu}_{hel})_{math,2}"                         ,    0)
#####models["cosMuW_reco_weighted" ]      = (ROOT.RDF.TH1DModel("cosMuW_reco_weighted"              , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
#####
##models["cosPhiDs_coll" ]             = (ROOT.RDF.TH1DModel("cosPhiDsColl"                      , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
###models["cosPhiDs_lhcb" ]             = (ROOT.RDF.TH1DModel("cosPhiDsLhcb"                      , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
##models["cosPhiDs_lhcb_alt" ]         = (ROOT.RDF.TH1DModel("cosPhiDsLhcbAlt"                   , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
##models["cosPhiDs_reco_1" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco1"                     , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
##models["cosPhiDs_reco_2" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco2"                     , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
#####
##models["cosPlaneBs_coll" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsColl"                    , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
###models["cosPlaneBs_lhcb" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsLhcb"                    , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
##models["cosPlaneBs_lhcb_alt" ]       = (ROOT.RDF.TH1DModel("cosPlaneBsLhcbAlt"                 , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
##models["cosPlaneBs_reco_1" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco1"                   , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
##models["cosPlaneBs_reco_2" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco2"                   , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
#####
##models["lxy_ds_sig" ]                = (ROOT.RDF.TH1DModel("lxy_ds_sig"                        , '',   21,      0,    10), r"L^{sig}_{xy}(D_{s})"                        ,    0)
##models["lxy_ds"     ]                = (ROOT.RDF.TH1DModel("lxy_ds"                            , '',   21,      0,     2), r"L_{xy}(D_{s})"                        ,    0)
###models["dxy_mu_sig" ]                = (ROOT.RDF.TH1DModel("dxy_mu_sig"                        , '',   31,    -40,    40), r"D^{sig}_{xy}"                               ,    0)
##models["cosPiK1" ]                   = (ROOT.RDF.TH1DModel("cosPiK1"                           , '',   21,     -1,     1), r"cos(#theta_{4})"                            ,    0)
####
#####Vertices
####
##models["pv_prob" ]                   = (ROOT.RDF.TH1DModel("pv_prob"                           , '',   30,    0.0,  1.0), r"Prob.(PV)"                                   ,    0)
##models["sv_prob" ]                   = (ROOT.RDF.TH1DModel("sv_prob"                           , '',   30,    0.0,  1.0), r"Prob.(B_{s} Vtx.)"                           ,    0)
##models["tv_prob" ]                   = (ROOT.RDF.TH1DModel("tv_prob"                           , '',   30,    0.0,  1.0), r"Prob.(D_{s} Vtx.)"                           ,    0)
##models["fv_prob" ]                   = (ROOT.RDF.TH1DModel("fv_prob"                           , '',   30,    0.0,  1.0), r"Prob.(#Phi Vtx.)"                            ,    0)
#
#
##models["rel_iso_03" ]                      = (ROOT.RDF.TH1DModel("rel_iso_03"                              , '',   30,    0.0,  1.0), r"Iso^{rel}_#mu"                      ,    0)
##models["rel_iso_03_pv" ]                   = (ROOT.RDF.TH1DModel("rel_iso_03_pv"                           , '',   30,    0.0,  1.0), r"Iso^{rel}_#mu(PV)"                  ,    0)
##models["rel_iso_03_sv" ]                   = (ROOT.RDF.TH1DModel("rel_iso_03_sv"                           , '',   30,    0.0,  1.0), r"Iso^{rel}_#mu(SV)"                  ,    0)
##models["rel_iso_03_tv" ]                   = (ROOT.RDF.TH1DModel("rel_iso_03_tv"                           , '',   30,    0.0,  1.0), r"Iso^{rel}_#mu(TV)"                  ,    0)
##models["rel_iso_03_ds_sv" ]                = (ROOT.RDF.TH1DModel("rel_iso_03_ds_sv"                        , '',   30,    0.0,  1.0), r"Iso^{rel}_Ds(SV)"                   ,    0)
##
##models["dsPhoton_m" ]                      = (ROOT.RDF.TH1DModel("dsPhoton_m"                              , '',   30,    2.0,  2.3), r"m(Ds + #gamma)"                     ,    0)
##models["photon_pt" ]                      = (ROOT.RDF.TH1DModel("photon_pt"                              , '',   30,    0.0,  1), r"p_{T}(#gamma) (GeV)"                     ,    0)
##models["ds_vtx_cosine_xy" ]                = (ROOT.RDF.TH1DModel("ds_vtx_cosine_xy"                        , '',   30,    -1.0, 1.0), r"Ds Vtx Cos"                         ,    0)
##models["ds_vtx_cosine_xy_pv" ]             = (ROOT.RDF.TH1DModel("ds_vtx_cosine_xy_pv"                     , '',   30,    -1.0, 1.0), r"Ds Vtx Cos PV"                      ,    0)
##models["rel_iso_03_pv" ]             = (ROOT.RDF.TH1DModel("rel_iso_03_pv"                     , '',   30,    0.0, 1.0), r"Ds Vtx Cos PV"                      ,    0)
##models["signed_decay_ip3d_mu_ds_sig_sv" ]             = (ROOT.RDF.TH1DModel("signed_decay_ip3d_mu_ds_sv"                     , '',   30,    -0.5, 0.5), r"IP3D"                      ,    0)
##models["rel_iso_03_pv" ]             = (ROOT.RDF.TH1DModel("rel_iso_03_pv"                     , '',   30,    0.0, 0.2), r"Iso"                      ,    0)
#
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



pastNN_models["class" ]                    = (ROOT.RDF.TH1DModel("class"                            , '',   6,       0,    5.9999), r"Class"                            ,    0)

#pastNN_models["score0" ]                   = (ROOT.RDF.TH1DModel("score0"                           , '',   31,       0,    1), r"Score 0"                            ,    0)
pastNN_models["score1" ]                   = (ROOT.RDF.TH1DModel("score1"                           , '',   31,       0,    1), r"Score 1"                            ,    0)
#pastNN_models["score2" ]                   = (ROOT.RDF.TH1DModel("score2"                           , '',   31,       0,    1), r"Score 2"                            ,    0)
#pastNN_models["score3" ]                   = (ROOT.RDF.TH1DModel("score3"                           , '',   31,       0,    1), r"Score 3"                            ,    0)
#pastNN_models["score4" ]                   = (ROOT.RDF.TH1DModel("score4"                           , '',   31,       0,    1), r"Score 4"                            ,    0)
#pastNN_models["score5" ]                   = (ROOT.RDF.TH1DModel("score5"                           , '',   31,       0,    1), r"Score 5"                            ,    0)
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


























