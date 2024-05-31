import ROOT

##############################################################################
## Define the structure of your histograms which are plotted in plotNTuples ##
##############################################################################



models = {}
## 
models["phiPi_m" ]                   = (ROOT.RDF.TH1DModel("phiPi_m"                           , '',   21,   1.91, 2.028), r"D_{s} mass (GeV)"                           ,    0)

# Basics
#models["gen_q2" ]                    = (ROOT.RDF.TH1DModel("gen_q2"                            , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2_coll" ]                   = (ROOT.RDF.TH1DModel("q2_coll"                           , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#models["q2_lhcb" ]                  = (ROOT.RDF.TH1DModel("q2_lhcb"                           , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#models["q2_lhcb_alt" ]               = (ROOT.RDF.TH1DModel("q2_lhcb_alt"                       , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#models["q2_reco_1" ]                 = (ROOT.RDF.TH1DModel("q2_reco_1"                         , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#models["q2_reco_2" ]                 = (ROOT.RDF.TH1DModel("q2_reco_2"                         , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
#models["q2_reco_weighted" ]          = (ROOT.RDF.TH1DModel("q2_reco_weighted"                  , '',   21,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
"""
models["gen_m2_miss" ]               = (ROOT.RDF.TH1DModel("gen_m2_miss"                       , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
models["m2_miss_coll" ]              = (ROOT.RDF.TH1DModel("m2_miss_coll"                      , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
#models["m2_miss_lhcb" ]             = (ROOT.RDF.TH1DModel("m2_miss_lhcb"                      , '',   21,      0,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)
models["m2_miss_lhcb_alt" ]          = (ROOT.RDF.TH1DModel("m2_miss_lhcb_alt"                  , '',   21,    -10,     6), r"m^{2}_{miss} (GeV^{2})"                     ,    0)

models["gen_pt_miss" ]               = (ROOT.RDF.TH1DModel("gen_pt_miss"                       , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
models["pt_miss_coll" ]              = (ROOT.RDF.TH1DModel("pt_miss_coll"                      , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
#models["pt_miss_lhcb" ]             = (ROOT.RDF.TH1DModel("pt_miss_lhcb"                      , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
models["pt_miss_lhcb_alt" ]          = (ROOT.RDF.TH1DModel("pt_miss_lhcb_alt"                  , '',   21,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)

models["gen_bs_pt" ]                 = (ROOT.RDF.TH1DModel("gen_bs_pt"                         , '',   21,      0,    90), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_coll" ]                = (ROOT.RDF.TH1DModel("bs_coll_pt"                        , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
#models["bs_pt_lhcb" ]               = (ROOT.RDF.TH1DModel("bs_pt_lhcb"                        , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_lhcb_alt" ]            = (ROOT.RDF.TH1DModel("bs_pt_lhcb_alt"                    , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_reco_1" ]              = (ROOT.RDF.TH1DModel("bs_pt_reco_1"                      , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_reco_2" ]              = (ROOT.RDF.TH1DModel("bs_pt_reco_2"                      , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_reco_weighted" ]       = (ROOT.RDF.TH1DModel("bs_pt_reco_weighted"               , '',   21,      0,    80), r"p_{T}(B_{s}) (GeV)"                         ,    0)

models["gen_e_star" ]                = (ROOT.RDF.TH1DModel("gen_e_star"                        , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
models["e_star_coll" ]               = (ROOT.RDF.TH1DModel("e_star_coll"                       , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
models["e_star_lhcb" ]              = (ROOT.RDF.TH1DModel("e_star_lhcb"                       , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
models["e_star_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("e_star_lhcb_alt"                   , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
models["e_star_reco_1" ]             = (ROOT.RDF.TH1DModel("e_star_reco_1"                     , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
models["e_star_reco_2" ]             = (ROOT.RDF.TH1DModel("e_star_reco_2"                     , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)
models["e_star_reco_weighted" ]      = (ROOT.RDF.TH1DModel("e_star_reco_weighted"              , '',   21,      0,     3), r"E*_{#mu} (GeV)"                             ,    0)

models["gen_e_gamma" ]               = (ROOT.RDF.TH1DModel("gen_e_gamma"                       , '',   21,      0,     4), r"E_{#gamma} (GeV)"                           ,    0)
models["e_gamma" ]                   = (ROOT.RDF.TH1DModel("e_gamma"                           , '',   21,      0,     4), r"E_{#gamma} (GeV)"                           ,    0)
#models["gen_mu_rel_iso" ]           = (ROOT.RDF.TH1DModel("gen_mu_rel_iso"                    , '',   21,      0,    20), r"Iso^{rel}_{#mu}"                            ,    0)
models["mu_rel_iso_03" ]             = (ROOT.RDF.TH1DModel("mu_rel_iso"                        , '',   21,      0,   0.3), r"Iso^{rel}_{#mu}"                            ,    0)
#models["disc_is_negative" ]         = (ROOT.RDF.TH1DModel("disc_is_negative"                  , '',   21,     -1,     2), r"Disc < 0"                                   ,    0)
#models["disc_negativity" ]          = (ROOT.RDF.TH1DModel("disc_negativity"                   , '',   21,      0,     6), r"Disc < 0"                                   ,    0)
# Angles
models["gen_cosMuW" ]                = (ROOT.RDF.TH1DModel("gen_cosMuWGen"                     , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
#models["cosMuW_coll" ]               = (ROOT.RDF.TH1DModel("cosMuWColl"                        , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
#models["cosMuW_lhcb" ]              = (ROOT.RDF.TH1DModel("cosMuWLhcb"                        , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("cosMuWLhcbAlt"                     , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_reco_1" ]             = (ROOT.RDF.TH1DModel("cosMuWReco1"                       , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_reco_2" ]             = (ROOT.RDF.TH1DModel("cosMuWReco2"                       , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_reco_weighted" ]      = (ROOT.RDF.TH1DModel("cosMuW_reco_weighted"       , '',   21,     -1,     1), r"cos(#theta_{1})"                            ,    0)

models["gen_cosPhiDs" ]              = (ROOT.RDF.TH1DModel("gen_cosPhiDsGen"                   , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_coll" ]             = (ROOT.RDF.TH1DModel("cosPhiDsColl"                      , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
#models["cosPhiDs_lhcb" ]            = (ROOT.RDF.TH1DModel("cosPhiDsLhcb"                      , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_lhcb_alt" ]         = (ROOT.RDF.TH1DModel("cosPhiDsLhcbAlt"                   , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_reco_1" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco1"                     , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_reco_2" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco2"                     , '',   21,     -1,     1), r"cos(#theta_{2})"                            ,    0)

models["gen_cosPlaneBs" ]            = (ROOT.RDF.TH1DModel("gen_cosPlaneBsGen"                 , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_coll" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsColl"                    , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
#models["cosPlaneBs_lhcb" ]          = (ROOT.RDF.TH1DModel("cosPlaneBsLhcb"                    , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_lhcb_alt" ]       = (ROOT.RDF.TH1DModel("cosPlaneBsLhcbAlt"                 , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_reco_1" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco1"                   , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_reco_2" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco2"                   , '',   21,     -1,     1), r"cos(#chi)"                                  ,    0)
#models["gen_lxy_ds_sig" ]           = (ROOT.RDF.TH1DModel("gen_lxy_ds_sig"                    , '',   21,     -1,     2), r"cos(#chi)"                                  ,    0)
models["lxy_ds_sig" ]                = (ROOT.RDF.TH1DModel("lxy_ds_sig"                        , '',   21,      0,    10), r"L^{sig}_{xy}(D_{s})"                        ,    0)
#models["gen_dxy_mu_sig" ]           = (ROOT.RDF.TH1DModel("gen_dxy_mu_sig"                    , '',   21,    -50,    50), r"cos(#chi)"                                  ,    0)
models["dxy_mu_sig" ]                = (ROOT.RDF.TH1DModel("dxy_mu_sig"                        , '',   31,    -40,    40), r"D^{sig}_{xy}"                               ,    0)
models["cosPiK1" ]                   = (ROOT.RDF.TH1DModel("cosPiK1"                           , '',   21,     -1,     1), r"cos(#theta_{4})"                            ,    0)
"""
"""
# Vertices
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


"""
