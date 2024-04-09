import ROOT

##############################################################################
## Define the structure of your histograms which are plotted in plotNTuples ##
##############################################################################



models = {}

##                                  (ROOT.RDF.TH1DModel("NAME"                             , '', BINS,  MIN,   MAX), "LABEL",               0)

# Basics
models["gen_q2" ]                     = (ROOT.RDF.TH1DModel("gen_q2"                            , '',   30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2_coll" ]                    = (ROOT.RDF.TH1DModel("q2_coll"                           , '',   30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2_lhcb" ]                    = (ROOT.RDF.TH1DModel("q2_lhcb"                           , '',   30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2_lhcb_alt" ]                = (ROOT.RDF.TH1DModel("q2_lhcb_alt"                       , '',   30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2_reco_1" ]                  = (ROOT.RDF.TH1DModel("q2_reco_1"                         , '',   30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2_reco_2" ]                  = (ROOT.RDF.TH1DModel("q2_reco_2"                         , '',   30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)

models["gen_m2_miss" ]                = (ROOT.RDF.TH1DModel("gen_m2_miss"                       , '',   30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)
models["m2_miss_coll" ]               = (ROOT.RDF.TH1DModel("m2_miss_coll"                      , '',   30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)
models["m2_miss_lhcb" ]               = (ROOT.RDF.TH1DModel("m2_miss_lhcb"                      , '',   30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)
models["m2_miss_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("m2_miss_lhcb_alt"                  , '',   30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)

models["gen_pt_miss" ]                = (ROOT.RDF.TH1DModel("gen_pt_miss"                       , '',   30,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
models["pt_miss_coll" ]               = (ROOT.RDF.TH1DModel("pt_miss_coll"                      , '',   30,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
models["pt_miss_lhcb" ]               = (ROOT.RDF.TH1DModel("pt_miss_lhcb"                      , '',   30,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)
models["pt_miss_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("pt_miss_lhcb_alt"                  , '',   30,      0,    30), r"p_{T}^{miss} (GeV)"                         ,    0)

models["gen_bs_pt" ]                  = (ROOT.RDF.TH1DModel("gen_bs_pt"                         , '',  30,      0,    90), r"p_{T}(B_{s}) (GeV)"                          ,    0)
models["bs_pt_coll" ]                 = (ROOT.RDF.TH1DModel("bs_coll_pt"                        , '',  30,      0,    90), r"p_{T}(B_{s}) (GeV)"                          ,    0)
models["bs_pt_lhcb" ]                 = (ROOT.RDF.TH1DModel("bs_pt_lhcb"                        , '',  30,      0,    90), r"p_{T}(B_{s}) (GeV)"                          ,    0)
models["bs_pt_lhcb_alt" ]             = (ROOT.RDF.TH1DModel("bs_pt_lhcb_alt"                    , '',  30,      0,    90), r"p_{T}(B_{s}) (GeV)"                          ,     0)
models["bs_pt_reco_1" ]               = (ROOT.RDF.TH1DModel("bs_pt_reco_1"                      , '',  30,      0,    90), r"p_{T}(B_{s}) (GeV)"                          ,    0)
models["bs_pt_reco_2" ]               = (ROOT.RDF.TH1DModel("bs_pt_reco2"                       , '',  30,      0,    90), r"p_{T}(B_{s}) (GeV)"                          ,    0)

models["gen_e_star" ]                 = (ROOT.RDF.TH1DModel("gen_e_star"                        , '', 50,      0,    3), r"E*_{#mu} (GeV)"                                ,    0)
models["e_star_coll" ]                = (ROOT.RDF.TH1DModel("e_star_coll"                       , '', 50,      0,    3), r"E*_{#mu} (GeV)"                                ,    0)
models["e_star_lhcb" ]                = (ROOT.RDF.TH1DModel("e_star_lhcb"                       , '', 50,      0,    3), r"E*_{#mu} (GeV)"                                ,    0)
models["e_star_lhcb_alt" ]            = (ROOT.RDF.TH1DModel("e_star_lhcb_alt"                   , '', 50,      0,    3), r"E*_{#mu} (GeV)"                                ,    0)
models["e_star_reco_1" ]              = (ROOT.RDF.TH1DModel("e_star_reco_1"                     , '', 50,      0,    3), r"E*_{#mu} (GeV)"                                ,    0)
models["e_star_reco_2" ]              = (ROOT.RDF.TH1DModel("e_star_reco_2"                     , '', 50,      0,    3), r"E*_{#mu} (GeV)"                                ,    0)

models["gen_e_gamma" ]                = (ROOT.RDF.TH1DModel("gen_e_gamma"                       , '', 50,      0,    7), r"E_{#gamma} (GeV)"                             ,    0)
models["e_gamma" ]                    = (ROOT.RDF.TH1DModel("e_gamma"                           , '', 50,      0,    7), r"E_{#gamma} (GeV)"                             ,    0)

models["gen_mu_rel_iso" ]             = (ROOT.RDF.TH1DModel("gen_mu_rel_iso"                    , '', 50,      0,    20), r"Iso^{rel}_{#mu}"                              ,    0)
models["mu_rel_iso" ]                 = (ROOT.RDF.TH1DModel("mu_rel_iso"                        , '', 50,      0,    20), r"Iso^{rel}_{#mu}"                              ,    0)

models["disc_is_negative" ]          = (ROOT.RDF.TH1DModel("disc_is_negative"                    , '', 5,      -1,   2), r"Disc < 0"                              ,    0)
#models["gen_disc_is_ngeative" ]      = (ROOT.RDF.TH1DModel("gen_disc_is_negative"                    , '', 5,      -1,   2), r"Disc < 0"                              ,    0)

# Angles
models["gen_cosMuW" ]                = (ROOT.RDF.TH1DModel("gen_cosMuWGen"                       , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_coll" ]               = (ROOT.RDF.TH1DModel("cosMuWColl"                      , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_lhcb" ]               = (ROOT.RDF.TH1DModel("cosMuWLhcb"                      , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("cosMuWLhcbAlt"                   , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_reco_1" ]             = (ROOT.RDF.TH1DModel("cosMuWReco1"                     , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuW_reco_2" ]             = (ROOT.RDF.TH1DModel("cosMuWReco2"                     , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)

models["gen_cosPhiDs" ]              = (ROOT.RDF.TH1DModel("gen_cosPhiDsGen"                     , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_coll" ]             = (ROOT.RDF.TH1DModel("cosPhiDsColl"                    , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_lhcb" ]             = (ROOT.RDF.TH1DModel("cosPhiDsLhcb"                    , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_lhcb_alt" ]         = (ROOT.RDF.TH1DModel("cosPhiDsLhcbAlt"                 , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_reco_1" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco1"                   , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDs_reco_2" ]           = (ROOT.RDF.TH1DModel("cosPhiDsReco2"                   , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)

models["gen_cosPlaneBs" ]            = (ROOT.RDF.TH1DModel("gen_cosPlaneBsGen"                   , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_coll" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsColl"                  , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_lhcb" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsLhcb"                  , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_lhcb_alt" ]       = (ROOT.RDF.TH1DModel("cosPlaneBsLhcbAlt"               , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_reco_1" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco1"                 , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBs_reco_2" ]         = (ROOT.RDF.TH1DModel("cosPlaneBsReco2"                 , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
"""
# Vertices
models["pv_x_gen" ]                  = (ROOT.RDF.TH1DModel("pv_x_gen"                        , '', 30,    0.0,  0.02), r"(PV)_{x} (cm)"                            ,    0)
models["pv_y_gen" ]                  = (ROOT.RDF.TH1DModel("pv_y_gen"                        , '', 30,   0.02,  0.06), r"(PV)_{y} (cm)"                            ,    0)
models["pv_z_gen" ]                  = (ROOT.RDF.TH1DModel("pv_z_gen"                        , '', 30,     -5,     5), r"(PV)_{z} (cm)"                            ,    0)
models["pv_x" ]                      = (ROOT.RDF.TH1DModel("pv_x"                            , '', 30,    0.0,  0.02), r"(PV)_{x} (cm)"                            ,    0)
models["pv_y" ]                      = (ROOT.RDF.TH1DModel("pv_y"                            , '', 30,   0.02,    12), r"(PV)_{y} (cm)"                            ,    0)
models["pv_z" ]                      = (ROOT.RDF.TH1DModel("pv_z"                            , '', 30,     -5,     5), r"(PV)_{z} (cm)"                            ,    0)

models["sv_x_gen" ]                  = (ROOT.RDF.TH1DModel("sv_x_gen"                        , '', 30,     -1,     1), r"(SV)_{x} (cm)"                            ,    0)
models["sv_y_gen" ]                  = (ROOT.RDF.TH1DModel("sv_y_gen"                        , '', 30,     -1,     1), r"(SV)_{y} (cm)"                            ,    0)
models["sv_z_gen" ]                  = (ROOT.RDF.TH1DModel("sv_z_gen"                        , '', 30,     -5,     5), r"(SV)_{z} (cm)"                            ,    0)
models["sv_x" ]                      = (ROOT.RDF.TH1DModel("sv_x"                            , '', 30,     -1,     1), r"(SV)_{x} (cm)"                            ,    0)
models["sv_y" ]                      = (ROOT.RDF.TH1DModel("sv_y"                            , '', 30,     -1,     1), r"(SV)_{y} (cm)"                            ,    0)
models["sv_z" ]                      = (ROOT.RDF.TH1DModel("sv_z"                            , '', 30,     -5,     5), r"(SV)_{z} (cm)"                            ,    0)

models["tv_x_gen" ]                  = (ROOT.RDF.TH1DModel("tv_x_gen"                        , '', 30,     -1,     1), r"(TV)_{x} (cm)"                            ,    0)
models["tv_y_gen" ]                  = (ROOT.RDF.TH1DModel("tv_y_gen"                        , '', 30,     -1,     1), r"(TV)_{y} (cm)"                            ,    0)
models["tv_z_gen" ]                  = (ROOT.RDF.TH1DModel("tv_z_gen"                        , '', 30,     -5,     5), r"(TV)_{z} (cm)"                            ,    0)
models["tv_x" ]                      = (ROOT.RDF.TH1DModel("tv_x"                            , '', 30,     -1,     1), r"(TV)_{x} (cm)"                            ,    0)
models["tv_y" ]                      = (ROOT.RDF.TH1DModel("tv_y"                            , '', 30,     -1,     1), r"(TV)_{y} (cm)"                            ,    0)
models["tv_z" ]                      = (ROOT.RDF.TH1DModel("tv_z"                            , '', 30,     -5,     5), r"(TV)_{z} (cm)"                            ,    0)
 
models["fv_x_gen" ]                  = (ROOT.RDF.TH1DModel("fv_x_gen"                        , '', 30,     -1,     1), r"(FV)_{x} (cm)"                            ,    0)
models["fv_y_gen" ]                  = (ROOT.RDF.TH1DModel("fv_y_gen"                        , '', 30,     -1,     1), r"(FV)_{y} (cm)"                            ,    0)
models["fv_z_gen" ]                  = (ROOT.RDF.TH1DModel("fv_z_gen"                        , '', 30,     -5,     5), r"(FV)_{z} (cm)"                            ,    0)
models["fv_x" ]                      = (ROOT.RDF.TH1DModel("fv_x"                            , '', 30,     -1,     1), r"(FV)_{x} (cm)"                            ,    0)
models["fv_y" ]                      = (ROOT.RDF.TH1DModel("fv_y"                            , '', 30,     -1,     1), r"(FV)_{y} (cm)"                            ,    0)
models["fv_z" ]                      = (ROOT.RDF.TH1DModel("fv_z"                            , '', 30,     -5,     5), r"(FV)_{z} (cm)"                            ,    0)



"""




