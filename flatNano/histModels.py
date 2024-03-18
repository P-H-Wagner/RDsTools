import ROOT

##############################################################################
## Define the structure of your histograms which are plotted in plotNTuples ##
##############################################################################



models = {}

##                                  (ROOT.RDF.TH1DModel("NAME"                             , '', BINS,  MIN,   MAX), "LABEL",               0)

# Basics
models["q2Gen" ]                    = (ROOT.RDF.TH1DModel("q2Gen"                           , '', 30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2Coll" ]                   = (ROOT.RDF.TH1DModel("q2Coll"                          , '', 30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2Lhcb" ]                   = (ROOT.RDF.TH1DModel("q2Lhcb"                          , '', 30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2LhcbAlt" ]                = (ROOT.RDF.TH1DModel("q2LhcbAlt"                       , '', 30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2Reco1" ]                  = (ROOT.RDF.TH1DModel("q2Reco1"                         , '', 30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)
models["q2Reco2" ]                  = (ROOT.RDF.TH1DModel("q2Reco2"                         , '', 30,      0,    12), r"q^{2} (GeV^{2})"                            ,    0)

models["m2missGen" ]                = (ROOT.RDF.TH1DModel("m2missGen"                       , '', 30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)
models["m2missColl" ]               = (ROOT.RDF.TH1DModel("m2missColl"                      , '', 30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)
models["m2missLhcb" ]               = (ROOT.RDF.TH1DModel("m2missLhcb"                      , '', 30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)
models["m2missLhcbAlt" ]            = (ROOT.RDF.TH1DModel("m2missLhcbAlt"                   , '', 30,      0,    8), r"m^{2}_{miss} (GeV^{2})"                      ,    0)

models["bs_gen_pt" ]                = (ROOT.RDF.TH1DModel("bs_gen_pt"                       , '', 30,      0,    50), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_coll" ]               = (ROOT.RDF.TH1DModel("bs_coll_pt"                      , '', 30,      0,    50), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_lhcb" ]               = (ROOT.RDF.TH1DModel("bs_pt_lhcb"                      , '', 30,      0,    50), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_lhcb_alt" ]           = (ROOT.RDF.TH1DModel("bs_pt_lhcb_alt"                  , '', 30,      0,    50), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_reco_1" ]             = (ROOT.RDF.TH1DModel("bs_pt_reco_1"                    , '', 30,      0,    50), r"p_{T}(B_{s}) (GeV)"                         ,    0)
models["bs_pt_reco_2" ]             = (ROOT.RDF.TH1DModel("bs_pt_reco2"                     , '', 30,      0,    50), r"p_{T}(B_{s}) (GeV)"                         ,    0)


# Angles
models["cosMuWGen" ]                = (ROOT.RDF.TH1DModel("cosMuWGen"                       , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuWColl" ]               = (ROOT.RDF.TH1DModel("cosMuWColl"                      , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuWLhcb" ]               = (ROOT.RDF.TH1DModel("cosMuWLhcb"                      , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuWLhcbAlt" ]            = (ROOT.RDF.TH1DModel("cosMuWLhcbAlt"                   , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuWReco1" ]              = (ROOT.RDF.TH1DModel("cosMuWReco1"                     , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)
models["cosMuWReco2" ]              = (ROOT.RDF.TH1DModel("cosMuWReco2"                     , '', 30,      -1,    1), r"cos(#theta_{1})"                            ,    0)

models["cosPhiDsGen" ]              = (ROOT.RDF.TH1DModel("cosPhiDsGen"                     , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDsColl" ]             = (ROOT.RDF.TH1DModel("cosPhiDsColl"                    , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDsLhcb" ]             = (ROOT.RDF.TH1DModel("cosPhiDsLhcb"                    , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDsLhcbAlt" ]          = (ROOT.RDF.TH1DModel("cosPhiDsLhcbAlt"                 , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDsReco1" ]            = (ROOT.RDF.TH1DModel("cosPhiDsReco1"                   , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)
models["cosPhiDsReco2" ]            = (ROOT.RDF.TH1DModel("cosPhiDsReco2"                   , '', 30,      -1,    1), r"cos(#theta_{2})"                            ,    0)

models["cosPlaneBsGen" ]            = (ROOT.RDF.TH1DModel("cosPlaneBsGen"                   , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBsColl" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsColl"                  , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBsLhcb" ]           = (ROOT.RDF.TH1DModel("cosPlaneBsLhcb"                  , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBsLhcbAlt" ]        = (ROOT.RDF.TH1DModel("cosPlaneBsLhcbAlt"               , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBsReco1" ]          = (ROOT.RDF.TH1DModel("cosPlaneBsReco1"                 , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)
models["cosPlaneBsReco2" ]          = (ROOT.RDF.TH1DModel("cosPlaneBsReco2"                 , '', 30,      -1,    1), r"cos(#chi)"                                  ,    0)

# Vertices
models["pv_x_gen" ]                 = (ROOT.RDF.TH1DModel("pv_x_gen"                        , '', 30,    0.0,  0.02), r"(PV)_{x} (cm)"                            ,    0)
models["pv_y_gen" ]                 = (ROOT.RDF.TH1DModel("pv_y_gen"                        , '', 30,   0.02,  0.06), r"(PV)_{y} (cm)"                            ,    0)
models["pv_z_gen" ]                 = (ROOT.RDF.TH1DModel("pv_z_gen"                        , '', 30,     -5,     5), r"(PV)_{z} (cm)"                            ,    0)
models["pv_x" ]                     = (ROOT.RDF.TH1DModel("pv_x"                            , '', 30,    0.0,  0.02), r"(PV)_{x} (cm)"                            ,    0)
models["pv_y" ]                     = (ROOT.RDF.TH1DModel("pv_y"                            , '', 30,   0.02,    12), r"(PV)_{y} (cm)"                            ,    0)
models["pv_z" ]                     = (ROOT.RDF.TH1DModel("pv_z"                            , '', 30,     -5,     5), r"(PV)_{z} (cm)"                            ,    0)

models["sv_x_gen" ]                 = (ROOT.RDF.TH1DModel("sv_x_gen"                        , '', 30,     -1,     1), r"(SV)_{x} (cm)"                            ,    0)
models["sv_y_gen" ]                 = (ROOT.RDF.TH1DModel("sv_y_gen"                        , '', 30,     -1,     1), r"(SV)_{y} (cm)"                            ,    0)
models["sv_z_gen" ]                 = (ROOT.RDF.TH1DModel("sv_z_gen"                        , '', 30,     -5,     5), r"(SV)_{z} (cm)"                            ,    0)
models["sv_x" ]                     = (ROOT.RDF.TH1DModel("sv_x"                            , '', 30,     -1,     1), r"(SV)_{x} (cm)"                            ,    0)
models["sv_y" ]                     = (ROOT.RDF.TH1DModel("sv_y"                            , '', 30,     -1,     1), r"(SV)_{y} (cm)"                            ,    0)
models["sv_z" ]                     = (ROOT.RDF.TH1DModel("sv_z"                            , '', 30,     -5,     5), r"(SV)_{z} (cm)"                            ,    0)

models["tv_x_gen" ]                 = (ROOT.RDF.TH1DModel("tv_x_gen"                        , '', 30,     -1,     1), r"(TV)_{x} (cm)"                            ,    0)
models["tv_y_gen" ]                 = (ROOT.RDF.TH1DModel("tv_y_gen"                        , '', 30,     -1,     1), r"(TV)_{y} (cm)"                            ,    0)
models["tv_z_gen" ]                 = (ROOT.RDF.TH1DModel("tv_z_gen"                        , '', 30,     -5,     5), r"(TV)_{z} (cm)"                            ,    0)
models["tv_x" ]                     = (ROOT.RDF.TH1DModel("tv_x"                            , '', 30,     -1,     1), r"(TV)_{x} (cm)"                            ,    0)
models["tv_y" ]                     = (ROOT.RDF.TH1DModel("tv_y"                            , '', 30,     -1,     1), r"(TV)_{y} (cm)"                            ,    0)
models["tv_z" ]                     = (ROOT.RDF.TH1DModel("tv_z"                            , '', 30,     -5,     5), r"(TV)_{z} (cm)"                            ,    0)

models["fv_x_gen" ]                 = (ROOT.RDF.TH1DModel("fv_x_gen"                        , '', 30,     -1,     1), r"(FV)_{x} (cm)"                            ,    0)
models["fv_y_gen" ]                 = (ROOT.RDF.TH1DModel("fv_y_gen"                        , '', 30,     -1,     1), r"(FV)_{y} (cm)"                            ,    0)
models["fv_z_gen" ]                 = (ROOT.RDF.TH1DModel("fv_z_gen"                        , '', 30,     -5,     5), r"(FV)_{z} (cm)"                            ,    0)
models["fv_x" ]                     = (ROOT.RDF.TH1DModel("fv_x"                            , '', 30,     -1,     1), r"(FV)_{x} (cm)"                            ,    0)
models["fv_y" ]                     = (ROOT.RDF.TH1DModel("fv_y"                            , '', 30,     -1,     1), r"(FV)_{y} (cm)"                            ,    0)
models["fv_z" ]                     = (ROOT.RDF.TH1DModel("fv_z"                            , '', 30,     -5,     5), r"(FV)_{z} (cm)"                            ,    0)








