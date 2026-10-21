import json
import json
import ROOT
import os
import sys
import numpy as np

selec= "HOOK_SELECTION"

# skim it!

#create rdf
chain = ROOT.TChain("tree")
chain.Add("HOOK_FILE_IN")

#output
destination = "HOOK_FILE_OUT"
print("saving to:", destination)
print("====> Create rdf")
df = ROOT.RDataFrame(chain)
print("====> rdf DONE")

#lead python to the correct header folder
ROOT.gInterpreter.AddIncludePath(
    "/work/pahwagne/environments/tf/lib/python3.10/site-packages/correctionlib/include"
)

ROOT.gSystem.Load(
    "/work/pahwagne/environments/tf/lib/python3.10/site-packages/correctionlib/lib/libcorrectionlib.so"
)

################################# 
# PU correction + sys.          #
################################# 

pu_weights           = np.loadtxt("/work/pahwagne/RDsTools/pileup/pu_weights.txt")
pu_weights_str       = ", ".join(map(str, pu_weights   )) # of type: '0.01, 0.34, 0.89'

pu_weights_up        = np.loadtxt("/work/pahwagne/RDsTools/pileup/pu_weights_up.txt")
pu_weights_up_str    = ", ".join(map(str, pu_weights_up   )) # of type: '0.01, 0.34, 0.89'

pu_weights_down      = np.loadtxt("/work/pahwagne/RDsTools/pileup/pu_weights_down.txt")
pu_weights_down_str  = ", ".join(map(str, pu_weights_down   )) # of type: '0.01, 0.34, 0.89'


ROOT.gInterpreter.Declare(rf'''double get_w_pu(int npv){{

    float weights_cpp[] = {{{pu_weights_str}}};
    return weights_cpp[npv];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_w_pu_up(int npv){{

    float weights_cpp[] = {{{pu_weights_up_str}}};
    return weights_cpp[npv];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_w_pu_down(int npv){{

    float weights_cpp[] = {{{pu_weights_down_str}}};
    return weights_cpp[npv];
    }}
    '''
    )


################################# 
# Trigger SF correction + sys.  #
################################# 

# NEW scale factors obtained with  morphing
scale_factors        = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_fullBPark_morphed_with_correct_l1_prescales.txt")
scale_factors_err_hi = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_hi_fullBPark_morphed_with_correct_l1_prescales.txt")
scale_factors_err_lo = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_lo_fullBPark_morphed_with_correct_l1_prescales.txt")
pt_edges             = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/pt_edges_fullBPark_morphed_with_correct_l1_prescales.txt")
dxy_edges            = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/dxy_edges_fullBPark_morphed_with_correct_l1_prescales.txt")

# old sclae factors 
scale_factors_old    = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_fullBPark_morphed.txt")
scale_factors_err_hi_old = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_hi_fullBPark_morphed.txt")
scale_factors_err_lo_old = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_lo_fullBPark_morphed.txt")
pt_edges_old             = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/pt_edges_fullBPark_morphed.txt")
dxy_edges_old            = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/dxy_edges_fullBPark_morphed.txt")

#scale_factors_old        = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_fullBPark.txt")
#scale_factors_err_hi_old = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_hi_fullBPark.txt")
#scale_factors_err_lo_old = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_lo_fullBPark.txt")
#pt_edges_old             = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/pt_edges_fullBPark.txt")
#dxy_edges_old            = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/dxy_edges_fullBPark.txt")

#scale_factors        = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_D_only.txt")
#scale_factors_err_hi = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_hi_D_only.txt")
#scale_factors_err_lo = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_lo_D_only.txt")
#scale_factors_err    = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/AM/trigger_sf_err_D_only.txt")
#pt_edges             = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/pt_edges_D_only.txt")
#dxy_edges            = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/dxy_edges_D_only.txt")


#debug
#scale_factors     = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/trigger_sf_D1.txt")
#scale_factors_err = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/trigger_sf_err_D1.txt")
#pt_edges          = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/pt_edges_D1.txt")
#dxy_edges         = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/dxy_edges_D1.txt")

#scale_factors     = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/trigger_sf_A1.txt")
#scale_factors_err = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/trigger_sf_err_A1.txt")
#pt_edges          = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/pt_edges_A1.txt")
#dxy_edges         = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/dxy_edges_A1.txt")


#convert into string such that it can be handed to the cpp process line
pt_edges_str  = ", ".join(map(str, pt_edges     )) # of type: '0.01, 0.34, 0.89'
dxy_edges_str = ", ".join(map(str, dxy_edges    )) # of type: '0.01, 0.34, 0.89'

#scale_factors is a matrix!
lenx = scale_factors.shape[0]
leny = scale_factors.shape[1]

scale_factors_str        = "{"
scale_factors_err_hi_str = "{"
scale_factors_err_lo_str = "{"
#scale_factors_err_str    = "{"

############################################
for i,row in enumerate(scale_factors):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  scale_factors_str  += sub

  if i != (len(scale_factors)-1): scale_factors_str += ","

scale_factors_str += "}"

############################################
for i,row in enumerate(scale_factors_err_hi):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  scale_factors_err_hi_str  += sub

  if i != (len(scale_factors_err_hi)-1): scale_factors_err_hi_str += ","

scale_factors_err_hi_str += "}"

############################################
for i,row in enumerate(scale_factors_err_lo):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  scale_factors_err_lo_str  += sub

  if i != (len(scale_factors_err_lo)-1): scale_factors_err_lo_str += ","

scale_factors_err_lo_str += "}"


print(scale_factors_str)
print(pt_edges_str)
print(dxy_edges_str)

#################################

pt_edges_str_old  = ", ".join(map(str, pt_edges_old     )) # of type: '0.01, 0.34, 0.89'
dxy_edges_str_old = ", ".join(map(str, dxy_edges_old    )) # of type: '0.01, 0.34, 0.89'

#scale_factors is a matrix!
lenx = scale_factors_old.shape[0]
leny = scale_factors_old.shape[1]

scale_factors_str_old        = "{"
scale_factors_err_hi_str_old = "{"
scale_factors_err_lo_str_old = "{"
#scale_factors_err_str    = "{"

############################################
for i,row in enumerate(scale_factors_old):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  scale_factors_str_old  += sub

  if i != (len(scale_factors_old)-1): scale_factors_str_old += ","

scale_factors_str_old += "}"

############################################
for i,row in enumerate(scale_factors_err_hi_old):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  scale_factors_err_hi_str_old  += sub

  if i != (len(scale_factors_err_hi_old)-1): scale_factors_err_hi_str_old += ","

scale_factors_err_hi_str_old += "}"

############################################
for i,row in enumerate(scale_factors_err_lo_old):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  scale_factors_err_lo_str_old  += sub

  if i != (len(scale_factors_err_lo_old)-1): scale_factors_err_lo_str_old += ","

scale_factors_err_lo_str_old += "}"


print(scale_factors_str_old)
print(pt_edges_str_old)
print(dxy_edges_str_old)

ROOT.gInterpreter.Declare(rf'''double get_trigger_sf(float mu_pt, float dxy_mu_sig_pv){{

    float scale_factors_cpp[{lenx}][{leny}] = {scale_factors_str};
    float pt_edges_cpp [{lenx + 1}] = {{{pt_edges_str}}};
    float dxy_edges_cpp[{leny + 1}] = {{{dxy_edges_str}}};

    int index_pt  = 0;
    int index_dxy = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((mu_pt > pt_edges_cpp[i]) && (mu_pt < pt_edges_cpp[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    float dxy_mu_sig_pv_abs = dxy_mu_sig_pv;
    if (dxy_mu_sig_pv_abs < 0){{dxy_mu_sig_pv_abs = -1 * dxy_mu_sig_pv_abs ;}}

    if ((dxy_mu_sig_pv_abs > dxy_edges_cpp[i]) && (dxy_mu_sig_pv_abs < dxy_edges_cpp[i+1])){{
        index_dxy = i;
        break;
      }}
    }}


    return scale_factors_cpp[index_pt][index_dxy];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_trigger_sf_err_hi(float mu_pt, float dxy_mu_sig_pv){{

    float scale_factors_err_hi_cpp[{lenx}][{leny}] = {scale_factors_err_hi_str};
    float pt_edges_cpp [{lenx + 1}] = {{{pt_edges_str}}};
    float dxy_edges_cpp[{leny + 1}] = {{{dxy_edges_str}}};

    int index_pt  = 0;
    int index_dxy = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((mu_pt > pt_edges_cpp[i]) && (mu_pt < pt_edges_cpp[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    float dxy_mu_sig_pv_abs = dxy_mu_sig_pv;
    if (dxy_mu_sig_pv_abs < 0){{dxy_mu_sig_pv_abs = -1 * dxy_mu_sig_pv_abs ;}}

    if ((dxy_mu_sig_pv_abs > dxy_edges_cpp[i]) && (dxy_mu_sig_pv_abs < dxy_edges_cpp[i+1])){{
        index_dxy = i;
        break;
      }}
    }}


    return scale_factors_err_hi_cpp[index_pt][index_dxy];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_trigger_sf_err_lo(float mu_pt, float dxy_mu_sig_pv){{

    float scale_factors_err_lo_cpp[{lenx}][{leny}] = {scale_factors_err_lo_str};
    float pt_edges_cpp [{lenx + 1}] = {{{pt_edges_str}}};
    float dxy_edges_cpp[{leny + 1}] = {{{dxy_edges_str}}};

    int index_pt  = 0;
    int index_dxy = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((mu_pt > pt_edges_cpp[i]) && (mu_pt < pt_edges_cpp[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    float dxy_mu_sig_pv_abs = dxy_mu_sig_pv;
    if (dxy_mu_sig_pv_abs < 0){{dxy_mu_sig_pv_abs = -1 * dxy_mu_sig_pv_abs ;}}

    if ((dxy_mu_sig_pv_abs > dxy_edges_cpp[i]) && (dxy_mu_sig_pv_abs < dxy_edges_cpp[i+1])){{
        index_dxy = i;
        break;
      }}
    }}


    return scale_factors_err_lo_cpp[index_pt][index_dxy];
    }}
    '''
    )


ROOT.gInterpreter.Declare(rf'''double get_trigger_sf_old(float mu_pt, float dxy_mu_sig_pv){{

    float scale_factors_cpp_old[{lenx}][{leny}] = {scale_factors_str_old};
    float pt_edges_cpp_old [{lenx + 1}] = {{{pt_edges_str_old}}};
    float dxy_edges_cpp_old[{leny + 1}] = {{{dxy_edges_str_old}}};

    int index_pt  = 0;
    int index_dxy = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((mu_pt > pt_edges_cpp_old[i]) && (mu_pt < pt_edges_cpp_old[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    float dxy_mu_sig_pv_abs = dxy_mu_sig_pv;
    if (dxy_mu_sig_pv_abs < 0){{dxy_mu_sig_pv_abs = -1 * dxy_mu_sig_pv_abs ;}}

    if ((dxy_mu_sig_pv_abs > dxy_edges_cpp_old[i]) && (dxy_mu_sig_pv_abs < dxy_edges_cpp_old[i+1])){{
        index_dxy = i;
        break;
      }}
    }}


    return scale_factors_cpp_old[index_pt][index_dxy];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_trigger_sf_err_hi_old(float mu_pt, float dxy_mu_sig_pv){{

    float scale_factors_err_hi_cpp_old[{lenx}][{leny}] = {scale_factors_err_hi_str_old};
    float pt_edges_cpp_old[{lenx + 1}] = {{{pt_edges_str_old}}};
    float dxy_edges_cpp_old[{leny + 1}] = {{{dxy_edges_str_old}}};

    int index_pt  = 0;
    int index_dxy = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((mu_pt > pt_edges_cpp_old[i]) && (mu_pt < pt_edges_cpp_old[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    float dxy_mu_sig_pv_abs = dxy_mu_sig_pv;
    if (dxy_mu_sig_pv_abs < 0){{dxy_mu_sig_pv_abs = -1 * dxy_mu_sig_pv_abs ;}}

    if ((dxy_mu_sig_pv_abs > dxy_edges_cpp_old[i]) && (dxy_mu_sig_pv_abs < dxy_edges_cpp_old[i+1])){{
        index_dxy = i;
        break;
      }}
    }}


    return scale_factors_err_hi_cpp_old[index_pt][index_dxy];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_trigger_sf_err_lo_old(float mu_pt, float dxy_mu_sig_pv){{

    float scale_factors_err_lo_cpp_old[{lenx}][{leny}] = {scale_factors_err_lo_str_old};
    float pt_edges_cpp_old [{lenx + 1}] = {{{pt_edges_str_old}}};
    float dxy_edges_cpp_old[{leny + 1}] = {{{dxy_edges_str_old}}};

    int index_pt  = 0;
    int index_dxy = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((mu_pt > pt_edges_cpp_old[i]) && (mu_pt < pt_edges_cpp_old[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    float dxy_mu_sig_pv_abs = dxy_mu_sig_pv;
    if (dxy_mu_sig_pv_abs < 0){{dxy_mu_sig_pv_abs = -1 * dxy_mu_sig_pv_abs ;}}

    if ((dxy_mu_sig_pv_abs > dxy_edges_cpp_old[i]) && (dxy_mu_sig_pv_abs < dxy_edges_cpp_old[i+1])){{
        index_dxy = i;
        break;
      }}
    }}


    return scale_factors_err_lo_cpp_old[index_pt][index_dxy];
    }}
    '''
    )


#ROOT.gInterpreter.Declare(rf'''double get_trigger_sf_err(float mu_pt, float dxy_mu_sig_pv){{
#
#    float scale_factors_err_cpp[{lenx}][{leny}] = {scale_factors_err_str};
#    float pt_edges_cpp [{lenx + 1}] = {{{pt_edges_str}}};
#    float dxy_edges_cpp[{leny + 1}] = {{{dxy_edges_str}}};
#
#    int index_pt  = 0;
#    int index_dxy = 0;
#
#    for(int i=0; i<{lenx}; ++i){{
#
#    if ((mu_pt > pt_edges_cpp[i]) && (mu_pt < pt_edges_cpp[i+1])){{
#        index_pt = i;
#        break;
#      }} 
#    }}
#
#    for(int i=0; i<{leny}; ++i){{
#
#    float dxy_mu_sig_pv_abs = dxy_mu_sig_pv;
#    if (dxy_mu_sig_pv_abs < 0){{dxy_mu_sig_pv_abs = -1 * dxy_mu_sig_pv_abs ;}}
#
#    if ((dxy_mu_sig_pv_abs > dxy_edges_cpp[i]) && (dxy_mu_sig_pv_abs < dxy_edges_cpp[i+1])){{
#        index_dxy = i;
#        break;
#      }}
#    }}
#
#
#    return scale_factors_err_cpp[index_pt][index_dxy];
#    }}
#    '''
#    )



################################# 
# Bs lifetime correction + sys. #
################################# 

# load Bs MC lifetime
with open("/work/pahwagne/RDsTools/corrections/bs_lifetime/lifetimeMC.json") as f:
    tau_mc = json.load(f)
    print(f"====> MC lifetime tau is: {tau_mc}")

#PDG lifetime of Bs
tau_pdg      = 1.526e-12
tau_pdg_up   = (1.526 + 0.015) * 1e-12
tau_pdg_down = (1.526 - 0.015) * 1e-12

ROOT.gInterpreter.Declare(r'''double get_time(float pv_x, float pv_y, float pv_z, float sv_x, float sv_y, float sv_z, float beta){

    //remark sv and pv quantities are in cm
    float distance = sqrt(pow(sv_x - pv_x,2) + pow(sv_y - pv_y,2) + pow(sv_z - pv_z,2));
  
    //convert distance to meters
    distance *= 0.01; 
  
    //remark boost is in atural units, i.e. a number in [0,1] (the beta)
    float gamma = 1.0 / sqrt(1 - pow(beta,2));
  
    // use c in SI units to get a time in seconds
    float c = 299792458.0;
    float t = distance / (gamma * beta * c);
  
    return t;

    }
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_w_bs_tau(double time){{

    float w = (1 / {tau_pdg}) * exp(-time/{tau_pdg}) / ( (1/{tau_mc}) * exp(-time/{tau_mc}) );
    return w;
    }}
    '''
    )
ROOT.gInterpreter.Declare(rf'''double get_w_bs_tau_up(double time){{
    
    float w = (1 / {tau_pdg_up}) * exp(-time/{tau_pdg_up}) / ( (1/{tau_mc}) * exp(-time/{tau_mc}));
    return w;
    }}
    '''
    )
ROOT.gInterpreter.Declare(rf'''double get_w_bs_tau_down(double time){{

    float w = (1 / {tau_pdg_down}) * exp(-time/{tau_pdg_down}) / ( (1/{tau_mc}) * exp(-time/{tau_mc})); 
    return w;
    }}
    '''
    )

#################
# Bs kinematics #
#################

#load the weights and edges

with open('/work/pahwagne/RDsTools/corrections/bs_kinematics/fonll_over_pythia.json') as f:
    bs_kin = json.load(f)
with open('/work/pahwagne/RDsTools/corrections/bs_kinematics/fonll_over_pythia_up.json') as f:
    bs_kin_up = json.load(f)
with open('/work/pahwagne/RDsTools/corrections/bs_kinematics/fonll_over_pythia_down.json') as f:
    bs_kin_down = json.load(f)

with open('/work/pahwagne/RDsTools/corrections/bs_kinematics/bins_pt.json') as f:
    pt_edges = json.load(f)
with open('/work/pahwagne/RDsTools/corrections/bs_kinematics/bins_y.json') as f:
    y_edges = json.load(f)


############################################
pt_edges_str = ", ".join(map(str, pt_edges     )) # of type: '0.01, 0.34, 0.89'
y_edges_str  = ", ".join(map(str, y_edges    )) # of type: '0.01, 0.34, 0.89'

bs_kin_str = "{"

for i,row in enumerate(bs_kin):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  bs_kin_str  += sub

  if i != (len(bs_kin)-1): bs_kin_str += ","

bs_kin_str += "}"

############################################

bs_kin_str_up = "{"

for i,row in enumerate(bs_kin_up):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  bs_kin_str_up  += sub

  if i != (len(bs_kin)-1): bs_kin_str_up += ","

bs_kin_str_up += "}"

############################################

bs_kin_str_down = "{"

for i,row in enumerate(bs_kin_down):
  sub = "{"

  for j,val in enumerate(row):

    sub += str(val)
    if j != len(row)-1: sub += ","

  sub += "}"
  bs_kin_str_down  += sub

  if i != (len(bs_kin)-1): bs_kin_str_down += ","

bs_kin_str_down += "}"

############################################

ROOT.gInterpreter.Declare(r'''double get_rapidity(float eta, float phi, float pt, float m){

  //change coords to x y z
  float p_x = pt * cos (phi);
  float p_y = pt * sin (phi);
  float p_z = pt * sinh(eta);

  float e   = sqrt(pow(p_x,2) + pow(p_y,2) + pow(p_z,2) + pow(m,2));

  float y   = 0.5 * log((e + p_z) / (e - p_z));

  return y;
  }'''
  )

ROOT.gInterpreter.Declare(rf'''double get_w_bs_kin(float b_pt, float b_abs_y){{

    float scale_factors_cpp[{lenx}][{leny}] = {bs_kin_str};
    float pt_edges_cpp [{lenx + 1}] = {{{pt_edges_str}}};
    float y_edges_cpp[{leny + 1}] = {{{y_edges_str}}};

    int index_pt  = 0;
    int index_y = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((b_pt > pt_edges_cpp[i]) && (b_pt < pt_edges_cpp[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    if ((b_abs_y > y_edges_cpp[i]) && (b_abs_y < y_edges_cpp[i+1])){{
        index_y = i;
        break;
      }}
    }}


    return scale_factors_cpp[index_pt][index_y];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_w_bs_kin_up(float b_pt, float b_abs_y){{

    float scale_factors_cpp[{lenx}][{leny}] = {bs_kin_str_up};
    float pt_edges_cpp [{lenx + 1}] = {{{pt_edges_str}}};
    float y_edges_cpp[{leny + 1}] = {{{y_edges_str}}};

    int index_pt  = 0;
    int index_y = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((b_pt > pt_edges_cpp[i]) && (b_pt < pt_edges_cpp[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    if ((b_abs_y > y_edges_cpp[i]) && (b_abs_y < y_edges_cpp[i+1])){{
        index_y = i;
        break;
      }}
    }}


    return scale_factors_cpp[index_pt][index_y];
    }}
    '''
    )

ROOT.gInterpreter.Declare(rf'''double get_w_bs_kin_down(float b_pt, float b_abs_y){{

    float scale_factors_cpp[{lenx}][{leny}] = {bs_kin_str_down};
    float pt_edges_cpp [{lenx + 1}] = {{{pt_edges_str}}};
    float y_edges_cpp[{leny + 1}] = {{{y_edges_str}}};

    int index_pt  = 0;
    int index_y = 0;

    for(int i=0; i<{lenx}; ++i){{

    if ((b_pt > pt_edges_cpp[i]) && (b_pt < pt_edges_cpp[i+1])){{
        index_pt = i;
        break;
      }} 
    }}

    for(int i=0; i<{leny}; ++i){{

    if ((b_abs_y > y_edges_cpp[i]) && (b_abs_y < y_edges_cpp[i+1])){{
        index_y = i;
        break;
      }}
    }}


    return scale_factors_cpp[index_pt][index_y];
    }}
    '''
    )

####################
# dxy err morphing #
####################



ROOT.gInterpreter.Declare(r'''
  #include <correction.h>
  #include <memory>

  auto corrSet = correction::CorrectionSet::from_file("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/quantile_morph_and_trigger_sf_14sep26/dxy_morphing/sigma_dxy_morph_splot_fixed_muonid_fixed_decay_length.json")->at("mu_bs_dxy_e_morph");

  double get_corr_dxy(double pt, double abs_eta, double dxy_err) {
      return corrSet->evaluate({pt, abs_eta, dxy_err});
  }
  
'''
)


print("====> Define branch and create snapshot")
"HOOK_NEW_BRANCH"

print("====> Snapshot DONE")
