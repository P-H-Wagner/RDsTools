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

scale_factors        = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_fullBPark.txt")
scale_factors_err_hi = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_hi_fullBPark.txt")
scale_factors_err_lo = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/trigger_sf_err_lo_fullBPark.txt")
pt_edges             = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/pt_edges_fullBPark.txt")
dxy_edges            = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/ricc/dxy_edges_fullBPark.txt")

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

############################################
#for i,row in enumerate(scale_factors_err):
#  sub = "{"
#
#  for j,val in enumerate(row):
#
#    sub += str(val)
#    if j != len(row)-1: sub += ","
#
#  sub += "}"
#  scale_factors_err_str += sub
#
#  if i != (len(scale_factors_err)-1): scale_factors_err_str += ","
#
#scale_factors_err_str += "}"
#


print(scale_factors_str)
print(pt_edges_str)
print(dxy_edges_str)

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

print("====> Define branch and create snapshot")
"HOOK_NEW_BRANCH"

print("====> Snapshot DONE")
