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

pu_weights    = np.loadtxt("/work/pahwagne/RDsTools/pileup/pu_weights.txt")
scale_factors = np.loadtxt("/work/pahwagne/RDsTools/data/trigger_scale_factors/scale_factor_matrix.txt")
pt_edges  = np.loadtxt("/work/pahwagne/RDsTools/data/trigger_scale_factors/pt_edges.txt")
dxy_edges = np.loadtxt("/work/pahwagne/RDsTools/data/trigger_scale_factors/dxy_edges.txt")

#convert into string such that it can be handed to the cpp process line
pu_weights_str    = ", ".join(map(str, pu_weights   )) # of type: '0.01, 0.34, 0.89'
pt_edges_str  = ", ".join(map(str, pt_edges     )) # of type: '0.01, 0.34, 0.89'
dxy_edges_str = ", ".join(map(str, dxy_edges    )) # of type: '0.01, 0.34, 0.89'

#scale_factors is a matrix!
lenx = scale_factors.shape[0]
leny = scale_factors.shape[1]

scale_factors_str = "{"

for i,row in enumerate(scale_factors):
  sub = "{"

  for val in row:

    sub += str(val)
    if val != row[-1]: sub += ","

  sub += "}"
  scale_factors_str  += sub

  if i != (len(scale_factors)-1): scale_factors_str += ","

scale_factors_str += "}"
print(scale_factors_str)
print(pt_edges_str)
print(dxy_edges_str)

ROOT.gInterpreter.Declare(rf'''double get_pu_weights(int npv){{

    float weights_cpp[] = {{{pu_weights_str}}};
    return weights_cpp[npv];
    }}
    '''
    )



ROOT.gInterpreter.Declare(rf'''double get_scale_factors(float mu_pt, float dxy_mu_sig_pv){{

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



print("====> Define branch and create snapshot")
"HOOK_NEW_BRANCH"

print("====> Snapshot DONE")
