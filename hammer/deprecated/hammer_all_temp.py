'''
Script that computes the form factor weights for the Bs -> Ds mu nu process
Starting from Kiselev FF (standard MC FF), this script reweights to BGL FF (updated)
It takes as input a flat ntupla (like output of inspector)
'''
import sys
sys.path.append("/work/pahwagne/RDsTools/hammer/Hammer-install/lib/python3.8/site-packages")

import ROOT
from itertools import product
from time import time
from datetime import datetime, timedelta
from hammer.hammerlib import Hammer, IOBuffer, Particle, Process, FourMomentum, WTerm
from hammer import hepmc, pdg
from root_pandas import read_root, to_root
import numpy as np
import os
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/hammer"))
from bgl_scalar_variations import variations as variations_scalar
from bgl_vector_variations import variations as variations_vector
from bgl_scalar_central_values import central_values as central_values_scalar
from bgl_vector_central_values import central_values as central_values_vector

model_name = "bglvar"
model_name_hammer = "BGLVar"


def prepare_hammer(decay_chain, name, name_hammer, model_name, model_name_hammer, n_coeff, central_values):
  
  ham = Hammer()
  fbBuffer = IOBuffer
  
  ham.include_decay([decay_chain])
  
  ff_input_scheme = dict()
  ff_input_scheme[name] = input_schemes[decay_chain]
  
  ham.set_ff_input_scheme(ff_input_scheme)
  
  ff_schemes  = dict()
  ff_schemes[model_name ] = {name:model_name_hammer}
 
  ham.set_options(name_hammer + model_name_hammer  + ": " + central_values)
 
  for i, j in product(range(n_coeff), ['up', 'down']):
    unc = f"e{i}{j}"
    ff_schemes[ f"{model_name}_{unc}"] = {name : f"{model_name_hammer}_{unc}"}
  
  for k, v in ff_schemes.items():
          ham.add_ff_scheme(k, v)
  
  ham.set_units("GeV")
  ham.init_run()
  
  for i, j in product(range(n_coeff), ['up', 'down']):
          unc = f"e{i}{j}"
          ham.set_ff_eigenvectors(name_hammer, f"{model_name_hammer}_{unc}", variations[name][f"e{i}"][j])

  return ham, ff_schemes
  
def get_tree_with_weights( hammer_instances, ff_schemes, input_tree):

  default = "BsDs*TauNu"

  fin = ROOT.TFile.Open(input_tree)
  tree = fin.Get('tree')
  maxevents = HOOK_MAX_EVENTS
  tree_df = read_root(input_tree, 'tree') #, where='(gen_sig == 1)')
  
  pids = []
  weights = dict()
  
  for k in ff_schemes[default].keys():
          weights[k] = []
  
  start = time()
  maxevents = maxevents if maxevents>=0 else tree.GetEntries() # total number of events in the files
  
  #sys.exit()
  for i, ev in enumerate(tree):
 
      #print("Event:", i) 

      if (i+1)>maxevents:
          break
          
      if i%1000==0:
          percentage = float(i)/maxevents*100.
          speed = float(i)/(time()-start)
          eta = datetime.now() + timedelta(seconds=(maxevents-i) / max(0.1, speed))
          print('\t===> processing %d / %d event \t completed %.1f%s \t %.1f ev/s \t ETA %s s' %(i, maxevents, percentage, '%', speed, eta.strftime('%Y-%m-%d %H:%M:%S')))

      if (ev.gen_sig == 1): 

        #print(" =====> This is a Ds Tau event!")

        decay = "BsDsTauNu"   
 
        #ham.init_event()
        hammer_instances[decay].init_event()
   
        #print(ev.gen_tau_pt, ev.gen_ds_pt, ev.gen_bs_pt, ev.gen_bs_pdgid, ev.gen_sig)
    
        # Bs -> Ds Tau Nu
        thebs_p4    = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_bs_pt,   ev.gen_bs_eta,  ev.gen_bs_phi,  5.366)
        thetau_p4   = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_tau_pt,  ev.gen_tau_eta, ev.gen_tau_phi, 1.776)#ev.gen_tau_m)
        theds_p4    = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_ds_pt,   ev.gen_ds_eta,  ev.gen_ds_phi,  1.968)#ev.gen_ds_m)
        thenu_p4    = thebs_p4 - thetau_p4 - theds_p4
    
        thebs       = Particle(FourMomentum(thebs_p4.e(),  thebs_p4.px(),  thebs_p4.py(),  thebs_p4.pz()),  ev.gen_bs_pdgid)
        thetau      = Particle(FourMomentum(thetau_p4.e(), thetau_p4.px(), thetau_p4.py(), thetau_p4.pz()), int(15*ev.gen_ds_charge))
        theds       = Particle(FourMomentum(theds_p4.e(),  theds_p4.px(),  theds_p4.py(),  theds_p4.pz()),  ev.gen_ds_pdgid)
        thenu       = Particle(FourMomentum(thenu_p4.e(),  thenu_p4.px(),  thenu_p4.py(),  thenu_p4.pz()),  int(16 * -1 * ev.gen_ds_charge))
     
    
        Bc2JpsiLNu = Process()
        
        # each of these add_particle operations returns an index, needed to define vertices 
        thebs_idx   = Bc2JpsiLNu.add_particle(thebs  )
        theds_idx   = Bc2JpsiLNu.add_particle(theds  )
        thetau_idx  = Bc2JpsiLNu.add_particle(thetau  )
        thenu_idx   = Bc2JpsiLNu.add_particle(thenu  )
    
        # define decay vertex
        Bc2JpsiLNu.add_vertex(thebs_idx, [theds_idx, thetau_idx, thenu_idx])
    
        # save process id to later retrieve the per-event weight
        #pid = ham.add_process(Bc2JpsiLNu)
        pid = hammer_instances[decay].add_process(Bc2JpsiLNu)
  
  
        pids.append(pid)
    
        hammer_instances[decay].process_event()
        for i,k in enumerate(ff_schemes[default].keys()):
                #weights[k].append(ham.get_weight(k))
                
                #print(k)
                if i < 17: 
                  #print("filling weight: ", hammer_instances[decay].get_weight(k) ); 
                  weights[k].append(hammer_instances[decay].get_weight(k))
                else: weights[k].append(-1.0)

      if (ev.gen_sig == 11): 

        #print(" =====> This is a Ds* Tau event!")

        decay = "BsDs*TauNu"   
 
        #ham.init_event()
        hammer_instances[decay].init_event()
   
        # Bs -> Ds Tau Nu
        thebs_p4    = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_bs_pt,       ev.gen_bs_eta,      ev.gen_bs_phi,      5.366)
        thetau_p4   = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_tau_pt,      ev.gen_tau_eta,     ev.gen_tau_phi,     1.776)#ev.gen_tau_m)
        thedsStar_p4= ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_dsStar_pt,   ev.gen_dsStar_eta,  ev.gen_dsStar_phi,  2.112)#ev.gen_ds_m)
        thenu_p4    = thebs_p4 - thetau_p4 - thedsStar_p4
    
        thebs       = Particle(FourMomentum(thebs_p4.e(),      thebs_p4.px(),      thebs_p4.py(),      thebs_p4.pz()),      ev.gen_bs_pdgid)
        thetau      = Particle(FourMomentum(thetau_p4.e(),     thetau_p4.px(),     thetau_p4.py(),     thetau_p4.pz()),     int(15*ev.gen_ds_charge))
        thedsStar   = Particle(FourMomentum(thedsStar_p4.e(),  thedsStar_p4.px(),  thedsStar_p4.py(),  thedsStar_p4.pz()),  433 * ev.gen_ds_charge)
        thenu       = Particle(FourMomentum(thenu_p4.e(),      thenu_p4.px(),      thenu_p4.py(),      thenu_p4.pz()),      int(16 * -1 * ev.gen_ds_charge))
     
    
        Bc2JpsiLNu = Process()
        
        # each of these add_particle operations returns an index, needed to define vertices 
        thebs_idx       = Bc2JpsiLNu.add_particle(thebs  )
        thedsStar_idx   = Bc2JpsiLNu.add_particle(thedsStar  )
        thetau_idx      = Bc2JpsiLNu.add_particle(thetau  )
        thenu_idx       = Bc2JpsiLNu.add_particle(thenu  )
    
        # define decay vertex
        Bc2JpsiLNu.add_vertex(thebs_idx, [thedsStar_idx, thetau_idx, thenu_idx])
    
        # save process id to later retrieve the per-event weight
        #pid = ham.add_process(Bc2JpsiLNu)
        pid = hammer_instances[decay].add_process(Bc2JpsiLNu)
  
  
        pids.append(pid)
    
        hammer_instances[decay].process_event()
        for k in ff_schemes[default].keys():
                #print(k)
                #print("filling weight: ", hammer_instances[decay].get_weight(k) ); 
                #weights[k].append(ham.get_weight(k))
                weights[k].append(hammer_instances[decay].get_weight(k))


      else:
        #for now append 0.0 for all other channels 
        for k in ff_schemes[default].keys():
                weights[k].append( -1.0 )

      if i>maxevents: break

  # the reduced tree has the same length (weights[k] has length maxevents)  
  reduced_tree = tree_df[:len(weights[k])].copy()
  print(reduced_tree)
  print(len(reduced_tree.index))
  print(np.nan_to_num(np.array(weights[k]))) 
  print(len(np.nan_to_num(np.array(weights[k])))) 
  for k in ff_schemes[default].keys():
          reduced_tree['hammer_'+k] = np.nan_to_num(np.array(weights[k])) 
  to_root(reduced_tree, 'HOOK_FILE_OUT', key='tree')

#########################
# CREATE HAMMER OBJECTS #
#########################

hammer_instances     = {}
ff_schemes           = {}
#
input_schemes               = {"BsDsTauNu": "ISGW2",           "BsDs*TauNu": "ISGW2", "BsDsMuNu": "CLN", "BsDs*MuNu": "CLN"}
variations                  = {"BsDs":      variations_scalar, "BsDs*":      variations_vector}

#####################
# Call it on Ds Tau #
#####################
decay       = "BsDsTauNu"
name        = "BsDs"
name_hammer = "BstoDs"

ham , ff_scheme      = prepare_hammer(decay,  name, name_hammer, model_name, model_name_hammer, 8,central_values_scalar)

hammer_instances[decay]     = ham
ff_schemes[decay]           = ff_scheme

##########################
# Call it on Ds Star Tau #
##########################
decay       = "BsDs*TauNu"
name        = "BsDs*"
name_hammer = "BstoDs*"

ham, ff_scheme      = prepare_hammer(decay,  name, name_hammer, model_name, model_name_hammer, 10, central_values_vector)

hammer_instances[decay]     = ham
ff_schemes[decay]           = ff_scheme


#input_tree = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/22_07_2024_11_23_05/all_signals_flatChunk_0.root"
#input_tree = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/29_08_2024_09_46_52/all_signals_flatChunk_0.root"

input_tree  = "HOOK_FILE_IN"
get_tree_with_weights( hammer_instances, ff_schemes, input_tree)


