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
from bgl_variations_scalar import variations

ham = Hammer()
fbBuffer = IOBuffer

# We don't add the tau decay, becasue that doesn't change the FF of the Bc
ham.include_decay(["BsDsTauNu"])

ff_input_scheme = dict()
ff_input_scheme["BsDs"] = "ISGW2"
ham.set_ff_input_scheme(ff_input_scheme)

ff_schemes  = dict()
ff_schemes['bglvar' ] = {'BsDs':'BGLVar' }
#for i, j in product(range(6), ['up', 'down']):
#        unc = 'e%d%s'%(i,j)
#        ff_schemes['bglvar_%s'%unc] = {'BsDs':'BGLVar_%s'%unc  }
key = "{"
key += "a0: ["  +  str(0.052255946347001495)    + ", " + str(-0.16027634967890908) + ", " + str(0.014141836205563255) + ", " + str(0.0) + "],";
#key += "a0: ["  +  str(0.0)    + ", " + str(0.0) + ", " + str(0.0) + ", " + str(0.0) + "],";
key += "ap: ["  +  str(0.0017893827864468802)   + ", " + str(-0.004691380424494185) + ", " + str( -0.015708534616906505) + ", " + str(0.0) + "],"
#key += "ap: ["  +  str(0.0)   + ", " + str(0.0) + ", " + str( 0.0) + ", " + str(0.0) + "],"
key += "}"
paras = "BstoDsBGLVar: " + key

ham.set_options(paras)

for k, v in ff_schemes.items():
        ham.add_ff_scheme(k, v)
ham.set_units("GeV")
ham.init_run()

#for i, j in product(range(6), ['up', 'down']):
#        unc = 'e%d%s'%(i,j)
#        ham.set_ff_eigenvectors('BstoDs', 'BGLVar_%s'%unc, variations['e%d'%i][j])

fname = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/dstau_BGLVar_13_12_2024_14_24_08/dstau_BGLVar_0.root"
fin = ROOT.TFile.Open(fname)
tree = fin.Get('tree')
maxevents = 3
tree_df = read_root(fname, 'tree', where='(gen_sig == 1)')

pids = []
weights = dict()

for k in ff_schemes.keys():
        weights[k] = []

start = time()
maxevents = maxevents if maxevents>=0 else tree.GetEntries() # total number of events in the files

#sys.exit()
for i, ev in enumerate(tree):

    if not (ev.gen_sig == 1): continue

    if (i+1)>maxevents:
        break
        
    if i%1000==0:
        percentage = float(i)/maxevents*100.
        speed = float(i)/(time()-start)
        eta = datetime.now() + timedelta(seconds=(maxevents-i) / max(0.1, speed))
        print('\t===> processing %d / %d event \t completed %.1f%s \t %.1f ev/s \t ETA %s s' %(i, maxevents, percentage, '%', speed, eta.strftime('%Y-%m-%d %H:%M:%S')))

    ham.init_event()
   
    print(ev.gen_tau_pt, ev.gen_ds_pt, ev.gen_bs_pt, ev.gen_bs_pdgid, ev.gen_sig)

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
    pid = ham.add_process(Bc2JpsiLNu)
    pids.append(pid)

    ham.process_event()
    for k in ff_schemes.keys():
            weights[k].append(ham.get_weight(k))
            print(k, "has weight: ", ham.get_weight(k) )    
    if i>maxevents: break

reduced_tree = tree_df[:len(weights[k])].copy()
for k in ff_schemes.keys():
        reduced_tree['hammer_'+k] = np.nan_to_num(np.array(weights[k])) 
to_root(reduced_tree, 'hammer_output_tau_v2.root', key='tree')
