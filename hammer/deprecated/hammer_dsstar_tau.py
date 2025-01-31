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
from bgl_variations_vector import variations


ham = Hammer()
fbBuffer = IOBuffer

# We don't add the tau decay, becasue that doesn't change the FF of the Bc
#ham.include_decay(["BstoDs*TauNu"])

##################################################################################################################
## BIG REMARK Nr. 1 !!!!!!!!!!!!!!!!!!!!!!!!                                                                    ##
## The string 'BsDs*' you give as key to ff_input_scheme etc shall NOT be the SAME as the one below from the    ##
## HAMMER src code. So it should not be 'BstoDs*' as in:                                                        ## 
## https://gitlab.com/mpapucci/Hammer/-/blob/master/src/FormFactors/BGL/FFBtoDstarBGL.cc?ref_type=heads#L41     ##
## Otherwise you get a seg dault!!                                                                              ##
##################################################################################################################

ff_input_scheme = dict()
ff_input_scheme["BsDs*"] = "ISGW2"
ham.set_ff_input_scheme(ff_input_scheme)

ff_schemes  = dict()
ff_schemes['bglvar' ] = {'BsDs*':'BGLVar' }

key = "{"
#cohen
#key += "avec: ["  +  str(0.004605664681641084)    + ", " + str(-0.002140593040599278) + ", " + str(0.15566982447466055) + "],";
#key += "bvec: ["  +  str(0.003303529928953319)    + ", " + str(-0.004284980385058838) + ", " + str(0.17791644334552834) + "],";
#key += "cvec: ["  +  str(-0.0018867020644757423)  + ", " + str(0.022525216948547932)  + ", " +                               "],";
#key += "dvec: ["  +  str(0.03980443778007538)     + ", " + str(-0.1872442367469107)   + ", " +                               "],";

#harrison
key += "avec: [0.026667, -0.048823, -0.001545],"; 
key += "bvec: [0.413130, -0.075637, -0.250136],";
key += "cvec: [1.206462, 2.327528],";
key += "dvec: [0.209480, -0.861254]";
key += "}"
paras = "BstoDs*BGLVar: " + key
#ham.set_options(paras)

#for i, j in product(range(6), ['up', 'down']):
#        unc = 'e%d%s'%(i,j)
#        ff_schemes['bglvar_%s'%unc] = {'BsDs*':'BGLVar_%s'%unc  }

for k, v in ff_schemes.items():
        ham.add_ff_scheme(k, v)
ham.set_units("GeV")
ham.init_run()

for i, j in product(range(6), ['up', 'down']):
        unc = 'e%d%s'%(i,j)

        ##################################################################################################################
        ## BIG REMARK Nr. 2 !!!!!!!!!!!!!!!!!!!!!!!!                                                                    ##
        ## The string 'BstoDs*' you give to the set_ff_eigenvectors() has to be the the SAME as in the HAMMER src code! ##
        ## See e.g.:                                                                                                    ## 
        ## https://gitlab.com/mpapucci/Hammer/-/blob/master/src/FormFactors/BGL/FFBtoDstarBGL.cc?ref_type=heads#L41     ##
        ## Otherwise all weights are the same!                                                                          ##
        ##################################################################################################################

        ham.set_ff_eigenvectors('BstoDs*', 'BGLVar_%s'%unc, variations['e%d'%i][j])

fname = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/22_07_2024_11_23_05/all_signals_flatChunk_0.root"
fname = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/29_08_2024_09_46_52/all_signals_flatChunk_0.root"
fname = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/16_11_2024_09_45_34/dsstartau_flatChunk_0.root"

fin = ROOT.TFile.Open(fname)
tree = fin.Get('tree')
maxevents = 10
tree_df = read_root(fname, 'tree', where='(gen_sig == 11)')

pids = []
weights = dict()

for k in ff_schemes.keys():
        weights[k] = []

start = time()
maxevents = maxevents if maxevents>=0 else tree.GetEntries() # total number of events in the files

#sys.exit()
for i, ev in enumerate(tree):

    if not (ev.gen_sig == 11): continue

    if (i+1)>maxevents:
        break
        
    if i%1000==0:
        percentage = float(i)/maxevents*100.
        speed = float(i)/(time()-start)
        eta = datetime.now() + timedelta(seconds=(maxevents-i) / max(0.1, speed))
        print('\t===> processing %d / %d event \t completed %.1f%s \t %.1f ev/s \t ETA %s s' %(i, maxevents, percentage, '%', speed, eta.strftime('%Y-%m-%d %H:%M:%S')))

    ham.init_event()
   
    #print(ev.gen_tau_pt, ev.gen_ds_pt, ev.gen_bs_pt, ev.gen_bs_pdgid, ev.gen_sig)

    # Bs -> Ds Tau Nu
    thebs_p4    = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_bs_pt,   ev.gen_bs_eta,  ev.gen_bs_phi,  5.366)
    thetau_p4   = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_tau_pt,  ev.gen_tau_eta, ev.gen_tau_phi, 1.776)#ev.gen_tau_m)
    theds_p4    = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(ev.gen_dsStar_pt,   ev.gen_dsStar_eta,  ev.gen_dsStar_phi,  2.112)#ev.gen_ds_m)
    thenu_p4    = thebs_p4 - thetau_p4 - theds_p4

    thebs       = Particle(FourMomentum(thebs_p4.e(),  thebs_p4.px(),  thebs_p4.py(),  thebs_p4.pz()),  ev.gen_bs_pdgid)
    thetau      = Particle(FourMomentum(thetau_p4.e(), thetau_p4.px(), thetau_p4.py(), thetau_p4.pz()), int(15*ev.gen_ds_charge))
    theds       = Particle(FourMomentum(theds_p4.e(),  theds_p4.px(),  theds_p4.py(),  theds_p4.pz()),  ev.gen_dsStar_pdgid)
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
#to_root(reduced_tree, 'hammer_output_tau_star_v2.root', key='tree')
