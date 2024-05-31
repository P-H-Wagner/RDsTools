import ROOT

bsMass_ = 5.36688

baseline = ' & '.join([
f'(dsMu_m < {bsMass_})',
'(k1_charge*k2_charge <0)',
'(mu_charge*pi_charge < 0)',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
'tv_prob > 0.1',
#'((cosPiK1 < -0.3) || (cosPiK1 > 0.3))',
'(fv_prob > 0.1)'
])

bkg = ' & '.join([
#'k1_charge*k2_charge > 0',
#'mu_charge*pi_charge > 0',
f'dsMu_m > {bsMass_}',
#'mu_rel_iso_03 > 0.3'
])

selections = {"baseline": baseline, "bkg": bkg}

# pick the skim selection
name = "HOOK_SELECTION"
selec = selections[name]
# skim it!

#create rdf
files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/HOOK_DATE_TIME/*"
chain = ROOT.TChain("tree")
chain.Add(files)

#output
destination = f"/scratch/pahwagne/skimmed_HOOK_SELECTION_HOOK_DATE_TIME.root"
print("saving to:", destination)
print("====> Create rdf")
df = ROOT.RDataFrame(chain)
print("====> rdf DONE")

print("====> Define branch and create snapshot")
"HOOK_NEW_BRANCH"

print("====> Snapshot DONE")
