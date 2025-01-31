import numpy as np
np.set_printoptions(precision=10, suppress=True, linewidth=200)
import yaml



from bgl_vector_cov_mat.py import cov
print(cov)

sig = []

for i in range(len(cov[0,:])):
  sig = np.sqrt(cov[i][i])
  
# Get Eigenvectors and eigenvalues of the covariance matrix
eVals, eVecs = np.linalg.eig(cov)    

for i,eig in enumerate(eVals):
  if eig < 0: 
    print(f"!! Eigenvalue is negative: {eig} ====> setting to 0.0")
    eVals[i] = 0.0      

# eigenvalues matrix
diag = np.identity(len(cov)) * eVals
print("EIGENVALUES")
print(eVals)
print("EIGENVECTOR")
print(eVecs)

#for i in range(8):
#  print("option P: eigenvector ", i, " is ", eVecs[:,i])

# rotate back to original basis
cov_new = eVecs.dot(diag.dot(eVecs.T))
rtol = 1e-5
print("closure test: rebuild the covariance matrix from the eigenvectors, eigenvalues")
print("\tpassed at {} level?".format(rtol), np.isclose(cov, cov_new, rtol=1e-5).all())

# sigmas should be squared root of the eigenvalues
eUncs = np.nan_to_num(np.sqrt(eVals))

# principal components

pc = np.zeros_like(cov_new)

for i,l in enumerate(eUncs):
  for j in range(len(eUncs)):
    pc[j,i] = eUncs[i]*eVecs[j,i]

  print("principal component ", i)
  print(l*eVecs[:,i])

print("\n\n")

variations = dict()
variations.clear()  

# create a structure, i.e. a dictionary that allows:
# for each movement along the direction of the j-th eigenvector
# define two possible ways, up and down

names_to_fill = ["e1","e2","e3","e4","e5","e6", "e7", "e8", "e9", "e10"] 
names_hammer = names_to_fill
names_to_fill = list(names_to_fill)


for name1 in names_to_fill:

    name1 = str(name1) #keep this for saving into yaml file!
    variations[name1] = dict()
    variations[name1]["up"  ] = dict()
    variations[name1]["down"] = dict()

    # for each of these, specify ho much distance needs to be travelled in
    # the x, y basis

    for name2 in names_to_fill:

        name2 = str(name2) #keep this for saving into yaml file!
        variations[name1]["up"  ]["delta_"+name2] = float(0.)
        variations[name1]["down"]["delta_"+name2] = float(0.)

print("PRINCIPAL COMP")
print(pc)

# now fill this dictionary
for i,name1 in enumerate(names_hammer):
  for j,name2 in enumerate(names_hammer):
    name1 = str(name1) #keep this for saving into yaml file!
    name2 = str(name2) #keep this for saving into yaml file!
   
    variations[name1]["up"  ]["delta_" + name2] = float(pc[j, i])
    variations[name1]["down"]["delta_" + name2] = float(-pc[j, i])

#print("AFTER") 
#print(json.dumps(variations, indent=4, sort_keys=True))
print("PRINCIPAL COMPONENT MATRIX") 
print(pc)

with open(f"bgl_vector_variations.yaml", "w") as fout:
  yaml.dump(variations, fout, default_flow_style=False)

