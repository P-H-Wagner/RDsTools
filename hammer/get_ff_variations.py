"""
All numbers from https://doi.org/10.1103/PhysRevD.100.094503
Use supplementary material: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.094503#supplemental
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import multivariate_normal
from itertools import combinations
import yaml
import json

np.set_printoptions(precision=5, suppress=True, linewidth=200)

def get_ff_variations( names_hammer, coeff, coeff_uncertainties, corr, model):

  print(names_hammer)
  print(coeff)
  print(coeff_uncertainties)

  isBGL = False

  if "bgl" in model and "scalar" in model:
    isBGL = True
  
    #remove zero columns/rows with index 3 and 7
    corr_new                = np.delete(np.delete(corr, [3, 7], axis=0), [3, 7], axis=1)
    coeff_new               = np.delete(coeff, [3,7])
    coeff_uncertainties_new = np.delete(coeff_uncertainties, [3,7])
    names_hammer_new        = np.delete(names_hammer, [3,7])

    corr = corr_new
    coeff = coeff_new
    coeff_uncertainties = coeff_uncertainties_new
    names_hammer = list(names_hammer_new)
  
  #name the eigenvectors like
  names_eigenvectors = names_hammer 

  # here we need elementwise multiplication, not matrix multiplication
  # need to use atleast_2d to really traspose the vector
  
  # Covariance = Correlation * sigma_x * sigma_y 
  cov = np.atleast_2d(coeff_uncertainties).T * corr * coeff_uncertainties
  print("CORRELATION MATRIX")
  print(corr)
  print("UNCERTAINTIES")
  print(coeff_uncertainties)
  print("COVARIANCE MATRIX")
  print(cov)
  
  # Get Eigenvectors and eigenvalues of the covariance matrix
  eVals, eVecs = np.linalg.eig(cov)    

  for i,eig in enumerate(eVals):
    if eig < 0: 
      #print(f"!! Eigenvalue is negative: {eig} ====> setting to 0.0")
      eVals[i] = 0.0      

  #print(eVals)
  # eigenvalues matrix
  diag = np.identity(len(coeff)) * eVals
  print("EIGENVALUES")
  print(eVals)
  print("EIGENVECTOR")
  print(eVecs)

  #for i in range(8):
  #  print("option P: eigenvector ", i, " is ", eVecs[:,i])
 
  # rotate back to original basis
  cov_new = eVecs.dot(diag.dot(eVecs.T))
  rtol = 1e-5
  #print("closure test: rebuild the covariance matrix from the eigenvectors, eigenvalues")
  #print("\tpassed at {} level?".format(rtol), np.isclose(cov, cov_new, rtol=1e-5).all())
  
  # sigmas should be squared root of the eigenvalues
  eUncs = np.nan_to_num(np.sqrt(eVals))
  
  # principal components
  #print("principal components")
  principal_comp = np.atleast_2d(eUncs*eVecs).T
  #print(principal_comp)
 
  print("\n\n")
  
  variations = dict()
  
  # create a structure, i.e. a dictionary that allows:
  # for each movement along the direction of the j-th eigenvector
  # define two possible ways, up and down

  names_to_fill = names_hammer.copy()
  names_to_fill = list(names_to_fill)

  if isBGL:
    names_to_fill.append("a03"); names_to_fill.append("ap3")

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
  print(principal_comp)
 
  # now fill this dictionary
  for i,name1 in enumerate(names_hammer):
    for j,name2 in enumerate(names_hammer):

      name1 = str(name1) #keep this for saving into yaml file!
      name2 = str(name2) #keep this for saving into yaml file!
     
      print("i is ", i, "j is ", j)
      print("setting: " + name1 + "up" + "delta_" + name2 + " to ", float(principal_comp[j, i]))
      variations[name2]["up"  ]["delta_" + name1] = float(principal_comp[j, i])
      variations[name2]["down"]["delta_" + name1] = float(-principal_comp[j, i])

  print("AFTER") 
  print(json.dumps(variations, indent=4, sort_keys=True))

  # save values of up/down variations
  with open(f"{model}_variations.py", "w") as fout:
      print("variations =", variations, file=fout)
      #yaml.dump(variations, fout, default_flow_style=False)
  
  # Plot 2D correlations
  
  for i,j in combinations(range(11),2): 
      #print("studying %s vs %s" %(names_hammer[i], names_hammer[j]))
      
      cval = np.array([
          coeff[i],
          coeff[j],
      ])
      
      uncs = np.array([
          coeff_uncertainties[i],
          coeff_uncertainties[j],
      ])
      
      # extract submatrix from cov matrix
      minicov = cov[np.ix_([i,j],[i,j])]
      
      # extract the relevant principal components
      mini_principal_comp = principal_comp.T[np.ix_([i,j],[i,j])]
      
  #     xmin = -2.*uncs[0] + cval[0]
  #     xmax =  2.*uncs[0] + cval[0]
  #     ymin = -2.*uncs[1] + cval[1]
  #     ymax =  2.*uncs[1] + cval[1]
  
      # need to keep the same scaling for x- and y-axis
      # because arrows, quiver... don"t scale.
  #     max_unc = max(uncs)
  #     xmin = -2.*max_unc + cval[0]
  #     xmax =  2.*max_unc + cval[0]
  #     ymin = -2.*max_unc + cval[1]
  #     ymax =  2.*max_unc + cval[1]
  
      min_unc = min(uncs)
      xmin = -6.*min_unc + cval[0]
      xmax =  6.*min_unc + cval[0]
      ymin = -6.*min_unc + cval[1]
      ymax =  6.*min_unc + cval[1]
  
      x = np.linspace(xmin, xmax, 500)
      y = np.linspace(ymin, ymax, 500)
      X, Y = np.meshgrid(x,y)
      pos = np.empty(X.shape + (2,))
      pos[:, :, 0] = X; 
      pos[:, :, 1] = Y
      rv = multivariate_normal(mean=cval, cov=minicov)
      
      # 2D plot
      plt.clf()
      fig = plt.figure(figsize=(5,5))
      aspect = abs(xmax-xmin)/abs(ymax-ymin)
      subplt = fig.add_subplot(111, box_aspect=1., aspect=aspect)
  
      Z = rv.pdf(pos)
  
      plt.imshow(Z, origin="lower", extent=[xmin, xmax, ymin, ymax])
  
      levels = [
          np.power(np.e, -9. ) * abs(np.max(Z)-np.min(Z))+np.min(Z),  # 3 sigma 
          np.power(np.e, -2. ) * abs(np.max(Z)-np.min(Z))+np.min(Z),  # 2 sigma 
          np.power(np.e, -0.5) * abs(np.max(Z)-np.min(Z))+np.min(Z),  # 1 sigma    
          abs(np.max(Z))                                           ,  # max
      ] 
      levels_str = [r"3 $\sigma$", r"2 $\sigma$", r"1 $\sigma$", r"0 $\sigma$"]
      contours = subplt.contour(X, Y, Z, levels=levels, colors="silver")
      fmt = {}
      for l, s in zip(contours.levels, levels_str):
          fmt[l] = s
      subplt.clabel(contours, contours.levels[:-1], inline=True, fmt=fmt)
      
      origin = np.array([np.ones(2)*cval[0], np.ones(2)*cval[1]])
      subplt.quiver(*origin, mini_principal_comp[:,0], mini_principal_comp[:,1], units="xy", color=["r","b"], angles="xy", scale_units="xy", scale=1.)
  
      plt.xlabel(names_hammer[i])
      plt.ylabel(names_hammer[j])
  
  #     for ix, iy in mini_principal_comp+origin:
  #         plt.text(ix, iy, "(%.2f, %.2f)"%(ix, iy))
  
      subplt.scatter(*cval)
      plt.text(*cval, "({:.1e}, {:.1e})".format(*cval), c="silver")
      
      plot_margin = 0.2
      plt.subplots_adjust(left=0+plot_margin, bottom=0+plot_margin)#, right=1, top=1, wspace=0, hspace=0)
      plt.savefig(f"{model}_contour_%s_vs_%s.pdf" %(names_hammer[i], names_hammer[j]))
      
      break
  
  # closure test
  # rebuild the covariance matrix from the principal components
  cov_from_pc = principal_comp.T.dot(np.identity(len(principal_comp)).dot(principal_comp))
  rtol = 1e-5
  #print("closure test: rebuild the covariance matrix from the principal components")
  #print("\tpassed at {} level?".format(rtol), np.isclose(cov, cov_from_pc, rtol=1e-5).all())
  
  # principal components + multivariate gaus mean
  #print("principal components + multivariate gaus mean")
  #print(principal_comp+coeff)
  
  return  
