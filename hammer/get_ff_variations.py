'''
All numbers from https://doi.org/10.1103/PhysRevD.100.094503
Use supplementary material: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.094503#supplemental
'''

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import multivariate_normal
from itertools import combinations

def get_ff_variations( names_hammer, coeff, coeff_uncertainties, corr, model):

  #name the eigenvectors like: "e0, e1, .., e9"
  names_eigenvectors = ['e%d' %i for i in range(len(names_hammer))]

  # here we need elementwise multiplication, not matrix multiplication
  # need to use atleast_2d to really traspose the vector
  
  # Covariance = Correlatin * sigma_x * sigma_y 
  cov = np.atleast_2d(coeff_uncertainties).T * corr * coeff_uncertainties
  print('covariance matrix')
  print(cov)
  
  # Get Eigenvectors and eigenvalues of the covariance matrix
  eVals, eVecs = np.linalg.eig(cov)    
  
  print('eigenvectors V (components wrt x,y)')
  for i, iev in enumerate(eVecs):
      print(i, iev)
  
  # eigenvalues matrix
  diag = np.identity(len(coeff)) * eVals
  #print('eigenvalues matrix W')
  #print(diag)
  
  # rotate back to original basis
  cov_new = eVecs.dot(diag.dot(eVecs.T))
  rtol = 1e-5
  print('closure test: rebuild the covariance matrix from the eigenvectors, eigenvalues')
  print('\tpassed at {} level?'.format(rtol), np.isclose(cov, cov_new, rtol=1e-5).all())
  
  # sigmas should be squared root of the eigenvalues
  eUncs = np.nan_to_num(np.sqrt(eVals))
  
  # principal components
  #print('principal components')
  principal_comp = np.atleast_2d(eUncs*eVecs).T
  #print(principal_comp)
  
  print('\n\n')
  
  variations = dict()
  
  # create a structure, i.e. a dictionary that allows:
  # for each movement along the direction of the j-th eigenvector
  # define two possible ways, up and down
  for iname in names_eigenvectors:
      variations[iname] = dict()
      variations[iname]['up'  ] = dict()
      variations[iname]['down'] = dict()
      # for each of these, specify ho much distance needs to be travelled in
      # the x, y basis
      for jname in names_hammer:
          variations[iname]['up'  ]['delta_'+jname] = 0.
          variations[iname]['down']['delta_'+jname] = 0.
  
  # now fill this dictionary
  for i in range(eUncs.shape[0]):
      for j in range(eVecs.shape[0]):
          variations[names_eigenvectors[j]]['up'  ]['delta_'+names_hammer[i]] =  principal_comp[j, i]
          variations[names_eigenvectors[j]]['down']['delta_'+names_hammer[i]] = -principal_comp[j, i]
  
  # save values of up/down variations
  with open(f'{model}_variations.py', 'w') as fout:
      print('variations =', variations, file=fout)
  
  #print ('='*80+'\n\n')
  #print(variations)
  #print ('='*80+'\n\n')
  
  
  # Plot 2D correlations
  
  for i,j in combinations(range(11),2): 
      print('studying %s vs %s' %(names_hammer[i], names_hammer[j]))
      
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
      # because arrows, quiver... don't scale.
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
  
      plt.imshow(Z, origin='lower', extent=[xmin, xmax, ymin, ymax])
  
      levels = [
          np.power(np.e, -9. ) * abs(np.max(Z)-np.min(Z))+np.min(Z),  # 3 sigma 
          np.power(np.e, -2. ) * abs(np.max(Z)-np.min(Z))+np.min(Z),  # 2 sigma 
          np.power(np.e, -0.5) * abs(np.max(Z)-np.min(Z))+np.min(Z),  # 1 sigma    
          abs(np.max(Z))                                           ,  # max
      ] 
      levels_str = [r'3 $\sigma$', r'2 $\sigma$', r'1 $\sigma$', r'0 $\sigma$']
      contours = subplt.contour(X, Y, Z, levels=levels, colors='silver')
      fmt = {}
      for l, s in zip(contours.levels, levels_str):
          fmt[l] = s
      subplt.clabel(contours, contours.levels[:-1], inline=True, fmt=fmt)
      
      origin = np.array([np.ones(2)*cval[0], np.ones(2)*cval[1]])
      subplt.quiver(*origin, mini_principal_comp[:,0], mini_principal_comp[:,1], units='xy', color=['r','b'], angles='xy', scale_units='xy', scale=1.)
  
      plt.xlabel(names_hammer[i])
      plt.ylabel(names_hammer[j])
  
  #     for ix, iy in mini_principal_comp+origin:
  #         plt.text(ix, iy, '(%.2f, %.2f)'%(ix, iy))
  
      subplt.scatter(*cval)
      plt.text(*cval, '({:.1e}, {:.1e})'.format(*cval), c='silver')
      
      plot_margin = 0.2
      plt.subplots_adjust(left=0+plot_margin, bottom=0+plot_margin)#, right=1, top=1, wspace=0, hspace=0)
      plt.savefig(f'{model}_contour_%s_vs_%s.pdf' %(names_hammer[i], names_hammer[j]))
      
      break
  
  # closure test
  # rebuild the covariance matrix from the principal components
  cov_from_pc = principal_comp.T.dot(np.identity(len(principal_comp)).dot(principal_comp))
  rtol = 1e-5
  print('closure test: rebuild the covariance matrix from the principal components')
  print('\tpassed at {} level?'.format(rtol), np.isclose(cov, cov_from_pc, rtol=1e-5).all())
  
  # principal components + multivariate gaus mean
  #print('principal components + multivariate gaus mean')
  #print(principal_comp+coeff)
  
  return  
