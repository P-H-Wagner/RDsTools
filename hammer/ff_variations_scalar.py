'''
All numbers from https://doi.org/10.1103/PhysRevD.100.094503
Use supplementary material: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.094503#supplemental
'''

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import multivariate_normal
from itertools import combinations
from get_ff_variations import get_ff_variations

# don't use numpy matrices
# https://numpy.org/doc/stable/reference/generated/numpy.matrix.html

names_bgl_coefficients = [
    'f0_a0' ,
    'f0_a1' ,
    'f0_a2' ,
    'f0_a3' ,
    'fp_a0' ,
    'fp_a1' ,
    'fp_a2' , 
    'fp_a3' 
]

names_bgl_coefficients_hammer = [
    'ap0',
    'ap1',
    'ap2',
    'ap3',
    'a00',
    'a01',
    'a02',
    'a03'
]

#name the eigenvectors like: "e0, e1, .., e9"
names_eigenvectors = ['e%d' %i for i in range(len(names_bgl_coefficients_hammer))]

#fill the coefficients according to the paper above 
# the referenced paper above only considers N = 2 extension, thus coeff a3,b3 and their uncertanties are missing!
bgl_coefficients = np.array([
 0.052255946347001495  , # f0  a0 --> Hammer ap0  
-0.16027634967890908   , # f0  a1 --> Hammer ap1 
 0.014141836205563255  , # f0  a2 --> Hammer ap2 
 0.0                   , # f0  a3 --> Hammer ap3 
 0.0017893827864468802 , # f+  b0 --> Hammer a00 
-0.004691380424494185  , # f+  b1 --> Hammer a01 
-0.015708534616906505  , # f+  b2 --> Hammer a02 
 0.0                   , # f+  b3 --> Hammer a03 
]).astype(np.float64)

bgl_coefficient_uncertainties = np.array([
0.0007906527485083701  , # f0  a0 --> Hammer ap0
0.016991421017241458   , # f0  a1 --> Hammer ap1
0.03169323419790213    , # f0  a2 --> Hammer ap2
0.0                    , # f0  a3 --> Hammer ap3
0.00004543889976730505 , # f+  b0 --> Hammer a00
0.0008662317680431299  , # f+  b1 --> Hammer a01
0.0047035296806361825  , # f+  b2 --> Hammer a02
0.0                    , # f+  b3 --> Hammer a03
]).astype(np.float64)


# correlation matrix (normalised covariance matrix)

import numpy as np
correlation_complete = np.array([
    [ 1.00000, -0.26188, -0.07456, 0.0,   0.43994, -0.00314,  0.13858, 0.0],
    [-0.26188,  1.00000, -0.18571, 0.0,   0.28968,  0.10064, -0.09532, 0.0],
    [-0.07456, -0.18571,  1.00000, 0.0,  -0.16234,  0.12522, -0.21066, 0.0],
    [     0.0,      0.0,      0.0, 0.0,       0.0,      0.0,      0.0, 0.0],    
    [ 0.43994,  0.28968, -0.16234, 0.0,   1.00000, -0.72825, -0.23200, 0.0],
    [-0.00314,  0.10064,  0.12522, 0.0,  -0.72825,  1.00000,  0.12053, 0.0],
    [ 0.13858, -0.09532, -0.21066, 0.0,  -0.23200,  0.12053,  1.00000, 0.0],
    [     0.0,      0.0,      0.0, 0.0,       0.0,      0.0,      0.0, 0.0]    
]).astype(np.float64)

corr  = correlation_complete

get_ff_variations(names_bgl_coefficients_hammer, bgl_coefficients, bgl_coefficient_uncertainties, corr, "bgl_scalar")

