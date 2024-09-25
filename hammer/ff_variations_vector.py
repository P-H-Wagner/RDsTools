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
    'g_a0' ,
    'g_a1' ,
    'g_a2' ,
    'f_a0' ,
    'f_a1' ,
    'f_a2' ,
#   'F1_a0',
    'F1_a1',
    'F1_a2',
    'F2_a0',
    'F2_a1',
    'F2_a2',
]

names_bgl_coefficients_hammer = [
    'a0',
    'a1',
    'a2',
    'b0',
    'b1',
    'b2',
#   'c0',  <== missing in Hammer
    'c1',
    'c2',
    'd0',
    'd1',
    'd2',
]

#name the eigenvectors like: "e0, e1, .., e9"
names_eigenvectors = ['e%d' %i for i in range(len(names_bgl_coefficients_hammer))]

#fill the coefficients according to the paper above
bgl_coefficients = np.array([
 0.004605664681641084   , # g  a0 --> Hammer a0
-0.002140593040599278   , # g  a1 --> Hammer a1
 0.15566982447466055    , # g  a2 --> Hammer a2
 0.003303529928953319   , # f  a0 --> Hammer b0
-0.004284980385058838   , # f  a1 --> Hammer b1
 0.17791644334552834    , # f  a2 --> Hammer b2
 #0.000500939732384485   , # F1 a0 --> Hammer c0 <=== MISSING FROM HAMMER
-0.0018867020644757423  , # F1 a1 --> Hammer c1
 0.022525216948547932   , # F1 a2 --> Hammer c2
 0.03980443778007538    , # F2 a0 --> Hammer d0
-0.1872442367469107     , # F2 a1 --> Hammer d1
# 0.004653366641100383   , # F2 a2 --> Hammer d2
]).astype(np.float64)

bgl_coefficient_uncertainties = np.array([
 0.001074224733603991   , # g  a0 --> Hammer a0
 0.020436404332036355   , # g  a1 --> Hammer a1
 0.15751856111042434    , # g  a2 --> Hammer a2
 0.00030889008194163943 , # f  a0 --> Hammer b0
 0.019749982460876937   , # f  a1 --> Hammer b1
 0.18234034323164144    , # f  a2 --> Hammer b2
 #0.00005167468489451729 , # F1 a0 --> Hammer c0 <=== MISSING FROM HAMMER
 0.003158525549523684   , # F1 a1 --> Hammer c1
 0.10378051996527263    , # F1 a2 --> Hammer c2
 0.008389888682776101   , # F2 a0 --> Hammer d0
 0.28296856464624665    , # F2 a1 --> Hammer d1
# 0.23536654041924412    , # F2 a2 --> Hammer d2
]).astype(np.float64)


# correlation matrix (normalised covariance matrix)

import numpy as np
correlation_complete = np.array([
    
#     a0        a1        a2        b0        b1        b2        c0        c1        c2        d0        d1        d2
    [ 1.00000,  0.33659, -0.29501, -0.02177, -0.00252,  0.04887,  0.09463,  0.05693, -0.07992,  0.02260, -0.00535, -0.00158],
    [ 0.33659,  1.00000,  0.27320,  0.06180,  0.07523,  0.09627, -0.02265,  0.09332,  0.20195,  0.14681,  0.11106, -0.07422],
    [-0.29501,  0.27320,  1.00000,  0.03314,  0.05118,  0.10339, -0.08024, -0.01450,  0.10580,  0.08645, -0.02625, -0.16763],
    [-0.02177,  0.06180,  0.03314,  1.00000,  0.98181,  0.00631,  0.23178,  0.21528, -0.11545,  0.21363, -0.02327, -0.08270],
    [-0.00252,  0.07523,  0.05118,  0.98181,  1.00000,  0.14260,  0.23435,  0.23815, -0.10307,  0.18737,  0.03022, -0.05477],
    [ 0.04887,  0.09627,  0.10339,  0.00631,  0.14260,  1.00000, -0.15393, -0.14880,  0.03450, -0.10226, -0.05931, -0.04357],
    [ 0.09463, -0.02265, -0.08024,  0.23178,  0.23435, -0.15393,  1.00000,  0.83136, -0.32258,  0.42373,  0.41520,  0.02785],
    [ 0.05693,  0.09332, -0.01450,  0.21528,  0.23815, -0.14880,  0.83136,  1.00000,  0.24205,  0.67425,  0.74508,  0.11953],
    [-0.07992,  0.20195,  0.10580, -0.11545, -0.10307,  0.03450, -0.32258,  0.24205,  1.00000,  0.46548,  0.51003,  0.11613],
    [ 0.02260,  0.14681,  0.08645,  0.21363,  0.18737, -0.10226,  0.42373,  0.67425,  0.46548,  1.00000,  0.22556, -0.32435],
    [-0.00535,  0.11106, -0.02625, -0.02327,  0.03022, -0.05931,  0.41520,  0.74508,  0.51003,  0.22556,  1.00000,  0.47691],
    [-0.00158, -0.07422, -0.16763, -0.08270, -0.05477, -0.04357,  0.02785,  0.11953,  0.11613, -0.32435,  0.47691,  1.00000]
]).astype(np.float64)


#Missing in Hammer: c0, d2 
# c0 missing concluded from: https://gitlab.com/mpapucci/Hammer/-/blob/master/src/FormFactors/BGL/FFBtoDstarBGL.cc
# together with the referenced paper: https://arxiv.org/pdf/1902.09553 a
# d2 missing not as clear... took it from: 
# https://github.com/ocerri/BPH_RDntuplizer/blob/a4b157f5a64473bf3db360019d55fa2217199015/plugins/HammerWeightsProducer.cc#L358 

corr  = np.array([
    
#     a0        a1        a2        b0        b1        b2         c1        c2        d0        d1        
    [ 1.00000,  0.33659, -0.29501, -0.02177, -0.00252,  0.04887,   0.05693, -0.07992,  0.02260, -0.00535] ,
    [ 0.33659,  1.00000,  0.27320,  0.06180,  0.07523,  0.09627,   0.09332,  0.20195,  0.14681,  0.11106] ,
    [-0.29501,  0.27320,  1.00000,  0.03314,  0.05118,  0.10339,  -0.01450,  0.10580,  0.08645, -0.02625] ,
    [-0.02177,  0.06180,  0.03314,  1.00000,  0.98181,  0.00631,   0.21528, -0.11545,  0.21363, -0.02327] ,
    [-0.00252,  0.07523,  0.05118,  0.98181,  1.00000,  0.14260,   0.23815, -0.10307,  0.18737,  0.03022] ,
    [ 0.04887,  0.09627,  0.10339,  0.00631,  0.14260,  1.00000,  -0.14880,  0.03450, -0.10226, -0.05931] ,
    [ 0.05693,  0.09332, -0.01450,  0.21528,  0.23815, -0.14880,   1.00000,  0.24205,  0.67425,  0.74508] ,
    [-0.07992,  0.20195,  0.10580, -0.11545, -0.10307,  0.03450,   0.24205,  1.00000,  0.46548,  0.51003] ,
    [ 0.02260,  0.14681,  0.08645,  0.21363,  0.18737, -0.10226,   0.67425,  0.46548,  1.00000,  0.22556] ,
    [-0.00535,  0.11106, -0.02625, -0.02327,  0.03022, -0.05931,   0.74508,  0.51003,  0.22556,  1.00000] 
]).astype(np.float64)


get_ff_variations(names_bgl_coefficients_hammer, bgl_coefficients, bgl_coefficient_uncertainties, corr, "bgl_vector")
