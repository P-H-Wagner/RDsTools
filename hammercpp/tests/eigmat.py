import numpy as np
np.set_printoptions(precision=10, suppress=True, linewidth=200)

###################################################### 
# Using https://arxiv.org/pdf/2111.09849, Table 48   #
# for central values, errors, and correlation matrix #
###################################################### 

# define BCL parameter names
parName = ['ap0','ap1', 'ap2', 'ap3','a00','a01','a02']
# define BCL central values 
parVal = np.array([0.374, -0.672, 0.07, 1.34, 0.2203, 0.089, 0.24])
# define BCL uncertainties 
parSig = np.array([0.12, 0.64, 0.31, 0.52, 0.68, 0.57, 0.23])

# define correlation matrix 
corrM =  np.array([
[1.0, 0.2471, -0.1715, -0.2396, 0.6445, 0.3791, 0.2857],
[0.2471, 1.0, 0.4198, 0.1724, 0.4626, 0.8183, 0.7948],
[-0.1715, 0.4198, 1.0, 0.8136, 0.3804, 0.7293, 0.7481],
[-0.2396, 0.1724, 0.8136, 1.0, 0.2823, 0.5120, 0.5529],
[0.6445, 0.4626, 0.3804, 0.2823, 1.0, 0.6570, 0.4837],
[0.3791, 0.8183, 0.7293, 0.5120, 0.6570, 1.0, 0.9220],
[0.2857, 0.7948, 0.7481, 0.5529, 0.4837, 0.9220, 1.0]
])

# get covariance mattrix from correlation matrix
covM = np.zeros_like(corrM)

for i in range(len(parName)): 
  for j in range(len(parName)): 
    covM[i][j] = corrM[i][j] * parSig[i] * parSig[j]

# diagonalize the covariance matrix and get eigenvalues and eigenvectors
# eigVec[:,i] is the i-th eigenvector - column-wise!
# eigenvectors are normalized by linalg
eigVal, eigVec = np.linalg.eig(covM)

# take the sqrt of the eigenvalues
eigSig = np.sqrt(eigVal)

#now fill the eigmat
eigmat = np.zeros_like(covM)

for i in range(len(parName)): 
  for j in range(len(parName)): 
    #also here, keep the scaled eigenvectors columnwise!
    eigmat[j][i] = eigSig[i] * eigVec[j][i] 

print("===> EIGMAT (SCALED EIGENVECTORS ARE COLUMN_WISE!)")
print(eigmat)
