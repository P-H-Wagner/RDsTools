import numpy as np

# LHCb-PAPER-2018-050/
# values taken from https://lhcbproject.web.cern.ch/lhcbproject/Publications/LHCbProjectPublic/Directory_LHCb-PAPER-2018-050/Table_4.pdf

bins = [

4,
5,
6,
7,
8,
9,
10,
11,
12, 
13,
14,
15, #artifical new bin [15,16] (instead of [14,16]
16,
18,
20, #extend to 45 to match CMS binning
23,
26,
29,
34,
45, 
]

values = [

0.324,
0.281,
0.257,
0.245,
0.227,
0.210,
0.194,
0.191,
0.172,
0.159,
0.165, 
0.165, # double for bin [15,16]
0.136,
0.126,
0.109, #
0.109, #
0.109, # all of those are copied for bins > 20
0.109, #
0.109  #

]

unc = [

np.sqrt(0.001**2 + 0.025**2),
np.sqrt(0.001**2 + 0.018**2),
np.sqrt(0.001**2 + 0.017**2),
np.sqrt(0.001**2 + 0.017**2),
np.sqrt(0.001**2 + 0.015**2),
np.sqrt(0.001**2 + 0.015**2),
np.sqrt(0.001**2 + 0.013**2),
np.sqrt(0.001**2 + 0.014**2),
np.sqrt(0.001**2 + 0.013**2),
np.sqrt(0.001**2 + 0.012**2),
np.sqrt(0.001**2 + 0.012**2),
np.sqrt(0.001**2 + 0.012**2),
np.sqrt(0.001**2 + 0.010**2),
np.sqrt(0.001**2 + 0.010**2),
np.sqrt(0.001**2 + 0.009**2),
np.sqrt(0.001**2 + 0.009**2),
np.sqrt(0.001**2 + 0.009**2),
np.sqrt(0.001**2 + 0.009**2),
np.sqrt(0.001**2 + 0.009**2)


]
