import numpy as np

# These values are taken from https://www.hepdata.net/record/ins3118755
# resp. https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/BPH-21-007/index.html
# resp. BPH-21-007

bins = [

#4, #extend the first bin to 4 (not starting from 11) to match lhcb
#5,
#6,
#7,
8,
9,
10,
11,
12,
13,
14,
15, 
16,
18,
20,
23,
26,
29,
34,
45, 
]

values = [

#0.2648,#
#0.2648,#
#0.2648,#
#0.2648,# all of those are duplicated for the first bin [4,12]
0.2648,#
0.2648,#
0.2648,#
0.2648,#
0.241,
0.2348,
0.2326,
0.2287,
0.2229,
0.2227,
0.2237,
0.2199,
0.2207,
0.2193,
0.2251,
]

unc = [

np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0056**2 + 0.0082**2),
np.sqrt(0.0039**2 + 0.0065**2),
np.sqrt(0.0031**2 + 0.0056**2),
np.sqrt(0.0028**2 + 0.006 **2),
np.sqrt(0.0018**2 + 0.0059**2),
np.sqrt(0.0018**2 + 0.0062**2),
np.sqrt(0.0016**2 + 0.0065**2),
np.sqrt(0.0018**2 + 0.0058**2),
np.sqrt(0.002 **2 + 0.007 **2),
np.sqrt(0.002 **2 + 0.0051**2),
np.sqrt(0.002 **2 + 0.0061**2),
np.sqrt(0.0029**2 + 0.0059**2),

]
