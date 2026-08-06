import numpy as np

# LHCb-PAPER-2018-050/
# values taken from https://lhcbproject.web.cern.ch/lhcbproject/Publications/LHCbProjectPublic/Directory_LHCb-PAPER-2018-050/Table_4.pdf

bins = [

4,
12, #merge all bins in one fat bin from [4,12] and take the average
13,
14,
16,
18,
20,
45 #extend the last bin to 45 
]

values = [

#0.324, #bin [4,5]
#0.281,
#0.257,
#0.245,
#0.227,
#0.210,
#0.194,
#0.191, ... bin [11,12]
todo, #artificial bin [4,12]

0.172,
0.159,
0.165,
0.136,
0.126,
0.109

]
