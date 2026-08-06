import numpy as np

# These values are taken from https://www.hepdata.net/record/ins3118755
# resp. https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/BPH-21-007/index.html
# resp. BPH-21-007

bins = [

4, #extend the first bin to 4 (not starting from 11)
12,
13,
14,
#15, #merge to one bin [14,16]
16,
18,
20,
#23,
#26,
#29,
#34,
45, #merge to one bin [20,45]
]
#higher bins in pt are omitted

values = [

0.2648,
0.241,
0.2348,
#0.2326,
#0.2287,
todo,  #take the average of the above
0.2229,
0.2227,
0.2237,
0.2199,
0.2207,
0.2193,
0.2251,
]


