import ROOT
import argparse
import os
import sys
import yaml
from datetime import datetime

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getRdf
from helper import * 

import numpy as np

now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")



# get chain and df of signal sample
ch, df = getRdf(f"skimmed/{sig_cons_25[0]}")

#selection
selection = base_wout_tv_25 + " && (mu7_ip4)"
selection += " && (gen_sig == 0)" #fit only dsmu

#fit the peak! :)
s, h  = getSigma(df, "phiPi_m", selection, dt)


