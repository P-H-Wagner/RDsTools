import uproot
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot(file, name):

  h = file["hist_scale_factor"]
  
  # get matrix (it is symmetric!)
  matrix = h.values()
  var    = h.variances()
  err    = np.sqrt(var)
  
  pt_edges  = h.axis(0).edges()
  dxy_edges = h.axis(1).edges() 
  
  #save into txt format to use later in skimmer
  np.savetxt(f"trigger_sf_{name}.txt"    , matrix   )
  np.savetxt(f"trigger_sf_err_{name}.txt", err      )
  np.savetxt(f"pt_edges_{name}.txt"      , pt_edges )
  np.savetxt(f"dxy_edges_{name}.txt"     , dxy_edges)
 
  plt.figure() 
  xtickpos = []
  xlabel   = []
  for i in range(len(pt_edges)-1):
      xtickpos.append(i)
      xlabel  .append(f"[{pt_edges[i]},{pt_edges[i+1]}]")
  
  ytickpos = []
  ylabel   = []
  for i in range(len(dxy_edges)-1):
      ytickpos.append(i)
      ylabel.append(f"[{dxy_edges[i]},{dxy_edges[i+1]}]")
  
  
  #imshow takes y first, then x, ANNOYING!!! --> transpose
  plt.imshow(matrix.T, origin="lower", aspect="auto", cmap="viridis")
  plt.colorbar(label="Scale factor")
  
  #bin number is always one smaller than edges :D
  for i in range(len(pt_edges)-1):
      for j in range(len(dxy_edges)-1):
          plt.text(i,j, f"{matrix[i,j]:.4f}\n" + rf"$\pm${err[i,j]:.4f}", ha="center", va="center", color="white", fontsize=5)
  
  plt.xticks(xtickpos,xlabel, rotation = -20)
  plt.yticks(ytickpos,ylabel)
  
  xlabel = "mu_pt"
  ylabel = "dxy_sig"
  
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.tight_layout()
  
  plt.savefig(f"scale_factors_{name}.pdf")
  plt.close()


file = uproot.open("scale_factors_fullBPark.root")
plot(file, "fullBPark")
file = uproot.open("scale_factors_D1.root")
plot(file, "D1")


