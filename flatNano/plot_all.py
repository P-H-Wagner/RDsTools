
###############################
## Plotting script           ##
###############################

import ROOT
import math
import numpy as np






# Turn off graphics display
ROOT.gROOT.SetBatch(True)

import argparse
parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help') 

parser.add_argument('filename')   
args = parser.parse_args()

#input files
files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.filename}/*"

#chain them
chain = ROOT.TChain("tree")
chain.Add(files)

#get branche names
names = [branch.GetName() for branch in chain.GetListOfBranches()]

#get number of entries
nEntries = chain.GetEntries("bs_pt_reco_1 == bs_pt_reco_1");
nEntriesMu  = chain.GetEntries("(bs_pt_reco_1 == bs_pt_reco_1) && (sig == 0 || sig == 5)");
nEntriesTau = chain.GetEntries("(bs_pt_reco_1 == bs_pt_reco_1) && (sig == 1 || sig == 6)");

#colors to multiVar
colors = [
    ROOT.kBlack,
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kGreen,
    ROOT.kOrange,
    ROOT.kMagenta,
    ROOT.kCyan,
    ROOT.kYellow,
    ROOT.kGray,
    ROOT.kViolet,
    ROOT.kTeal,
    ROOT.kSpring,
    ROOT.kAzure,
    ROOT.kPink,
    ROOT.kSpring + 6,
    ROOT.kTeal + 6,
    ROOT.kAzure + 6,
    ROOT.kPink + 6
]

labelVar={"q2": "Q^{2} ", 
        "m2_miss": r"m^{2}_#mathrm{miss}",
        "cosMuW": r"cos(#mu W)",
        "cosPhiDs": r"cos(#phi D_{s})",
        "cosPiDs": r"cos(#pi D_{s})",
}

labelMethod={
        "coll": r"Collinear",
        #"lhcb": r"LHCb z",
        "lhcb_alt": r"LHCb xyz",
        "reco_1": r"Math Sol. 1", 
        "reco_2": r"Math Sol. 2",
        "gen"   : r"Gen",
        "Coll": r"Collinear",
        #"lhcb": r"LHCb z",
        "LhcbAlt": r"LHCb xyz",
        "Reco1": r"Math Sol. 1", 
        "Reco2": r"Math Sol. 2",
        "Gen"   : r"Gen",

}

#signal indices

#   Ds mu  Ds tau  Ds* mu  Ds* tau
sigs = [0,      1,      5,      6]

#signal names
sigNames = {"0": r"D_{s} #mu","1": r"D_{s} #tau", "5":r"D_{s}* #mu", "6":r"D_{s}* #tau", "gen": "Gen"}

#disable title and stats box
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

################################################################################################
#function which plots the desired variables into one canvas 

def multiVar(variables, bins, begin, end, name, norm = True, xLabel = "[GeV]", yLabel = "a.u.",color = colors,log = False):

  """
  variables        = list of strings, the variables to plot
  begin, en d      = floats, begin and end of range
  bins             = int, nr of bins
  name             = str, the name under which the plot should be saved (name.pdf)
  norm             = bool, tells if we normalize the histograms
  xLabel, yLabel   = strings, axis labels
  """

  #holds histograms
  histos = {}
  addresses = {}

  #define histograms
  for i,var in enumerate(variables):
    histos[var] = ROOT.TH1F(var,var,bins,begin,end)
    histos[var].SetLineColor(color[i])
    histos[var].SetLineWidth(2)

  histos["best"] = ROOT.TH1F("best","best",bins,begin,end)
  histos["best"].SetLineColor(ROOT.kOrange + 7)
  histos["best"].SetLineWidth(2)

  #fill all histograms

  print("Start filling " + name + " histogram; Events to fill: " + str(nEntries))


  reco2Solved = 0
  sign = 1
  for i in range(nEntries):
    if (i%20000 == 0):
      print(" ==> Processing " + str(i) + "th Event")
    chain.GetEntry(i)

    reco1Val = np.nan
    reco2Val = np.nan
    genVal   = np.nan
 
    for var in variables:
      if "cosMuW" in var: 
        sign = -1
      if (getattr(chain,"sig") == 0):
        #if not (math.isnan(getattr(chain,"bs_pt_reco_1"))):
          histos[var].Fill(sign*getattr(chain,var))

          #Get the recos and gen methods
          if ("reco_1" in var) or ("Reco1" in var): reco1Val = getattr(chain,var)
          if ("reco_2" in var) or ("Reco2" in var): reco2Val = getattr(chain,var)
          if ("gen" in var):                        genVal   = getattr(chain,var)
       
      if ("reco_1" in var) or ("Reco1" in var): reco1Str = var
      if ("reco_2" in var) or ("Reco2" in var): reco2Str = var
      if ("gen" in var) or ("Gen" in var): genStr = var 
  
  reco1Solved = chain.GetEntries("abs(" + reco1Str + "-" + genStr + ") < abs(" + reco2Str + "-" +genStr + ") && (bs_pt_reco_1 == bs_pt_reco_1)")
  print(f"Reco 1 was the best solution {reco1Solved/nEntries} of the time, reco2 solved it {1-reco1Solved/nEntries} of the time")

  """
    if (reco1Val == reco1Val) and (reco2Val == reco2Val):

      # get the better solution (closer to gen)
      recos = [reco1Val, reco2Val]
      diffs = [abs(reco1Val - genVal), abs(reco2Val - genVal)]
   
      # argmin is either 0 (reco1) or 1 (reco2). We want to count:
      reco2Solved += np.argmin(diffs) # increase counter if reco2 solved it
  """
  #print(f"Reco 1 was the best solution {1 - reco2Solved/nEntries} of the time, reco2 solved it {reco2Solved/nEntries} of the time")
  
  # Fillt histo with the better solution 
  for var in variables:
    if ("reco_1" in var) or ("Reco1" in var):
      hDummy1 = histos[var].Clone()
      #reco1 is the better solution 95% of the time
      histos["best"].Add(hDummy1, reco1Solved/nEntries * 1/histos[var].Integral())
    if ("reco_2" in var) or ("Reco2" in var):
      hDummy2 = histos[var].Clone()
      #reco2 is the better soltuion 5% of the time
      histos["best"].Add(hDummy2, 1-reco1Solved/nEntries * 1/histos[var].Integral())

 
  #scale histograms
  if norm:
    for var in histos.keys():
      print(var)
      try: histos[var].Scale(1/histos[var].Integral())
      except: continue
  #get maximumm of y axis 
  yMax = max([histos[var].GetMaximum() for var in histos.keys()]) 

  #Draw histos
  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  canvas.cd()

  if log: canvas.SetLogy()

  if len(variables) > 2:
    legend = ROOT.TLegend(0.25, 0.75, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
    legend.SetNColumns(3)

  else:
    legend = ROOT.TLegend(0.55, 0.8, 0.9, 0.9);
    legend.SetNColumns(2)


  legend.SetTextSize(0.03)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)

  for i,var in enumerate(variables):
    #if ("reco_1" in var) or ("Reco1" in var) or ("reco_2" in var) or ("Reco2" in var): continue #comment this out to not plot recos
    #multiVar
    if i == 0:
      histos[var].SetMaximum(yMax*1.4)
      histos[var].GetYaxis().SetTitle(yLabel)
      histos[var].GetXaxis().SetTitle(xLabel)
      histos[var].GetXaxis().SetTitleOffset(1.3)

      histos[var].Draw("HIST")
    else:
      histos[var].Draw("HIST SAME")


    # legend
    varLegend = r""

    #loop over methods
    for iName in labelMethod.keys():
      if iName in var: varLegend += labelMethod[iName]
    if ("lhcb" in var and "lhcb_alt" not in var) or ("Lhcb" in var and "LhcbAlt" not in var): varLegend += "LHCb z" #special bc lhcb also in lhcb_alt

    if "cosPiK" in var and "Gen" not in var: varLegend += "Reco"

    legend.AddEntry(histos[var], varLegend, "l")

   
  histos["best"].Draw("HIST SAME") ##CHANGE
  legend.AddEntry(histos["best"], "Weighted Math Sol.", "l")
  legend.SetTextSize(0.04)
  legend.Draw("SAME")

  #saving
  canvas.SaveAs("./plots/" + name + "_weighted.pdf")



################################################################################################
#function which plots the desired signals into one canvas 
def multiSig(var, bins, begin, end, name, signals = sigs, norm = True, xLabel = r" Q^{2} [GeV]", yLabel = "a.u.", color = colors, title = ""):

  """
  variables        = list of strings, the variables to plot, the secon one is gen
  begin, en d      = floats, begin and end of range
  bins             = int, nr of bins
  name             = str, the name under which the plot should be saved (name.pdf)
  norm             = bool, tells if we normalize the histograms
  xLabel, yLabel   = strings, axis labels
  """

  #holds histograms
  histos = {}
  addresses = {}

  toPlotSig = [str(sig) for sig in signals]
  toPlotGen = ["gen_" + str(sig) for sig in signals]

  plotGen = False

  if len(var) > 1: 
    plotGen = True
    toPlot = toPlotSig + toPlotGen

  else: toPlot = toPlotSig 

  print(toPlot) 
  #define histograms
  for i,sig in enumerate(toPlot):
  
    histos[str(sig)] = ROOT.TH1F(str(sig),str(sig),bins,begin,end)
    histos[str(sig)].SetLineColor(color[i])
    histos[str(sig)].SetLineWidth(2)
 
  if "reco" in var[0]:
    for sig in toPlotSig:

      histos["reco2_"+sig] = ROOT.TH1F("reco"+sig,"reco"+sig,bins,begin,end)
      histos["reco2_"+sig].SetLineColor(ROOT.kViolet - 9)
      histos["reco2_"+sig].SetLineWidth(2)

      histos["best_"+sig] = ROOT.TH1F("best"+sig,"best"+sig,bins,begin,end)
      histos["best_"+sig].SetLineColor(ROOT.kViolet - 9)
      histos["best_"+sig].SetLineWidth(2)

  print(histos)
  #fill all histograms

  print("Start filling " + name + " histogram; Events to fill: " + str(nEntries))



  for i in range(nEntries):
    if (i%20000 == 0):
      print(" ==> Processing " + str(i) + "th Event")
    chain.GetEntry(i)
    for sig in signals:
      if (sig == getattr(chain,"sig")):
        histos[str(sig)].Fill(getattr(chain,var[0]))

        if "reco" in var[0]:
          histos["reco2_" + str(sig)].Fill(getattr(chain,var[2]))


        if plotGen:
          histos["gen_" + str(sig)].Fill(getattr(chain,var[1]))
   
  if "reco" in var[0]:

    for sig in toPlotSig:

      histos["best_" + sig].Add(histos[sig], 2/3 * 1/histos[sig].Integral())    
      histos["best_" + sig].Add(histos["reco2_" + sig], 1/3 * 1/histos["reco2_" + sig].Integral())    

 
  #scale histograms
  if norm:
    for sig in toPlot:
      histos[str(sig)].Scale(1/histos[str(sig)].Integral())

  #get maximumm of y axis 
  yMax = max([histos[str(sig)].GetMaximum() for sig in signals]) 

 
  

  if plotGen:

    if "reco" in var[0]:
      print("I AM HERE!")
      for key1,key2 in zip(toPlotSig,toPlotGen):
        #Draw histos
        canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
        canvas.cd()
        legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetNColumns(2)
        legend.SetTextSize(0.04)
        
        histos["best_"+key1].SetMaximum(yMax*1.2)
        histos["best_"+key1].GetXaxis().SetTitle(xLabel)
        histos["best_"+key1].GetYaxis().SetTitle(yLabel)
        histos["best_"+key1].Draw("HIST")
     
        histos[key2].SetMaximum(yMax*1.2)
        histos[key2].GetXaxis().SetTitle(xLabel)
        histos[key2].GetYaxis().SetTitle(yLabel)
        histos[key2].Draw("HIST SAME")
     
        histos["best_"+key1].GetXaxis().SetTitleOffset(1.3)
        histos[key2].GetXaxis().SetTitleOffset(1.3)
  
        legend.AddEntry(histos["best_"+key1], sigNames[key1] , "l")
        legend.AddEntry(histos[key2], "Gen", "l")
  
        legend.Draw("SAME")
        canvas.SaveAs(f"./plots/{name}_{key1}_weighted.pdf")
  


    else:
      for key1,key2 in zip(toPlotSig,toPlotGen):
  
        #Draw histos
        canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
        canvas.cd()
        legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetNColumns(2)
        legend.SetTextSize(0.04)
        
        histos[key1].SetMaximum(yMax*1.2)
        histos[key1].GetXaxis().SetTitle(xLabel)
        histos[key1].GetYaxis().SetTitle(yLabel)
        histos[key1].Draw("HIST")
     
        histos[key2].SetMaximum(yMax*1.2)
        histos[key2].GetXaxis().SetTitle(xLabel)
        histos[key2].GetYaxis().SetTitle(yLabel)
        histos[key2].Draw("HIST SAME")
     
        histos[key1].GetXaxis().SetTitleOffset(1.3)
        histos[key2].GetXaxis().SetTitleOffset(1.3)
  
        legend.AddEntry(histos[key1], sigNames[key1] , "l")
        legend.AddEntry(histos[key2], "Gen", "l")
  
        legend.Draw("SAME")
        canvas.SaveAs(f"./plots/{name}_{key1}.pdf")

  else:

    #Draw histos
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    canvas.cd()
    legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetNColumns(2)
    legend.SetTextSize(0.04)


    for i,sig in enumerate(toPlot):



  
      #multiVar
      if i == 0:
        if title != "": histos[str(sig)].SetTitle(title)
        histos[str(sig)].SetMaximum(yMax*1.2)
        histos[str(sig)].GetXaxis().SetTitle(xLabel)
        histos[str(sig)].GetYaxis().SetTitle(yLabel)
        histos[str(sig)].Draw("HIST")
      else:
        histos[str(sig)].Draw("HIST SAME")
  
      #legend
      legend.AddEntry(histos[str(sig)], sigNames[str(sig)], "l")
  
    legend.Draw("SAME")
  
    #saving
    canvas.SaveAs("./plots/" + name + ".pdf")



def compareFamilies(var, bins, begin, end, name, gen1 = [0,5], gen2 = [1,6], genName = [r"#mu - signals", r"#tau - signals" ], norm = True, xLabel = r" Q^{2} [GeV]", yLabel = "a.u.",log = False):

  """
  variables        = list of strings, the variables to plot
  begin, en d      = floats, begin and end of range
  bins             = int, nr of bins
  name             = str, the name under which the plot should be saved (name.pdf)
  norm             = bool, tells if we normalize the histograms
  xLabel, yLabel   = strings, axis labels
  """

  #holds histograms
  histos = {}
  addresses = {}

  generations = [1,2]
  colors = [ROOT.kAzure, ROOT.kPink] 
 
  #define histograms
  for i,sig in enumerate(generations):

    if len(var) > 1:
      
      histos[str(sig)+"reco1"] = ROOT.TH1F(str(sig)+"reco1",str(sig)+"reco1",bins,begin,end)
      histos[str(sig)+"reco2"] = ROOT.TH1F(str(sig)+"reco2",str(sig)+"reco2",bins,begin,end)
      histos[str(sig)+"best"]  = ROOT.TH1F(str(sig)+"best",str(sig)+"best",bins,begin,end)
      histos[str(sig)+"best"].SetLineColor(colors[i])
      histos[str(sig)+"best"].SetLineWidth(2)

    else:

      histos[str(sig)] = ROOT.TH1F(str(sig),str(sig),bins,begin,end)
      histos[str(sig)].SetLineColor(colors[i])
      histos[str(sig)].SetLineWidth(2)
  

  #fill all histograms

  print("Start filling " + name + " histogram; Events to fill: " + str(nEntries))

  sign = 1
  if "cosMuW" in var[0]: sign = -1
  
  if len(var) == 1:
    for i in range(1000):
      if (i%20000 == 0):
        print(" ==> Processing " + str(i) + "th Event")
      chain.GetEntry(i)

           
      if getattr(chain,"sig") in gen1:
        histos["1"].Fill(sign*getattr(chain,var[0]))
      else:
        histos["2"].Fill(sign*getattr(chain,var[0]))
 
      #scale histograms
    if norm:
      for sig in generations:
        histos[str(sig)].Scale(1/histos[str(sig)].Integral())
 
  else:
    #reco case
    for i in range(1000):
      if (i%20000 == 0):
        print(" ==> Processing " + str(i) + "th Event")

      chain.GetEntry(i)

      if getattr(chain,"sig") in gen1:
          #first generation, fill with reco 1,2
          histos["1reco1"].Fill(sign*getattr(chain,var[0]))
          histos["1reco2"].Fill(sign*getattr(chain,var[1]))
          #histos["gen1"].Fill(sign*getattr(chain,var[2]))

      else:
          #second generation gill with reco 1,2
          histos["2reco1"].Fill(sign*getattr(chain,var[0])) 
          histos["2reco2"].Fill(sign*getattr(chain,var[1])) 
          #histos["gen2"].Fill(sign*getattr(chain,var[2]))

    #out of event loop
    # get the variable names
    if ("reco_1" in var[0]) or ("Reco1" in var[0]): reco1Str = var[0]
    if ("reco_2" in var[1]) or ("Reco2" in var[1]): reco2Str = var[1]
    if ("gen" in var[2]) or ("Gen" in var[2]): genStr = var[2] 


    reco1Solved = {}
    nEntriesGen= {}
    if gen1 == [0,5]: #mu vs tau, mu is gneration 1!
      reco1Solved["1"] = chain.GetEntries("abs(" + reco1Str + "-" + genStr + ") < abs(" + reco2Str + "-" +genStr + ") && (bs_pt_reco_1 == bs_pt_reco_1) && ((sig == 0) || (sig == 5)) ") #mu
      reco1Solved["2"] = chain.GetEntries("abs(" + reco1Str + "-" + genStr + ") < abs(" + reco2Str + "-" +genStr + ") && (bs_pt_reco_1 == bs_pt_reco_1) && ((sig == 1)  || (sig == 6)) ") #tau

      nEntriesGen["1"] = nEntriesMu
      nEntriesGen["2"] = nEntriesTau
     
      #print(f"gen1: Reco 1 was the best solution {reco1Solved[\"1\"]/nEntries} of the time, reco2 solved it {1-reco1Solved[\"1\"]/nEntries} of the time")
      #print(f"gen2: Reco 1 was the best solution {reco1Solved[\"2\"]/nEntries} of the time, reco2 solved it {1-reco1Solved[\"2\"]/nEntries} of the time")

    hDummy1 = {}
    hDummy2 = {}
    for sig in generations:
      #loop over generation 1
      for ivar in var[0:2]:
        #loop over  reco1 and reco2
        if ("reco_1" in ivar) or ("Reco1" in ivar):
          # reco 1 case
          hDummy1[str(sig)] = histos[str(sig)+"reco1"].Clone()
          histos[str(sig)+"best"].Add(hDummy1[str(sig)], reco1Solved[str(sig)]/nEntriesGen[str(sig)] * 1/hDummy1[str(sig)].Integral())
  
        if ("reco_2" in ivar) or ("Reco2" in ivar):
          # reco2 case
          hDummy2[str(sig)] = histos[str(sig)+"reco2"].Clone()
          histos[str(sig)+"best"].Add(hDummy2[str(sig)],(1-reco1Solved[str(sig)]/nEntriesGen[str(sig)]) * 1/hDummy2[str(sig)].Integral())
  
    

  if len(var) == 1:

    #calculate gini coeff
    gini_tot = histos["1"].Integral() + histos["2"].Integral()
    area_1 = 0
    area_2 = 0
  
    for i in range(bins):
      #for the first histogram
      if (histos["1"].GetBinContent(i+1) > histos["2"].GetBinContent(i+1)):
        area_1 += histos["1"].GetBinContent(i+1)-histos["2"].GetBinContent(i+1) 
        area_2 += histos["2"].GetBinContent(i+1) 
  
      else: 
        area_1 += histos["1"].GetBinContent(i+1) 
        area_2 += histos["2"].GetBinContent(i+1)-histos["1"].GetBinContent(i+1) 

  else:

    #calculate gini coeff
    gini_tot = histos["1best"].Integral() + histos["2best"].Integral()
    area_1 = 0
    area_2 = 0
  
    for i in range(bins):
      #for the first histogram
      if (histos["1best"].GetBinContent(i+1) > histos["2best"].GetBinContent(i+1)):
        area_1 += histos["1best"].GetBinContent(i+1)-histos["2best"].GetBinContent(i+1) 
        area_2 += histos["2best"].GetBinContent(i+1) 
  
      else: 
        area_1 += histos["1best"].GetBinContent(i+1) 
        area_2 += histos["2best"].GetBinContent(i+1)-histos["1best"].GetBinContent(i+1) 
  
  gini_coeff = round((area_1 + area_2)/gini_tot,2)
  giniText = ROOT.TPaveText(0.2, 0.8, 0.5, 0.9, 'nbNDC')
  giniText.AddText(f'Gini Coeff = {gini_coeff}')
  giniText.SetTextFont(42)
  giniText.SetTextSize(0.04)
  giniText.SetFillColor(0)
  giniText.SetFillStyle(0)
  giniText.SetBorderSize(0)

  if len(var) ==1:
    #get maximumm of y axis 
    yMax = max([histos[str(sig)].GetMaximum() for sig in generations]) 

  else:
    yMax = max([histos[str(sig)+"best"].GetMaximum() for sig in generations])

  #Draw histos
  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  canvas.cd()
  if log: canvas.SetLogy()

  legend = ROOT.TLegend(0.5, 0.8, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
  legend.SetTextSize(0.03)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  legend.SetNColumns(2)
  legend.SetTextSize(0.04)

  for i,sig in enumerate(generations):
    if len(var) == 1:
      #multiVar
      if i == 0:
        histos[str(sig)].SetMaximum(yMax*1.2)
        histos[str(sig)].GetYaxis().SetTitle(yLabel)
        histos[str(sig)].GetXaxis().SetTitle(xLabel)
        histos[str(sig)].Draw("HIST")
        histos[str(sig)].GetXaxis().SetTitleOffset(1.3)
  
      else:
        histos[str(sig)].Draw("HIST SAME")
  
      #legend
      legend.AddEntry(histos[str(sig)], genName[i], "l")

    else:
      #multiVar
      if i == 0:
        histos[str(sig)+"best"].SetMaximum(yMax*1.2)
        histos[str(sig)+"best"].GetYaxis().SetTitle(yLabel)
        histos[str(sig)+"best"].GetXaxis().SetTitle(xLabel)
        histos[str(sig)+"best"].Draw("HIST")
        histos[str(sig)+"best"].GetXaxis().SetTitleOffset(1.3)
  
      else:
        histos[str(sig)+"best"].Draw("HIST SAME")
  
      #legend
      legend.AddEntry(histos[str(sig)+"best"], genName[i], "l")


  legend.Draw("SAME")
  giniText.Draw("EP SAME")
  #saving
  canvas.SaveAs("./plots/" + name + ".pdf")


################################################################################################
#function which plots the desired signals into one canvas 
def testFit(variables, bins, begin, end, name, signals = sigs, norm = True, xLabel = r" p_{T}(#mu) [GeV]", yLabel = "a.u.", color = colors, zoom = False):

  """
  variables        = [before,after,gen]
  begin, end       = floats, begin and end of range
  bins             = int, nr of bins
  name             = str, the name under which the plot should be saved (name.pdf)
  norm             = bool, tells if we normalize the histograms
  xLabel, yLabel   = strings, axis labels
  """

  #holds histograms
  histos = {}
  addresses = {}

  fitTestNames = {"before": "Before Fit", "after": "After Fit", "gen": "Gen"}
  
  #define histograms
  for i,key in enumerate(["before","after","gen"]):
    histos[key] = ROOT.TH1F(key,key,bins,begin,end)
    histos[key].SetLineColor(color[i])
    histos[key].SetLineWidth(2)

  #fill all histograms

  print("Start filling " + key + " histogram; Events to fill: " + str(nEntries))

  for i in range(nEntries):
    if (i%20000 == 0):
      print(" ==> Processing " + str(i) + "th Event")
    chain.GetEntry(i)
    histos["before"].Fill(getattr(chain,variables[0]))
    histos["after"].Fill(getattr(chain,variables[1]))
    histos["gen"].Fill(getattr(chain,variables[2]))

  #scale histograms
  if norm:
    for key in ["before","after","gen"]:
      histos[key].Scale(1/histos[key].Integral())

  #get maximumm of y axis 
  yMax = max([histos[key].GetMaximum() for key in ["before","after","gen"] ]) 

  #Draw histos
  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 800)

  canvas.Divide(1, 2)
  canvas.cd(1)

  ROOT.gPad.SetPad(0, 0.25, 1, 1)

  legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
  legend.SetTextSize(0.03)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  legend.SetTextSize(0.04)

  for i,key in enumerate(["before","after","gen"]):

    #multiVar
    if i == 0:
      histos[key].SetMaximum(yMax*1.2)
      histos[key].GetYaxis().SetTitle(yLabel)
      #histos[key].GetXaxis().SetTitle(xLabel)
      #h1.GetXaxis().SetTickLength(0)
      histos[key].GetXaxis().SetLabelSize(0)

      histos[key].Draw("HIST")
    else:
      histos[key].Draw("HIST SAME")

    #legend
    legend.AddEntry(histos[key], fitTestNames[key], "l")

  legend.Draw("SAME")
 
  canvas.cd(2)

  ROOT.gPad.SetPad(0, 0.05, 1, 0.31)
  #ROOT.gPad.SetTopMargin(0.03) 
  ROOT.gPad.SetBottomMargin(0.27)
  ratio1 = histos["before"].Clone()
  ratio2 = histos["after"].Clone()



  for i in range(bins):

    if (histos["gen"].GetBinContent(i+1) != 0):
      print("im here")
      print(ratio1.GetBinContent(i+1) / histos["gen"].GetBinContent(i+1))
      #ratio1.Divide(histos["gen"])
      ratio1.SetBinContent(i+1,ratio1.GetBinContent(i+1) / histos["gen"].GetBinContent(i+1))
      #if (histos["after"].GetBinContent(i+1) != 0): 
      #ratio2.Divide(histos["gen"])
      ratio2.SetBinContent(i+1,ratio2.GetBinContent(i+1) / histos["gen"].GetBinContent(i+1)) 

    else: 
      ratio1.SetBinContent(i+1,0)
      ratio2.SetBinContent(i+1,0)

  ymax = 1.2 * max([ratio1.GetMaximum(),ratio2.GetMaximum()])
  ymin = 0.8 * min([ratio1.GetMinimum(),ratio2.GetMinimum()])

  ratio1.GetXaxis().SetTitle(xLabel)
  ratio1.GetXaxis().SetTitleSize(0.1)
  ratio1.GetXaxis().SetTickLength(0.1)
  ratio1.GetXaxis().SetLabelSize(0.08)


  ratio1.GetYaxis().SetRangeUser(ymin,ymax)
  if zoom: ratio1.GetYaxis().SetRangeUser(0.8,1.2)

  ratio1.GetYaxis().SetTitle(r"Reco / Gen")
  ratio1.GetYaxis().SetTitleSize(0.1)
  ratio1.GetYaxis().SetLabelSize(0.08)
  ratio1.GetYaxis().SetTitleOffset(0.5)
  ratio1.GetYaxis().SetNdivisions(8)
  ratio1.GetYaxis().CenterTitle() 
  ratio1.GetXaxis().SetTitleOffset(1.3)

  ratio1.SetMarkerStyle(8)
  ratio2.SetMarkerColor(color[0])

  ratio2.SetMarkerStyle(22)
  ratio2.SetMarkerColor(color[1])

  ratio1.Draw("EP")
  ratio2.Draw("SAME EP")

  zoomed = ""
  if zoom: zoomed = "_zoomed"

  #saving
  canvas.SaveAs(f"./plots/{name}{zoomed}.pdf")

def fitVtx(variables, bins, begin, end, name, labelVtx, signals = sigs, norm = True, xLabel = r" p_{T}(#mu) [GeV]", yLabel = "a.u.", color = colors):

  """
  variables        = [before,after,gen]
  begin, end       = floats, begin and end of range
  bins             = int, nr of bins
  name             = str, the name under which the plot should be saved (name.pdf)
  norm             = bool, tells if we normalize the histograms
  xLabel, yLabel   = strings, axis labels
  """

  #holds histograms
  histos = {}
  addresses = {}

  fitTestNames = {"fit": labelVtx, "gen": "Gen"}
  
  #define histograms
  for i,key in enumerate(["fit","gen"]):
    histos[key] = ROOT.TH1F(key,key,bins,begin,end)
    histos[key].SetLineColor(color[i])
    histos[key].SetLineWidth(2)

  #fill all histograms

  print("Start filling " + key + " histogram; Events to fill: " + str(nEntries))

  for i in range(nEntries):
    if (i%20000 == 0):
      print(" ==> Processing " + str(i) + "th Event")
    chain.GetEntry(i)
    histos["fit"].Fill(getattr(chain,variables[0]))
    histos["gen"].Fill(getattr(chain,variables[1]))

  #scale histograms
  if norm:
    for key in ["fit","gen"]:
      histos[key].Scale(1/histos[key].Integral())

  #get maximumm of y axis 
  yMax = max([histos[key].GetMaximum() for key in ["fit","gen"]]) 

  #Draw histos
  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)

  #canvas.Divide(1, 2)
  canvas.cd()

  #ROOT.gPad.SetPad(0, 0.25, 1, 1)

  legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
  legend.SetTextSize(0.03)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  legend.SetTextSize(0.04)

  for i,key in enumerate(["fit","gen"]):

    #multiVar
    if i == 0:
      histos[key].SetMaximum(yMax*1.2)
      histos[key].GetYaxis().SetTitle(yLabel)
      histos[key].GetXaxis().SetTitle(xLabel)
      histos[key].GetXaxis().SetTitleOffset(1.3)

      histos[key].Draw("HIST")
    else:
      histos[key].Draw("HIST SAME")

    #legend
    legend.AddEntry(histos[key], fitTestNames[key], "l")

  legend.Draw("SAME")

  """ 
  canvas.cd(2)

  ROOT.gPad.SetPad(0, 0.05, 1, 0.31)
  #ROOT.gPad.SetTopMargin(0.03) 
  ROOT.gPad.SetBottomMargin(0.22)
  ratio1 = histos["fit"].Clone()

  ymax = 1
  ymin = 1

  for i in range(bins):

    if (histos["gen"].GetBinContent(i+1) != 0):
      ratio1.SetBinContent(i+1,ratio1.GetBinContent(i+1) / histos["gen"].GetBinContent(i+1))

    else: 
      ratio1.SetBinContent(i+1,0)

    if ymax < ratio1.GetBinContent(i+1) : ymax = ratio1.GetBinContent(i+1)
    if ymin > ratio1.GetBinContent(i+1) : ymin = ratio1.GetBinContent(i+1)

  #print("found maximum",ratio1.GetMaximum())
  if ymin == 0: ymin = 0.001

  ratio1.GetXaxis().SetTitle(xLabel)
  ratio1.GetXaxis().SetTitleSize(0.1)
  ratio1.GetXaxis().SetTickLength(0.1)
  ratio1.GetXaxis().SetLabelSize(0.08)

  ratio1.GetYaxis().SetRangeUser(ymin,ymax)

  ratio1.GetYaxis().SetTitle(r"Reco / Gen")
  ratio1.GetYaxis().SetTitleSize(0.1)
  ratio1.GetYaxis().SetLabelSize(0.08)
  ratio1.GetYaxis().SetTitleOffset(0.5)
  ratio1.GetYaxis().SetNdivisions(8)
  ratio1.GetYaxis().CenterTitle() 


  ratio1.SetMarkerStyle(8)

  ratio1.Draw("EP")
  """

  
  #saving
  canvas.SaveAs(f"./plots/{name}.pdf")





###################################
## Call Plotting function        ##
###################################

#multiSig("lxy_bs" ,50,0,3   ,"lxy_bs" ,xLabel = "[cm]")
#multiSig("lxy_ds" ,50,0,1   ,"lxy_ds" ,xLabel = "[cm]")
#multiSig("lxy_phi",50,0,1.8 ,"lxy_phi",xLabel = "[cm]")

#multiSig("dxy_mu" ,50,0,0.3 ,"dxy_mu" ,xLabel = "[cm]")
#multiSig("dxy_pi" ,50,0,0.6 ,"dxy_pi" ,xLabel = "[cm]")
#multiSig("dxy_k1" ,50,0,0.4 ,"dxy_k1" ,xLabel = "[cm]")

#multiSig("dz_mu"  ,50,0,0.3 ,"dz_mu" ,xLabel = "[cm]")
#multiSig("dz_pi"  ,50,0,0.7 ,"dz_pi" ,xLabel = "[cm]")
#multiSig("dz_k1"  ,50,0,0.5 ,"dz_k1" ,xLabel = "[cm]")

#reco vs all signals and gen
recoColors = colors[0:4] + 4*[colors[5]]
#multiSig(["bs_pt_reco_1","bs_gen_pt","bs_pt_reco_2"]  ,50,0,100 ,"multiSigBsPtReco"      ,xLabel = "p_{T}(B_{s}) [GeV]", color = recoColors,title = "Math. Sol. on all signals")
#multiSig(["bs_pt_lhcb_alt","bs_gen_pt"]  ,50,0,100 ,"multiSigBsPtLhcbAlt" ,xLabel = "p_{T}(B_{s}) [GeV]", color = recoColors, title = "LHCb method on all signals")

# mu vs tau plots (default)
compareFamilies(["q2_lhcb_alt"]  ,50,0,12,"muVsTauQ2LhcbAlt" ,xLabel = r"Q^{2} [GeV]")
compareFamilies(["q2_lhcb"]  ,50,0,12,"muVsTauQ2Lhcb" ,xLabel = r"Q^{2} [GeV]")
compareFamilies(["q2_reco_1"]  ,50,0,12,"muVsTauQ2Reco1" ,xLabel = r"Q^{2} [GeV]")
compareFamilies(["q2_reco_2"]  ,50,0,12,"muVsTauQ2Reco2" ,xLabel = r"Q^{2} [GeV]")
compareFamilies(["q2_reco_1","q2_reco_2","q2_gen"]  ,50,0,12,"muVsTauQ2RecoWeighted" ,xLabel = r"Q^{2} [GeV]")

#compareFamilies("q2_coll"  ,50,0,12,"muVsTauQ2Coll" ,xLabel = r"Q^{2} [GeV]")

#compareFamilies("bs_pt_lhcb_alt"  ,50,0,100,"muVsTauBsPt" ,xLabel = "p_{T}(B_{s}) [GeV]")
#compareFamilies("bs_boost_lhcb_alt"  ,50,0.9,1,"muVsTauBsBoost" ,xLabel = r"#beta(B_{s})")
#compareFamilies("dxy_mu_sig"  ,50,0,20,"muVsTauIPSigmu" ,xLabel = r"IP^{xy}_{sig} (#mu)")
#compareFamilies("dz_mu_sig"  ,50,0,20,"muVsTauIPzSigmu" ,xLabel = r"IP^{z}_{sig} (#mu)")
#compareFamilies("dxy_mu"  ,50,0,0.2,"muVsTauIPmu" , xLabel = r"IP^{xy} (#mu) [cm]",log = True)
#compareFamilies("dz_mu"   ,50,0,0.5,"muVsTauIPzmu" ,xLabel = r"IP^{z} (#mu) [cm]")
#compareFamilies("lxy_bs"  ,50,0,1.8,"muVsTauLxyBs" ,xLabel = r"L_{xy} (B_{s}) [cm]")
#compareFamilies("cosMuWLhcbAlt"  ,50,-1,1,"muVsTauCosMuW" ,xLabel = r"cos(#theta_{1})")


# ds vs ds star plots
#compareFamilies("bs_pt_lhcb_alt"  ,50,0,100,"dsVsStarBsPt" ,gen1 = [0,1], gen2 = [5,6], genName = ["D_{s} - signals", "D_{s}* - signals"] ,xLabel = r" p_{T}(B_{s}) [GeV]")
#compareFamilies("ds_fitted_pt"  ,50,0,50,"dsVsStarDsPt" ,gen1 = [0,1], gen2 = [5,6], genName = ["D_{s} - signals", "D_{s}* - signals"] ,xLabel = r" p_{T}(D_{s}) [GeV]")
#compareFamilies("ds_fitted_boost"  ,50,0.9,1,"dsVsStarDsBoost" ,gen1 = [0,1], gen2 = [5,6], genName = ["D_{s} - signals", "D_{s}* - signals"] ,xLabel = "#beta(D_{s})")
#compareFamilies("lxy_ds"  ,50,0,1,"dsVsStarLxyDs" ,gen1 = [0,1], gen2 = [5,6], genName = ["D_{s} - signals", "D_{s}* - signals"] ,xLabel = "L_{xy}  [cm]", log = True)
#compareFamilies("cosMuWLhcbAlt"  ,50,-1,1,"dsVsStarCosMuW" ,gen1 = [0,1], gen2 = [5,6], genName = ["D_{s} - signals", "D_{s}* - signals"] ,xLabel = "cos(#theta_{1})")
#bs momentum comparison
#multiVar(["bs_pt_coll", "bs_pt_lhcb", "bs_pt_lhcb_alt", "bs_pt_reco_1", "bs_pt_reco_2", "bs_gen_pt"], 50,0.0,150.0, "bsPtComparison", xLabel = r"p_{T}(B_{s}) [GeV]")

#multiVar(["q2_coll", "q2_lhcb", "q2_lhcb_alt", "q2_reco_1", "q2_reco_2", "q2_gen"], 50,0.0,12.0, "q2Comparison",xLabel = r"Q^{2} [GeV^{2}]")
#multiVar(["m2_miss_coll", "m2_miss_lhcb", "m2_miss_lhcb_alt","m2_miss_reco_1", "m2_miss_reco_2",  "m2_miss_gen"], 50,0.0,8.0, "m2MissComparison", xLabel = r" m^{2}_{miss} [GeV^{2}]", log = True)

#multiVar(["bs_eta_coll", "bs_eta_lhcb", "bs_eta_lhcb_alt", "bs_eta_reco_1", "bs_eta_reco_2", "bs_gen_eta"], 50,-3,3.0, "bsEtaComparison", xLabel = r"#eta(B_{s})")
#multiVar(["bs_phi_coll", "bs_phi_lhcb", "bs_phi_lhcb_alt", "bs_phi_reco_1", "bs_phi_reco_2", "bs_gen_phi"], 50,-3.2,3.2, "bsPhiComparison", xLabel = r"#phi(B_{s})")
#multiVar(["bs_boost_coll", "bs_boost_lhcb", "bs_boost_lhcb_alt", "bs_boost_reco_1", "bs_boost_reco_2", "bs_boost_gen"], 50,0.9,1, "bsBoostComparison", xLabel = r"#beta(B_{s})")
#multiVar(["bs_boost_coll_pt", "bs_boost_lhcb_pt", "bs_boost_lhcb_alt_pt", "bs_boost_reco_1_pt", "bs_boost_reco_2_pt", "bs_boost_gen_pt"], 50,0,1, "bsPtBoostComparison")
#multiVar(["bs_boost_coll_eta", "bs_boost_lhcb_eta", "bs_boost_lhcb_alt_eta", "bs_boost_reco_1_eta", "bs_boost_reco_2_eta", "bs_boost_gen_eta"], 50,-3,3, "bsEtaBoostComparison")
#multiVar(["bs_boost_coll_phi", "bs_boost_lhcb_phi", "bs_boost_lhcb_alt_phi", "bs_boost_reco_1_phi", "bs_boost_reco_2_phi", "bs_boost_gen_phi"], 50,-3.2,3.2, "bsPhiBoostComparison")

#test the fitter
#testFit(["mu_pt", "mu_refitted_pt", "mu_gen_pt"], 50,7,30, "fitTestMuPt",xLabel = "p_{T}(#mu) [GeV]",zoom = True) 
#testFit(["phiPi_pt", "ds_fitted_pt", "ds_gen_pt"], 50,0,30, "fitTestDsPt",xLabel = "p_{T}(D_{s}) [GeV]",zoom = True) 
#testFit(["kk_pt", "phi_fitted_pt", "phi_gen_pt"], 50,0,20, "fitTestPhiPt",xLabel = "p_{T}(#phi) [GeV]",zoom = True) 
#testFit(["pi_pt", "pi_refitted_pt", "pi_gen_pt"], 50,0,30, "fitTestPiPt",xLabel = "p_{T}(#phi) [GeV]",zoom = True) 
#testFit(["k1_pt", "k1_refitted_pt", "k1_gen_pt"], 50,0,30, "fitTestK1Pt",xLabel = "p_{T}(#phi) [GeV]",zoom = True) 
#testFit(["mu_pt", "mu_refitted_pt", "mu_gen_pt"], 50,7,30, "fitTestMuPt",xLabel = "p_{T}(#mu) [GeV]") 
#testFit(["phiPi_pt", "ds_fitted_pt", "ds_gen_pt"], 50,0,30, "fitTestDsPt",xLabel = "p_{T}(D_{s}) [GeV]",zoom = True) 
#testFit(["kk_pt", "phi_fitted_pt", "phi_gen_pt"], 50,0,20, "fitTestPhiPt",xLabel = "p_{T}(#phi) [GeV]",zoom = True) 
#testFit(["pi_pt", "pi_refitted_pt", "pi_gen_pt"], 50,0,30, "fitTestPiPt",xLabel = "p_{T}(#phi) [GeV]",zoom = True) 
#testFit(["k1_pt", "k1_refitted_pt", "k1_gen_pt"], 50,0,30, "fitTestK1Pt",xLabel = "p_{T}(#phi) [GeV]",zoom = True) 
#hel angles
#collinear approx spoils this distribution as it ds mu system is back to back in bs rest frame --> dont plot it
#multiVar(["cosMuWLhcb", "cosMuWLhcbAlt", "cosMuWReco1", "cosMuWReco2", "cosMuWGen"], 30,-1,1, "cosMuWComparison",xLabel = r"cos(#theta_{1})", color = colors[1:-1]) 
#multiVar(["cosMuWLhcb", "cosMuWLhcbAlt", "cosMuWReco1", "cosMuWReco2", "cosMuWGen"], 30,-1,1, "cosMuWComparison",xLabel = "cos(#mu W)", color = colors[1:-1]) 
#multiVar(["cosMuWGen", "cosMuWGenLhcb", "cosMuWGenReco1", "cosMuWGenReco2"], 30,-1,1, "cosMuWGenComparison",norm = False, xLabel = "") 

#multiVar(["cosPhiDsColl","cosPhiDsLhcb", "cosPhiDsLhcbAlt", "cosPhiDsReco1", "cosPhiDsReco2", "cosPhiDsGen"], 50,-1,1, "cosPhiDsComparison",xLabel = r"cos(#theta_{2})") 
#multiVar(["cosPiDsLhcb", "cosPiDsLhcbAlt", "cosPiDsReco1", "cosPiDsReco2", "cosPiDsGen", "cosPiDsGenLhcb"], 50,-1,1, "cosPiDsComparison",xLabel = "") 
#multiVar(["cosPlaneBsColl","cosPlaneBsLhcb", "cosPlaneBsLhcbAlt", "cosPlaneBsReco1", "cosPlaneBsReco2", "cosPlaneBsGen"], 30,-1,1, "cosPlaneBs",xLabel = r"cos(#theta_{3})")
#multiVar(["cosPlaneDsColl","cosPlaneDsLhcb", "cosPlaneDsLhcbAlt", "cosPlaneDsReco1", "cosPlaneDsReco2", "cosPlaneDsGen"], 30,-1,1, "cosPlaneDs",xLabel = "cos Plane 2")
"""
multiVar(["cosPiK1","cosPiK1Gen"],30,-1,1,"cosPiK1", xLabel = r"cos(#theta_{4})",color = [ROOT.kBlack, ROOT.kMagenta])

#compare prefit, postfit and gen momenta

fitVtx(["pv_x",    "pv_x_gen"],                     50,0.0,0.02, "pv_x","Selected Vtx", xLabel = r"(PV)_{x} [cm]", color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["pv_y",    "pv_y_gen"],                     50,0.02,0.06, "pv_y","Selected Vtx", xLabel = r"(PV)_{y} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["pv_z",    "pv_z_gen"],                     50,-5,5, "pv_z","Selected Vtx", xLabel = r"(PV)_{z} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])

fitVtx(["sv_x",    "sv_x_gen"],                     50,-1,1, "sv_x", "Fitted Vtx",xLabel = r"(SV)_{x} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["sv_y",    "sv_y_gen"],                     50,-1,1, "sv_y", "Fitted Vtx",xLabel = r"(SV)_{y} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["sv_z",    "sv_z_gen"],                     50,-5,5, "sv_z", "Fitted Vtx",xLabel = r"(SV)_{z} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])

fitVtx(["tv_x",    "tv_x_gen"],                     50,-1,1, "tv_x", "Fitted Vtx",xLabel = r"(TV)_{x} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["tv_y",    "tv_y_gen"],                     50,-1,1, "tv_y", "Fitted Vtx",xLabel = r"(TV)_{y} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["tv_z",    "tv_z_gen"],                     50,-5,5, "tv_z", "Fitted Vtx",xLabel = r"(TV)_{z} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])

fitVtx(["fv_x",    "fv_x_gen"],                     50,-1,1, "fv_x", "Fitted Vtx",xLabel = r"(FV)_{x} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["fv_y",    "fv_y_gen"],                     50,-1,1, "fv_y", "Fitted Vtx",xLabel = r"(FV)_{y} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
fitVtx(["fv_z",    "fv_z_gen"],                     50,-5,5, "fv_z", "Fitted Vtx",xLabel = r"(FV)_{z} [cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])

#multiVar(["k1_pt",    "k1_refitted_pt",    "k1_gen_pt"],                     30,0,30, "pre_vs_postfit_k1")
#multiVar(["k2_pt",    "k2_refitted_pt",    "k2_gen_pt"],                     30,0,30, "pre_vs_postfit_k2")
#multiVar(["pi_pt",    "pi_refitted_pt",    "pi_gen_pt"],                     50,0,30, "pre_vs_postfit_pi")
#multiVar(["mu_pt",    "mu_refitted_pt",    "mu_gen_pt"],                       50,0,30, "pre_vs_postfit_mu_pt", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta])
#multiVar(["mu_eta",    "mu_refitted_eta",    "mu_gen_eta"],                    50,-3,3, "pre_vs_postfit_mu_eta", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta])
#multiVar(["mu_phi",    "mu_refitted_phi",    "mu_gen_phi"],                    50,-8,8, "pre_vs_postfit_mu_phi", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta])

#multiVar(["kk_pt",    "phi_refitted_pt",   "phi_gen_pt", "phi_fitted_pt"],     50,0,30, "pre_vs_postfit_phi_pt", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6 ])
#multiVar(["kk_eta",    "phi_refitted_eta",   "phi_gen_eta", "phi_fitted_eta"], 50,-3,3, "pre_vs_postfit_phi_eta", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])
#multiVar(["kk_phi",    "phi_refitted_phi",   "phi_gen_phi", "phi_fitted_phi"], 50,-8,8, "pre_vs_postfit_phi_phi", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])

#multiVar(["phiPi_pt", "ds_refitted_pt",    "ds_gen_pt", "ds_fitted_pt"],       50,0,30, "pre_vs_postfit_ds_pt", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])
#multiVar(["phiPi_eta", "ds_refitted_eta",    "ds_gen_eta", "ds_fitted_eta"],   50,-3,3, "pre_vs_postfit_ds_eta", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])
#multiVar(["phiPi_phi", "ds_refitted_phi",    "ds_gen_phi", "ds_fitted_phi"],   50,-8,8, "pre_vs_postfit_ds_phi", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])
"""

