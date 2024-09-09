import pickle
import numpy as np
import ROOT
import argparse
import glob
import re
import tensorflow as tf
import matplotlib.pyplot as plt
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


parser = argparse.ArgumentParser()
parser.add_argument('model')
args = parser.parse_args()

model = args.model

# import training and validation sets

with open(f'./outputs/{model}/xx_val.pck', 'rb') as f:
    xx_val = pickle.load(f)

with open(f'./outputs/{model}/xx_train.pck', 'rb') as f:
    xx_train = pickle.load(f)

with open(f'./outputs/{model}/y_val.pck', 'rb') as f:
    y_val = pickle.load(f)

with open(f'./outputs/{model}/y_train.pck', 'rb') as f:
    y_train = pickle.load(f)


# get models 

all_models = glob.glob(f"./outputs/{model}/" + "/fold*.h5")

#find number of folds
pattern = r'fold_(\d+)'
strings = [ re.search(pattern,model).group(1) for model in all_models     ] 
ints    = [ int(re.search(pattern,model).group(1)) for model in all_models ]
argmax  = ints.index(max(ints))
nfolds  = strings[argmax]
nfolds = int(nfolds) + 1 #f.e.: if the strings go from 0 to 4 -> means we have 5 folds

#for each fold, find best model (the newest)
models = {}

for n in range(nfolds):
  all_models = glob.glob(f"./outputs/{model}/"+ f"/fold_{n}*.h5")
  pattern = rf'fold_{n}_saved-model-(\d+)'
  #take float here to avoid leading 0 problem
  strings   = [ re.search(pattern,model).group(1) for model in all_models]
  floats    = [ float(re.search(pattern,model).group(1)) for model in all_models]
  argmax    = floats.index(max(floats))
  nmodel    = strings[argmax]

  full_path = glob.glob(f"./outputs/{model}/" + f'/fold_{n}_saved-model-{nmodel}*.h5')[0]
  pattern   = r'saved-model-(.+?)\.h5' 

  model_label = re.search(pattern, full_path).group(1)
  print(f"For fold {n} we pick model {re.search(pattern, full_path).group(1)}")

  model_filename = f'./outputs/{model}/fold_{n}_saved-model-{model_label}.h5'
  model_tf = tf.keras.models.load_model(model_filename)

  models[n] = model_tf


def plotScoreOneVsAll(models, xx_train, y_train, xx_val, y_val, sig, selec = False):
  '''
    Plot the score of all channels class sig
  '''

  channels = [0,1,2,3,4,5]
  col ={0:ROOT.kCyan,1:ROOT.kGreen,2:ROOT.kMagenta,3:ROOT.kRed,4:ROOT.kBlue,5:ROOT.kGray}
  linestyles = {0:'solid',1:'dotted',2:'dashed',3:'dashdot',4:(0, (1, 10)),5:(0, (3, 5, 1, 5, 1, 5))}

  hist_content_val = {}
  hist_content_train = {}

  for chan in channels:

    dummy_val = []
    dummy_train = []

    for n in range(nfolds):

      #class predictions 1D
      #is of shape [1,2,5,4,3,2,4,2,1, ....] 
      true_1d_val   = np.argmax(y_val[n], axis=1)
      true_1d_train = np.argmax(y_train[n], axis=1)

      x_chan_val    = xx_val[n][ true_1d_val == chan ]
      x_chan_train  = xx_train[n][ true_1d_train == chan ]

 
      #...and predict their score! (data is already scaled!)
      score_val    = models[n].predict(x_chan_val) 
      score_train  = models[n].predict(x_chan_train) 

      if selec: 
        score_val   = score_val[score_val[:,5] < 0.3]
        score_train = score_train[score_train[:,5] < 0.3]

      y_chan_val    = score_val[:,sig] 
      y_chan_train  = score_val[:,sig] 

      # Plot data into histo
      hist_val   = plt.hist(y_chan_val, bins=np.arange(0,1.025,0.025))
      hist_train = plt.hist(y_chan_train, bins=np.arange(0,1.025,0.025))

      #append the bin content (at index 0) for every fold! 
      dummy_val.append(hist_val[0])
      dummy_train.append(hist_train[0])

    #average over the folds to get the average histo
    hist_content_val[chan]   = sum(dummy_val)/nfolds
    hist_content_train[chan] = sum(dummy_train)/nfolds


  #holds all histos excect the one which has chan = sig
  the_rest_val   = [0] * len(hist_content_val[0]) 
  the_rest_train = [0] * len(hist_content_train[0])

  fig = plt.figure()

  for chan in channels:

    if chan != sig:
      the_rest_val   = [a+b for a, b in zip(the_rest_val,   hist_content_val[chan])]
      the_rest_train = [a+b for a, b in zip(the_rest_train, hist_content_train[chan])]


  #bin edges are at index 1 (the same for all, can take the last one)
  bin_edges = hist_val[1]
  bin_width = np.diff(bin_edges) 
  n_bins = len(bin_edges) - 1 


  sig_val = ROOT.TH1D("","",n_bins, bin_edges)
  sig_train = ROOT.TH1D("","",n_bins, bin_edges)
  back_val = ROOT.TH1D("","",n_bins, bin_edges)
  back_train = ROOT.TH1D("","",n_bins, bin_edges)

  #fill root histo for plotting and KS
  for i in range(n_bins):
    sig_val.SetBinContent(i,hist_content_val[sig][i])
    sig_train.SetBinContent(i,hist_content_train[sig][i])
    back_val.SetBinContent(i,the_rest_val[i])
    back_train.SetBinContent(i,the_rest_train[i])

  c1=ROOT.TCanvas()
  if sig_val.Integral()!=0: sig_val.Scale(1./sig_val.Integral())
  if sig_train.Integral()!=0: sig_train.Scale(1./sig_train.Integral())
  if back_val.Integral()!=0: back_val.Scale(1./back_val.Integral())
  if back_train.Integral()!=0: back_train.Scale(1./back_train.Integral())

  c1.Draw()

  #label axes
  sig_train.GetXaxis().SetTitle("score")
  sig_train.GetYaxis().SetTitle("a.u.")
  sig_train.GetYaxis().SetRangeUser(0, max([sig_val.GetMaximum(),sig_train.GetMaximum(),back_val.GetMaximum(),back_train.GetMaximum()])*1.6)


  #signal channel is like given
  sig_train.SetFillColor(col[sig])
  sig_train.SetLineColor(col[sig])
  sig_train.SetFillStyle(3345)
  sig_train.Draw('HIST')

  sig_val.SetLineColor(col[sig])
  sig_val.SetMarkerColor(col[sig])
  sig_val.SetMarkerStyle(8)
  sig_val.Draw('EP SAME')


  back_train.SetFillColor(ROOT.kBlack)
  back_train.SetLineColor(ROOT.kBlack)
  back_train.SetFillStyle(3345)
  back_train.Draw('HIST SAME')

  back_val.SetLineColor(ROOT.kBlack)
  back_val.SetMarkerColor(ROOT.kBlack)
  back_val.SetMarkerStyle(8)
  back_val.Draw('EP SAME')

  ks_score_sig  = sig_val.KolmogorovTest(sig_train)
  print(ks_score_sig)
  ks_score_back = back_val.KolmogorovTest(back_train)
  print(ks_score_back)
  ks_value = ROOT.TPaveText(0.25, 0.76, 0.88, 0.80, 'nbNDC')
  ks_value.AddText(f'KS score of class {sig} = {round(ks_score_sig,3)} vs KS score of the rest = {round(ks_score_back,3)}')
  ks_value.SetFillColor(0)
  ks_value.Draw('EP SAME')

  leg = ROOT.TLegend(.25,.82,.83,.88)
  leg.AddEntry(sig_train ,f'Train Class {sig}' ,'F' )
  leg.AddEntry(sig_val ,f'Test Class {sig}' ,'EP' )
  leg.AddEntry(back_train ,f'Train - the rest' ,'F' )
  leg.AddEntry(back_val ,f'Test - the rest' ,'EP' )
  leg.SetNColumns(4)
  leg.Draw("SAME")


  c1.SaveAs(f"./outputs/{model}/"  + f'/plotter_one_vs_all_KS_test_{sig}.pdf')
  c1.SaveAs(f"./outputs/{model}/"  + f'/plotter_one_vs_all_KS_test_{sig}.png')
  #print('KS score: ',ks_score, len(train_pred),len(test_pred))


plotScoreOneVsAll(models, xx_train, y_train, xx_val, y_val, 1, selec = False)

