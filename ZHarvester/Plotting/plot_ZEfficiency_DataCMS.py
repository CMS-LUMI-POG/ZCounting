import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas
import numpy as np
import shutil
import pdb

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cms", type=str, help="give the CMS csv per Measurement as input", required=True)
parser.add_argument("-s", "--saveDir", default='./', type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Data Acquisition ##########

data = pandas.read_csv(str(args.cms), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','tdate_begin','tdate_end'])

# remove rows with invalid Z rates and low statistics
#data = data[np.isfinite(data['ZRate'])]
#data = data[data['z_relstat'] < 0.05]

data['tdate'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2

meta = dict()

########## Plot ##########
fills = data.drop_duplicates('fill')['fill'].values

features = (
    ('ZBBeff_mc','corrected Z-BB-Reconstruction efficiency', 0.7, 1.0),
    ('ZBEeff_mc','corrected Z-BE-Reconstruction efficiency', 0.7, 1.0),
    ('ZEEeff_mc','corrected Z-EE-Reconstruction efficiency', 0.7, 1.0),
    ('ZBBeff'  ,'Z-BB-Reconstruction efficiency', 0.7, 1.0),
    ('ZBEeff'  ,'Z-BE-Reconstruction efficiency', 0.7, 1.0),
    ('ZEEeff'  ,'Z-EE-Reconstruction efficiency', 0.7, 1.0),
    ('HLTeffB' ,'Muon HLT-B efficiency',0.8, 1.0 ),
    ('HLTeffE' ,'Muon HLT-E efficiency',0.8, 1.0),
    ('SeleffB' ,'Muon Sel-B efficiency',0.8, 1.0),
    ('SeleffE' ,'Muon Sel-E efficiency',0.8, 1.0),
    ('GloeffB' ,'Muon Glo to Trk -B efficiency',0.9, 1.0),
    ('GloeffE' ,'Muon Glo to Trk -E efficiency',0.9, 1.0),
    ('GloToStaeffB' ,'Muon Glo to Sta -B efficiency',0.98, 1.01),
    ('GloToStaeffE' ,'Muon Glo to Sta -E efficiency',0.98, 1.01),
    #('StaeffB' ,'Muon Sta-B efficiency',0.9, 1.0),
    #('StaeffE' ,'Muon Sta-E efficiency',0.9, 1.0),
    #('TrkeffB' ,'Muon Trk-B efficiency',0.95,1.01),
    #('TrkeffE' ,'Muon Trk-E efficiency',0.95,1.01),
    ('zYieldBB_purity'    ,'Z BB purity',0.9,1.0),
    ('zYieldBE_purity'    ,'Z BE purity',0.9,1.0),
    ('zYieldEE_purity'    ,'Z EE purity',0.9,1.0),
    )

##### loop over Fills and produce fill specific plots
for fill in fills:
    dFill = data.loc[data['fill'] == fill]

    subDir = outDir+"/PlotsFill_"+str(fill)
    if not os.path.isdir(subDir):
        os.mkdir(subDir)

    ### Efficiency ###

    for eff, name, ymin, ymax in features:
        graph_Zeff = ROOT.TGraph(len(dFill),dFill['tdate'].values,dFill[eff].values )
        graph_Zeff.SetName("graph_Zeff")
        graph_Zeff.SetMarkerStyle(22)
        graph_Zeff.SetMarkerColor(ROOT.kOrange+8)
        graph_Zeff.SetFillStyle(0)
        graph_Zeff.SetMarkerSize(1.5)
        graph_Zeff.GetXaxis().SetTimeDisplay(1)
        graph_Zeff.SetTitle(name+", Fill "+str(fill))
        graph_Zeff.GetYaxis().SetTitle("Efficiency")
        if(eff == 'Zfpr'):
            graph_Zeff.GetYaxis().SetTitle("b/(s + b)")
        graph_Zeff.GetYaxis().SetTitleSize(0.07)
        graph_Zeff.GetYaxis().SetTitleOffset(1.1)
        graph_Zeff.GetXaxis().SetTitle("Time")
        graph_Zeff.GetXaxis().SetTitleSize(0.06)
        graph_Zeff.GetXaxis().SetTitleOffset(0.75)
        graph_Zeff.GetXaxis().SetLabelSize(0.05)
        graph_Zeff.GetYaxis().SetLabelSize(0.05)
        graph_Zeff.GetYaxis().SetRangeUser(ymin,ymax)

        c1=ROOT.TCanvas("c1","c1",1000,600)
        c1.cd(1)
        graph_Zeff.Draw("AP")
        legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
        legend.AddEntry(graph_Zeff,"CMS","pe")
        #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        #text1.SetNDC()
        #text1.Draw()
        c1.SaveAs(subDir+"/"+str(eff)+"_"+str(fill)+".png")
        c1.Close()
        if eff in meta.keys():
            meta[eff].append(np.mean(dFill[eff].values))
        else:
            meta[eff] = [np.mean(dFill[eff].values)]

### Efficiency of all fills###

for eff, name, ymin, ymax in features:
    graph_meta = ROOT.TGraph(len(fills),fills.astype(float),np.array(meta[eff]))
    graph_meta.SetName("graph_meta")
    graph_meta.SetMarkerStyle(22)
    graph_meta.SetMarkerColor(ROOT.kOrange+8)
    graph_meta.SetFillStyle(0)
    graph_meta.SetMarkerSize(1.5)
    graph_meta.SetTitle(name+", Fill averages ")
    graph_meta.GetYaxis().SetTitle("Efficiency")
    if(eff == 'Zfpr'):
        graph_meta.GetYaxis().SetTitle("b/(s + b)")
    graph_meta.GetYaxis().SetTitleSize(0.07)
    graph_meta.GetYaxis().SetTitleOffset(1.1)
    graph_meta.GetXaxis().SetTitle("Fill")
    graph_meta.GetXaxis().SetTitleSize(0.06)
    graph_meta.GetXaxis().SetTitleOffset(0.75)
    graph_meta.GetXaxis().SetLabelSize(0.05)
    graph_meta.GetYaxis().SetLabelSize(0.05)
    graph_meta.GetYaxis().SetRangeUser(ymin,ymax)

    c1=ROOT.TCanvas("c1","c1",1000,600)
    c1.cd(1)
    graph_meta.Draw("AP")
    legend=ROOT.TLegend(0.65,0.65,0.9,0.9)
    legend.AddEntry(graph_meta,"CMS","pe")
    #text1=ROOT.TText(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text1.SetNDC()
    #text1.Draw()
    c1.SaveAs(outDir+"/allSummary_"+str(eff)+".png")
    c1.Close()
