import os,sys
import ROOT
from array import array
import argparse
from datetime import datetime
import pandas as pd
import numpy as np
import shutil

latex = ROOT.TLatex()
latex.SetNDC()

import pdb

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import to_RootTime

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("--rates", required=True, type=str, help="csv file with z rates per measurement")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
parser.add_argument("-y","--year",  default=2017,  type=int, help="give output dir")
args = parser.parse_args()

year = args.year

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)


########## Data Acquisition ##########

# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
data = data.sort_values(['fill','tdate_begin','tdate_end'])

data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2

#sort out points with low statistics
#data = data[data['z_relstat'] < 0.05]

data = data.replace([np.inf, -np.inf], np.nan).dropna().dropna()


########## Plot ##########

### chi2 values of each category
for category in ('HLTeffB_chi2pass', 'HLTeffB_chi2fail', 'HLTeffE_chi2pass', 'HLTeffE_chi2fail',
          'SeleffB_chi2pass', 'SeleffB_chi2fail', 'SeleffE_chi2pass', 'SeleffE_chi2fail',
          'GloeffB_chi2pass', 'GloeffB_chi2fail', 'GloeffE_chi2pass', 'GloeffE_chi2fail',
          'GloToStaeffB_chi2pass', 'GloToStaeffB_chi2fail', 'GloToStaeffE_chi2pass', 'GloToStaeffE_chi2fail',
          'zYieldBB_chi2', 'zYieldBE_chi2', 'zYieldEE_chi2'
         ):

    graph_chi2 = ROOT.TGraph(len(data),data['time'].values,data[category].values)
    graph_chi2.SetName("graph_chi2")
    graph_chi2.SetMarkerStyle(23)
    graph_chi2.SetMarkerColor(ROOT.kAzure-4)
    graph_chi2.SetMarkerSize(1.5)
    graph_chi2.SetTitle(category)

    graph_chi2.GetYaxis().SetTitle("#chi^{2}/ndf")
    graph_chi2.GetXaxis().SetTitle("Time")
    graph_chi2.GetXaxis().SetTimeDisplay(1)
    graph_chi2.GetXaxis().SetTimeOffset(0,"gmt")
    graph_chi2.GetXaxis().SetTitleSize(0.06)
    graph_chi2.GetYaxis().SetTitleSize(0.06)
    graph_chi2.GetXaxis().SetTitleOffset(0.72)
    graph_chi2.GetYaxis().SetTitleOffset(1.1)
    graph_chi2.GetXaxis().SetLabelSize(0.05)
    graph_chi2.GetYaxis().SetLabelSize(0.05)
    #graph_chi2.GetYaxis().SetRangeUser(-0.01,0.01)
    c3=ROOT.TCanvas("c3_"+category,"c3 "+category,1000,600)
    c3.SetGrid()

    # mean, where outlier with sigma > 1 are rejected
    avg_chi2 = np.mean(data[category][abs(data[category] - np.mean(data[category])) < np.std(data[category])])

    graph_chi2.Draw("AP")

    legend=ROOT.TLegend(0.2,0.8,0.4,0.9)
    legend.AddEntry(graph_chi2,"Measurement","p")
    legend.Draw("same")

    #text=ROOT.TLatex(0.3,0.83,"CMS Automatic, produced: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    #text.SetNDC()
    #text.Draw()
    text2=ROOT.TLatex(0.2,0.23,"avg chi2: "+str(avg_chi2))
    text2.SetNDC()
    text2.Draw()

    c3.SaveAs(outDir+"/"+category+".png")
    c3.Close()

    ndf = 8
    nBins = 20
    xmin = 0.
    xmax = 2.0*avg_chi2

    th1_chi2 = ROOT.TH1F("h1_"+category, "h1 title "+category, nBins, xmin, xmax)
    th1_chi2.GetXaxis().SetTitle("#chi^{2}/ndf")
    th1_chi2.GetYaxis().SetTitle(" ")
    th1_chi2.SetTitle(" ")
    th1_chi2.GetYaxis().SetTitleFont(42)
    th1_chi2.GetXaxis().SetTitleFont(42)
    th1_chi2.SetLineColor(1)
    th1_chi2.SetLineWidth(2)
    th1_chi2.SetFillColor(0)

    for v in data[category].values:
        th1_chi2.Fill(v)

    #add bin content for underflow/overflow bins to first/last bin
    nxbins = th1_chi2.GetNbinsX()
    th1_chi2.SetBinContent(1, th1_chi2.GetBinContent(1) + th1_chi2.GetBinContent(0))
    th1_chi2.SetBinContent(nxbins, th1_chi2.GetBinContent(nxbins+1) + th1_chi2.GetBinContent(nxbins))

    # normalize to area of one
    th1_chi2.Scale(1./th1_chi2.Integral())

    f_chi2 = ROOT.TF1("fchi2","ROOT::Math::chisquared_pdf(x*{0},{0})".format(ndf), xmin, xmax)
    f_chi2.SetLineWidth(2)
    f_chi2.SetLineColor(2)

    legend=ROOT.TLegend(0.8,0.8,0.9,0.9)
    legend.SetTextSize(0.04)
    legend.SetTextAlign(12)
    legend.SetTextFont(42)
    legend.AddEntry(f_chi2,"chi2(x,{0})".format(ndf),"l")
    legend.AddEntry(th1_chi2,"fits","f")

    canvas=ROOT.TCanvas("canvas_"+category,"canvas_"+category,1000,600)
    canvas.SetLeftMargin(0.1)
    canvas.SetRightMargin(0.03)
    canvas.SetTopMargin(0.055)
    canvas.SetBottomMargin(0.1)

    textsize = 24./(canvas.GetWh()*canvas.GetAbsHNDC())

    th1_chi2.GetXaxis().SetTitleSize(textsize*1.2)
    th1_chi2.GetYaxis().SetTitleSize(textsize*1.2)
    th1_chi2.GetXaxis().SetTitleOffset(1.)
    th1_chi2.GetYaxis().SetTitleOffset(0.75)
    th1_chi2.GetXaxis().SetLabelSize(textsize*1.2)
    th1_chi2.GetYaxis().SetLabelSize(textsize*1.2)

    th1_chi2.Draw('HIST')
    f_chi2.Draw("same")
    legend.Draw("same")

    latex.SetTextFont(42)
    latex.SetTextSize(textsize*1.2)
    latex.SetTextAlign(31)

    latex.DrawLatex(0.97, 0.95, str(year))

    latex.SetTextAlign(11)
    latex.DrawLatex(0.1, 0.95, category.replace("_"," ").replace("chi2"," "))
    latex.DrawLatex(0.20, 0.87, "Preliminary")

    latex.SetTextFont(62)
    latex.DrawLatex(0.13, 0.87, 'CMS')

    canvas.SaveAs(outDir+"/hist_"+category+".png")
    canvas.Close()
