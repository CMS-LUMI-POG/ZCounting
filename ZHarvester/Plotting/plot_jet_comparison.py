import uproot
import os,sys
import ROOT
import argparse
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import uncertainties as unc
import pdb
from scipy.stats import norm, shapiro, chisquare, anderson

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from python.utils import to_DateTime, load_input_csv
from python.corrections import apply_muon_prefire, apply_ECAL_prefire

parser = argparse.ArgumentParser()

parser.add_argument("-r","--rates", required=True, nargs='+', help="Nominator csv file with z rates per measurement")
parser.add_argument("--label",  default='Work in progress',  type=str, help="specify label ('Work in progress', 'Preliminary', )")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# --- settings
run2 = True
secPerLS=float(23.3)
labelsize = 12.5
textsize = 15

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino",],
    "font.size": textsize,
    'text.latex.preamble': [r"""\usepackage{bm}"""]
})

mpl.rcParams.update({
    "legend.fontsize" : "medium",
    "axes.labelsize" : "medium",
    "axes.titlesize" : "medium",
    "xtick.labelsize" : "medium",
    "ytick.labelsize" : "medium",
})

# --- PHYSICS luminosity

lumi_2016 = 36.33
lumi_2017 = 38.48 # unprescaled 42.04

# --- uncertainties on PHYSICS luminosity
include_unc_PHYSICS = run2
unc_2016 = np.sqrt((0.012)**2 + (0.017)**2 - 2*0.012*0.017*0.26)
unc_2017 = np.sqrt((0.023)**2 + (0.017)**2 - 2*0.023*0.017*0.76)
unc_2018 = np.sqrt((0.025)**2 + (0.017)**2 - 2*0.025*0.017*0.43)

print("relative uncertainties attributed to PHYSICS: ")
print("2016: "+str(unc_2016))
print("2017: "+str(unc_2017))
print("2018: "+str(unc_2018))

# --- if averages should be plotted
plot_averages = run2

########## Data Acquisition ##########

    
# --- z luminosity
print("get Z luminosity")
data = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.rates], ignore_index=True, sort=False)

if data['recZCount'].dtype==object:
    data['recZCount'] = data['recZCount'].apply(lambda x: unc.ufloat_fromstr(x).n)

apply_muon_prefire(data)
apply_ECAL_prefire(data)
    
data['zLumi'] = data['recZCount']

data['timeDown'] = data['beginTime'].apply(lambda x: to_DateTime(x))
data['timeUp'] = data['endTime'].apply(lambda x: to_DateTime(x))

# bring them in format to sort and plot them
data['timeDown'] = mpl.dates.date2num(data['timeDown'])
data['timeUp'] = mpl.dates.date2num(data['timeUp'])

# center of each time slice
data['time'] = data['timeDown'] + (data['timeUp'] - data['timeDown'])/2

data = data[data['recLumi'] > 0.]
data = data[data['zLumi'] > 0.]

data['zLumi'] = data['zLumi'] / sum(data['zLumi']) * sum(data['recLumi'])

data['zLumi_to_dLRec'] = data['zLumi'] / data['recLumi']

invalid_runs = {
    275657, 275658, 275659, # Outliers in all those runs of 2016. HFOC was used -> problem there?
    278017, 278018          # More outliers, not clear from where
}

print("sort out invalid runs")
for run in invalid_runs:
    data = data.loc[data['run'] != run]

data['weightLumi'] = data['recLumi']

print("analyze {0} fb^-1 of data (reference lumi)".format(data['weightLumi'].sum()/1000.))
print("analyze {0} fb^-1 of data (z lumi)".format(data['zLumi'].sum()/1000.))
print("ratio: z lumi/ ref. lumi = {0}".format(data['zLumi'].sum()/data['weightLumi'].sum()))

print("Outliers:")
data_out = data.loc[abs(data['zLumi_to_dLRec']-1) > 0.1]
print(data_out[["recLumi","run","fill", "measurement","zLumi_to_dLRec","recZCount"]])


# sort out outliers
data = data.loc[abs(data['zLumi_to_dLRec']-1) < 0.1]

data["recLumi"] = data["recLumi"] / 1000

zrate = data.groupby("run")[["recZCount","recLumi"]].agg({"recZCount":"sum", "recLumi":"sum"}).reset_index()

zrate["z_xsec"] = zrate["recZCount"] / zrate["recLumi"]

# normalize to 1
zrate["z_xsec"] = zrate["z_xsec"] / (zrate["recZCount"].sum() / zrate["recLumi"].sum())

# jet info
resource_dir = "/nfs/dust/cms/user/dwalter/ZCounting/CMSSW_12_4_9/src/ZCounting/ZHarvester/res/Run2_Luminosity_studies/Zcounting"

rundict = {
    275376: "Run2016B",       
    278820: "Run2016G",
    294645: "Run2017",          
    315252: "Run2018",  
    320394: "Run2018D",  
    355064: "Run2022",
    359021: "Run2022D",
}

# combined jetlumi info
with uproot.open(f"{resource_dir}/zcountingV2.root:counts") as jettree:
    jetdata = pd.DataFrame({key: val.array(library="pd") for key, val in jettree.items()})

data = pd.merge(jetdata, zrate, on='run', how='inner')

data["jetlum"] = data["jetlum"]/data["reclumj"]

for jet_type in ("jetlum",): # ("jetlum", "jtlum"):

    #if jet_type=="jtlum":
    df_map = data["jetcount"]>4*10e3
    rangey = [0.9,1.1]
    #else:
    #    df_map = data["zcount"]>10e4
    #    rangey = [0.5,1.5]

    df = data[df_map]

    for norm in (True, False):

        for mode in ("lumi", "run"):
            if mode == "lumi":
                xkey="recLumi"
                xx = df[xkey].cumsum().values 
                xlabel = "Luminosity [fb]"
            else:
                xkey="run"
                xx = df["run"].values
                xlabel = "Run number"

            rangex = [min(xx), max(xx)]

            plt.clf()
            fig = plt.figure(figsize=(10.0,4.0))
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.16)

            if norm:
                for jet, color in (
                    ("", "red"), ("1000", "magenta"), ("500", "red"), ("100","blue")
                    ):
                    if f"{jet_type}{jet}" in df.keys():
                        yy = df[f"{jet_type}{jet}"]/df["z_xsec"]
                        ax.scatter(xx, yy, marker='.', color=color, zorder=1, label=f"Jet {jet}")
            else:
                for jet, color in (
                    ("", "red"), ("1000", "magenta"), ("500", "red"), ("100","blue")
                    ):
                    if f"{jet_type}{jet}" in df.keys():                       
                        yy = df[f"{jet_type}{jet}"]
                        ax.scatter(xx, yy, marker='.', color=color, zorder=1, label=f"Jet {jet}")
                
                ax.scatter(xx, df["z_xsec"], marker='.', color='green', zorder=1, label="Z boson")

            ax.text(0.02, 0.96, "{\\bf{CMS}} "+"\\emph{"+args.label+"} \n", verticalalignment='top', transform=ax.transAxes)

            # indicate different eras
            for run, label in rundict.items():  
                sel = df["run"] > run

                if sum(sel) <= 0:
                    continue
                if label == "Run2016B":
                    continue

                if mode == "lumi":
                    tag = df[~sel]["recLumi"].sum()
                    tag_pos = tag+(rangex[1]-rangex[0])*0.002
                    tag_pos_label = tag+(rangex[1]-rangex[0])*0.007
                else:
                    tag = df[sel]["run"].values[0]
                    tag_pos = tag-(rangex[1]-rangex[0])*0.002
                    tag_pos_label = tag+(rangex[1]-rangex[0])*0.003

                ax.plot(np.array([tag_pos,tag_pos]), np.array(rangey), 'k--', linewidth=1, zorder=3)
                ax.text(tag_pos_label, rangey[0], label.replace("Run20",""), verticalalignment='bottom', horizontalalignment='left')

            ax.set_xlabel(xlabel)
            if norm:
                ax.set_ylabel("$\sigma/<\sigma> / \sigma_\mathrm{Z}/<\sigma_\mathrm{Z}>$")
            else:
                ax.set_ylabel("$\sigma/<\sigma>$")

            ax.set_xlim(rangex)
            ax.set_ylim(rangey)

            ax.plot(np.array(rangex), np.array([1.,1.]), 'k--', linewidth=1, zorder=3)

            ax.legend(loc="best", ncol=2, markerscale=3, scatteryoffsets=[0.5], frameon=False)

            ax.tick_params(axis='both',direction='in', which='both', top=True, right=True)
            ax.xaxis.set_major_locator(ticker.AutoLocator())
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.yaxis.set_major_locator(ticker.AutoLocator())
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

            outname = f"/scatter_{mode}_z_jets"
            if jet_type == "jetlum":
                outname += "_unprescaled"
            if norm:
                outname += "_norm"

            print("save histogram as {0}".format(f"{outDir}/{outname}"))

            plt.savefig(f"{outDir}/{outname}.png")
            plt.savefig(f"{outDir}/{outname}.pdf")
            plt.close()

        for year, df_year in df.groupby("year"):

            if year == 2015:
                yearlabel = "2016preVFP"
            elif year == 2016:
                yearlabel = "2016postVFP"
            else:
                yearlabel = str(year)

            xx = df_year[f"{jet_type}"].values

            # xx500 = df_year[f"{jet_type}500"].values 
            # if jet_type != "jtlum":
            #     xx100 = df_year[f"{jet_type}100"].values 
            #     xx1000 = df_year[f"{jet_type}1000"].values 

            xlabel = "$\sigma_\mathrm{jet} / <\sigma_\mathrm{jet}>$"

            if norm:
                xx = xx / df_year["z_xsec"].values
                # if jet_type != "jtlum":
                #     xx100 = xx100 / df_year["z_xsec"].values
                #     xx1000 = xx1000 / df_year["z_xsec"].values
                xlabel += " / $\sigma_\mathrm{Z} / <\sigma_\mathrm{Z}>$"

            nBins = 50
            hist_range = [0.5, 1.5]

            xx = np.array([min(max(x, hist_range[0]), hist_range[1]) for x in xx])
            # if jet_type != "jtlum":
            #     xx100 = np.array([min(max(x, hist_range[0]), hist_range[1]) for x in xx100])
            #     xx1000 = np.array([min(max(x, hist_range[0]), hist_range[1]) for x in xx1000])

            plt.clf()
            fig = plt.figure()
            fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.125)
            ax = fig.add_subplot(111)

            nEntries, bins, _ = ax.hist(xx, bins=nBins, range=hist_range, color="red", histtype="step")
            # if jet_type != "jtlum":
            #     nEntries100, bins, _ = ax.hist(xx100, bins=nBins, range=hist_range, color="blue", histtype="step")
            #     nEntries1000, bins, _ = ax.hist(xx1000, bins=nBins, range=hist_range, color="magenta", histtype="step")
            ax.set_ylabel("Number of entries")#, fontsize=ts)
            # mean, std = xx.mean(), xx.std()

            if False:
                # # plot a gaussian function with mean and std from distribution for comparison
                hist_integral = sum(nEntries * (bins[1:] - bins[:-1]))
                x = np.linspace(hist_range[0], hist_range[1], 100)
                plt.plot(x, hist_integral*norm.pdf(x,mean,std), color="red", linestyle="solid")
                
            ax.set_xlabel(xlabel)
            ax.text(0.03, 0.97, "{\\bf{CMS}} "+"\\emph{"+args.label+"} \n", verticalalignment='top', transform=ax.transAxes)

            ax.text(0.97, 0.91, "$\\mu$ = {0} ".format(round(xx.mean(),3)), color="red",
                verticalalignment='top', horizontalalignment="right", transform=ax.transAxes)
            # if jet_type != "jtlum":
            #     ax.text(0.97, 0.97, "$\\mu$ = {0}".format(round(xx100.mean(),3)), color="blue",
            #         verticalalignment='top', horizontalalignment="right", transform=ax.transAxes)
            #     ax.text(0.97, 0.85, "$\\mu$ = {0} ".format(round(xx1000.mean(),3)), color="magenta",
            #         verticalalignment='top', horizontalalignment="right", transform=ax.transAxes)

            ax.set_ylim(0,max(nEntries)*1.2)
            ax.set_xlim(hist_range)
            ax.tick_params(axis='both',direction='in', which='both', top=True, right=True)

            ax.yaxis.set_label_coords(-0.12, 0.5)
            ax.xaxis.set_major_locator(ticker.AutoLocator())
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.yaxis.set_major_locator(ticker.AutoLocator())
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

            histname = f"/hist_{yearlabel}"
            if jet_type == "jetlum":
                histname += "_unprescaled"
            if norm:
                histname += "_norm"

            print("save histogram as {0}".format(outDir+histname))

            plt.savefig(outDir+histname+".png")
            plt.savefig(outDir+histname+".pdf")
            plt.close()
