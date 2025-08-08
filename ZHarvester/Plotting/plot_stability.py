import os,sys
import argparse
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import uncertainties as unc
import pdb
from scipy.stats import norm    # for gauss function
import matplotlib.ticker as ticker

import mplhep as hep
hep.style.use(hep.style.ROOT)

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

sys.path.append(os.getcwd())
os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))

from common.utils import to_DateTime
from common.corrections import apply_muon_prefire, apply_ECAL_prefire

from common import parsing, plotting, logging, utils

parser = parsing.parser_plot()
parser.add_argument("-r","--rates", required=True, nargs='+', help="Nominator csv file with z rates per measurement")
args = parser.parse_args()
log = logging.setup_logger(__file__, args.verbose)

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# --- settings
colors, textsize, labelsize, markersize = plotting.set_matplotlib_style()

if args.year > 2018:
    energy=13.6
else:
    energy=13

# --- if averages should be plotted
plot_averages = True

########## Data Acquisition ##########    
# --- z luminosity
log.info("get Z luminosity")

data = utils.load_result(args.rates)

data['weight'] = data['recLumi']

# # sort out runs with low stat
# data = data.loc[data['recLumi'] > 10]

def make_hist(
    df,
    parameter,
    run_range=None,
    sumN=50,    # make averages of sumN measurements
    label="parameter",
    saveas="zcount",
    title=None,
    legend='upper right',
    rangey="auto", 
    color_hist="green", 
    color_scatter="blue"
):
    print("0")

    if df[parameter].dtype==object:
        df[parameter] = df[parameter].apply(lambda x: unc.ufloat_fromstr(x).n)
    
    lefttitle = "${"+str(energy)+"\,\mathrm{TeV}$"
    if title:
        lefttitle += " $(\mathrm{"+title+"})$"

    lefttitle=""
    if run_range:
        data = df.loc[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]
        if len(df) ==0:
            return
    else:
        # skip 2017 low pileup runs
        data = df.loc[(df["run"] < 306828) | (df["run"] >= 307083)]

    data = data.sort_values(['run','time'])

    # --- sum up each sumN rows
    time = data['time'].groupby(data.index // sumN).mean()
    run = data['run'].groupby(data.index // sumN).mean()
    fill = data['fill'].groupby(data.index // sumN).mean()
    param = data[parameter].groupby(data.index // sumN).mean()

    # deadtime = data['deadtime'].groupby(data.index // sumN).mean()
    weight = data['weight'].groupby(data.index // sumN).sum()

    # --- make histogram
    # # mean and std without outliers
    # mean = data[parameter][(data[parameter]<2.0) & (data[parameter]>0.5)].mean()
    # std = data[parameter][(data[parameter]<2.0) & (data[parameter]>0.5)].std()
    # mean and std with outliers
    mean = data[parameter].mean()
    std = data[parameter].std()

    width = 3*std
    rangex = (mean - width, mean + width)
    nBins = 60

    # plt.rcParams['font.size'] = '17.5'

    # --- make histogram
    # include overflow and underflow in last and first bin
    xx = np.array([min(max(v, mean-width), mean+width) for v in data[parameter].values])
    for weighted in (False, True):
        plt.clf()
        fig = plt.figure()
        fig.subplots_adjust(left=0.15, right=0.99, top=0.95, bottom=0.125)
        ax = fig.add_subplot(111)

        if weighted:
            nEntries, bins, _ = ax.hist(xx, weights=data['weight'].values, bins=nBins, range=rangex, color=color_hist)
            ax.set_ylabel("Integrated luminosity [pb$^{-1}$]")
        else:
            nEntries, bins, _ = ax.hist(xx, bins=nBins, range=rangex, color=color_hist)
            ax.set_ylabel("Number of entries")
        
        if True:
            # # plot a gaussian function with mean and std from distribution for comparison
            hist_integral = sum(nEntries * (bins[1:] - bins[:-1]))
            x = np.linspace(rangex[0], rangex[1], 100)
            plt.plot(x, hist_integral*norm.pdf(x,mean,std), color="red", linestyle="solid")
            
        ax.set_xlabel(label)
        hep.cms.label(
                ax=ax, 
                # fontsize=20, 
                label=args.label, 
                data=True, 
                rlabel=f"{round(data['weight'].sum()/1000,1)}"+"$\,\mathrm{fb}^{-1},$"+f"{args.year} ({energy}\, TeV)")
        ax.text(0.97, 0.97, "$\\mu$ = {0} \n $\\sigma$ = {1}".format(round(mean,3), round(std,3)), 
            verticalalignment='top', horizontalalignment="right", transform=ax.transAxes)

        ax.set_xlim(rangex)

        weighted_str = "_weighted" if weighted else ""
        outputname = f"{outDir}/hist_{parameter}{weighted_str}_{saveas}"

        log.info(f"save histogram as `{outputname}`")

        plt.savefig(f"{outputname}.png")
        plt.savefig(f"{outputname}.pdf")
        plt.close()

    plt.rcParams['font.size'] = '17.5'

    # --- make scatter
    for xx, xxSum, xlabel, suffix1 in (
        # (data['time'].values, time.values, "Time", "time"),
        (data['fill'].values, fill.values, "Fill number", "fill"),
        # (data['run'].values, run.values, "Run number", "run"),
        (data['weight'].cumsum().values/1000, weight.cumsum().values/1000., "Integrated luminosity [fb$^{-1}$]", "lumi"),
    ):
        rangex = min(xx)-(max(xx)-min(xx))*0.01, max(xx)+(max(xx)-min(xx))*0.01

        for yy, yySum, ylabel, suffix in (
            (data[parameter].values, param.values, label, "lumi"),
        ):

            plt.clf()
            fig = plt.figure(figsize=(10.0,4.0))
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.1, right=0.99, top=0.9, bottom=0.15)

            ax.scatter(xx, yy, s=data['weight'].values, marker='.', color=color_scatter, zorder=1, label="Measurement")

            # ww = data['weight'].values,
            # rms = (sum((ww*yy)**2)/sum(ww))**0.5
            # log.info(f"RMS = {rms}")

            if suffix1 in ("lumi", ) and plot_averages:
                # average lumi bars at centered at half of the lumi in each bar
                xxNew = np.array([xx[0]/2., ])
                xxNew = np.append(xxNew, xx[:-1]+(xx[1:] - xx[:-1])/2.)
                xx = xxNew

                xxErr = np.array([xxSum[0]/2., ])
                xxErr = np.append(xxErr, (xxSum[1:] - xxSum[:-1])/2.)
                                
                xxNew = np.array([xxSum[0]/2., ])
                xxNew = np.append(xxNew, xxSum[:-1]+(xxSum[1:] - xxSum[:-1])/2.)

                ax.errorbar(xxNew, yySum, xerr=(xxErr,xxErr), linestyle="", ecolor=color_hist, color=color_hist, zorder=2, label="Average")

            #hep.cms.label(ax=ax, fontsize=20, label=args.label, year=args.year, data=True)
            #hep.cms.label(ax=ax, fontsize=20, label=args.label, data=True, rlabel=f"{round(data['weight'].sum()/1000,1)}"+"$\,\mathrm{fb}^{-1},$"+f"{args.year} ({energy} TeV)")
            hep.cms.label(
                    ax=ax, 
                    #fontsize=20, 
                    label=args.label, 
                    data=True, 
                    rlabel=f"{round(data['weight'].sum()/1000,1)}"+"$\,\mathrm{fb}^{-1},$"+f"{args.year} ({energy} TeV)")

            ax.set_xlabel(xlabel)#, fontsize=20)
            ax.set_ylabel(ylabel)#, fontsize=20)

            if rangey == "auto":
                mean = np.mean(yy)
                std = np.std(yy)
                width = 4*std

                rangey_ = (mean - width, mean + width)
            else:
                rangey_ = rangey
            
            ax.set_ylim(rangey_)
            ax.set_xlim(rangex)

            # if suffix1 in ("lumi", "fill", "run", "time"):
            #     # plot horizontal line at 1
            #     ax.plot(np.array(rangex), np.array([1.,1.]), 'k--', linewidth=1, zorder=3)

            if suffix1 in ("fill", ):
                # plot vertical lines

                ax.plot(np.array([6166,6166]), np.array(rangey_), 'b-', linewidth=1, zorder=3)
                ax.text(6168, rangey_[1], "Start 8b4e scheme", rotation='vertical', verticalalignment='top',
                    fontsize=textsize*0.6, horizontalalignment='left', color="blue")
                ax.plot(np.array([6267,6267]), np.array(rangey_), 'b-', linewidth=1, zorder=3)
                ax.text(6269, rangey_[0], "Start leveling", rotation='vertical', verticalalignment='bottom',
                    fontsize=textsize*0.6, horizontalalignment='left', color="blue")
                # ax.plot(np.array([6070,6070]), np.array(rangey_), 'r-', linewidth=1, zorder=3)
                # ax.plot(np.array([6170,6170]), np.array(rangey_), 'r-', linewidth=1, zorder=3)

                def get_fill(x):
                    if len(data.loc[data['run'] > x]) * len(data.loc[data['run'] < x]) > 0:
                        return data.loc[data['run'] > x]['fill'].values[0]
                    else:
                        None 

                # trigger versions
                ax.plot(np.array([get_fill(296070),get_fill(296070)]), np.array(rangey_), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(297099),get_fill(297099)]), np.array(rangey_), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(297557),get_fill(297557)]), np.array(rangey_), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(299368),get_fill(299368)]), np.array(rangey_), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(300079),get_fill(300079)]), np.array(rangey_), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(302026),get_fill(302026)]), np.array(rangey_), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(306416),get_fill(306416)]), np.array(rangey_), 'r--', linewidth=1, zorder=3)

                # different eras
                for fill, fill_label in (
                    # 2017
                    (5838, "B"),
                    (5961, "C"),
                    (6146, "D"),
                    (6238, "E"),
                    (6297, "F"),
                    #2022
                    (get_fill(355100), "B"),
                    (get_fill(355862), "C"),
                    (get_fill(356426), "C (Muon)"),
                    (get_fill(357538), "D"),
                    (get_fill(359022), "E"),
                    (get_fill(360390), "F"),
                ):  
                    if fill==None:
                        continue

                    ax.plot(np.array([fill,fill]), np.array(rangey_), 'k--', linewidth=1, zorder=3)
                    ax.text(fill+2, rangey_[0], fill_label, verticalalignment='bottom', fontsize=textsize, horizontalalignment='left')

            if suffix1 in ("time", ):
                xfmt = mdates.DateFormatter('%Y-%m-%d')
                ax.xaxis.set_major_formatter(xfmt)

            #plt.xticks(fontsize = labelsize)
            #plt.yticks(fontsize = labelsize)

            if "lower" in legend:
                ncol = 3
            else:
                ncol = 2

            ax.legend(
                    loc=legend, 
                    #fontsize=20,
                    ncol=ncol, 
                    markerscale=3, 
                    scatteryoffsets=[0.5])
            ax.xaxis.set_label_coords(0.5, -0.1)

            outputname = f"{outDir}/scatter_{parameter}_{suffix}_{suffix1}_{saveas}"
            log.info(f"save scatter as `{outputname}`")

            plt.savefig(f"{outputname}.png")
            plt.savefig(f"{outputname}.pdf")
            plt.close()


# plot efficiencies
#for p in ["HLT", ]:#"ID", "Glo", "Sta"]:
for p in ["HLT","ID","Glo","Sta"]:

    param = "eff"+p

    label = "$\epsilon$ ("+p+")"

    color_hist="orange"
    color_scatter="red"

    if p=="HLT":
        rangey=(0.846,0.934)
    else:
        rangey="auto"

    make_hist(data, param, rangey=rangey, label=label, saveas="zcount", title=f"{args.year}", legend="lower right", color_hist=color_hist, color_scatter=color_scatter)
    # make_hist(data, param, label=label, saveas="2022BCD_zcount", title="2022(B,C,D)",  legend="lower right", run_range=(355065, 359021), color_hist=color_hist, color_scatter=color_scatter)

exit()

# plot signal model parameters
for c in ["HLT2", "HLT1", "IDFail", "GloFail", "GloPass", "StaFail", "StaPass"]:
    for p in ["mean", "sigma"]:

        log.info(f"Start param = {p+c}")
        param = p+c

        if param.startswith("mean"):
            label = "$\mu$("+param.replace("mean","")+")    "
        elif param.startswith("sigma"):
            label = "$\sigma$("+param.replace("sigma","")+")    "

        make_hist(data, param, label=label, saveas="zcount", title=f"{args.year}",  legend="lower right")
        # make_hist(data, param, label=label, saveas="2022BCD_zcount", title="2022(B,C,D)",  legend="lower right", run_range=(355065, 359021))
        # make_hist(data, param, label=label, saveas="2022EFG_zcount", title="2022(E,F,G)",  legend="lower right", run_range=(359022, 362760))

# plot background fractions
for p in ["fracGloPass","fracGloFail","fracHLT1","fracHLT2","fracIDFail","fracStaPass","fracStaFail"]:
    log.info(f"Start param {p}")

    # compute background fraction from signal and background yields
    k_sig = p.replace("frac","Nsig")
    k_bkg = p.replace("frac","Nbkg")

    if data[k_bkg].dtype==object:
        data[k_bkg] = data[k_bkg].apply(lambda x: unc.ufloat_fromstr(x).n)
    if data[k_sig].dtype==object:
        data[k_sig] = data[k_sig].apply(lambda x: unc.ufloat_fromstr(x).n)

    data[p] = data[k_bkg] / (data[k_bkg]+data[k_sig])
    data[p] = data[p].replace(np.nan,0)

    label = p
    make_hist(data, p, label=label, saveas="zcount", title=f"{args.year}", legend="lower right")	
