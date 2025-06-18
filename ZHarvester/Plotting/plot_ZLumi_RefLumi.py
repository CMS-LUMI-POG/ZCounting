import argparse
import pandas as pd
import uncertainties as unc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import os, sys
import pdb

import math

import mplhep as hep
hep.style.use(hep.style.ROOT)

sys.path.append(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from common.corrections import apply_muon_prefire, apply_ECAL_prefire
from common.utils import to_DateTime, load_input_csv

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

from common import parsing, plotting, logging

from ZUtils.python.utils import linear, pol2
from scipy.optimize import curve_fit

#mpl.rcParams["text.usetex"] = False

parser = parsing.parser_plot()
parser.add_argument("-r", "--rates", required=True, type=str, help="csv file with z rates per measurement")
parser.add_argument("-l", "--refLumi", default="", type=str, help="give a ByLs.csv as input for additional reference Luminosity")
parser.add_argument("-x", "--xsec", default="", type=str, help="csv file with z rates per measurement for absolute scale. Or number for cross section in pb to be used")
parser.add_argument("-f", "--fill", nargs="*",  type=int, default=[], help="specify a single fill to plot")
#parser.add_argument('--year', default=2024, type=int, help='Year of analysis (e.g., 2024)')
parser.add_argument("--rrange", nargs=2,  type=float, default=[0.961,1.039], help="Specify range in ratio plot")
parser.add_argument("--noFit", action="store_true", default=False, help="Don't do a fit")
parser.add_argument("--online", action="store_true", default=False, help="Run with online luminosity")
args = parser.parse_args()
log = logging.setup_logger(__file__, args.verbose)

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Settings ##########

secPerLS=float(23.3)

fmt = "png"

# plotting options
colors, textsize, labelsize, markersize = plotting.set_matplotlib_style(False)

xlabel = "LHC runtime [h]"
if args.online:
    ylabelLumi = "Inst. online luminosity [nb$^{-1}$s$^{-1}$]"
else:
    ylabelLumi = "Inst. ref. luminosity [nb$^{-1}$s$^{-1}$]"

ylabelEff = "Efficiency"

# Store the slops of the fits
slopes = {}   
slopes_pileup = {}   

# Store the lumi of each fill
yRef_lumi_per_fill = {}

 
if args.year > 2018:
    energy=13.6
else:
    energy=13

########## Data Acquisition ##########

if args.refLumi != "":
    # --- aditional external reference luminosity
    extLumi = True 

    data_ref = load_input_csv(str(args.refLumi))

    if args.fill != []:
        data_ref = data_ref.loc[data_ref['fill'].isin(args.fill)]
    else:
        log.info("Plot all fills")

    data_ref['time'] = data_ref['time'].apply(lambda x: to_DateTime(x, string_format = "mm/dd/yy"))
    data_ref['time'] = mpl.dates.date2num(data_ref['time'])

    data_ref['recorded(/nb)'] = data_ref['recorded(/pb)'] * 1000  # convert into /nb 
    data_ref['dLRec(/nb)'] = data_ref['recorded(/nb)']/secPerLS
    #Keep only what you need
    data_ref = data_ref[['fill','dLRec(/nb)','time', 'recorded(/nb)', 'avgpu']]		

else:
    extLumi = False 
  

# --- absolute scale for z luminosity
if args.xsec != "":
    if os.path.isfile(args.xsec):
        log.info("get Z cross section")
        data_xsec = pd.read_csv(args.xsec, sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
        data_xsec['zDelI'] = data_xsec['zDel'].apply(lambda x: unc.ufloat_fromstr(x).n)

        apply_muon_prefire(data_xsec)
        apply_ECAL_prefire(data_xsec)

        log.info("apply prefire corrections - done")

        xsecI = sum(data_xsec['zDelI'])/sum(data_xsec['lumiRec'])

    else:
        xsecI = int(args.xsec)
else:
    xsecI = 1


# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False) #, skiprows=[1,2,3,4,5])
if args.fill != []:
    data = data.loc[data['fill'].isin(args.fill)]

data['timeDown'] = data['beginTime'].apply(lambda x: to_DateTime(x))
data['timeUp'] = data['endTime'].apply(lambda x: to_DateTime(x))

# bring them in format to sort and plot them
data['timeDown'] = mpl.dates.date2num(data['timeDown'])
data['timeUp'] = mpl.dates.date2num(data['timeUp'])

# center of each time slice
data['time'] = data['timeDown'] + (data['timeUp'] - data['timeDown'])/2

data = data.sort_values(['fill','time'])

def unorm(x):
    # for counting experiments: define ufloat with poisson uncertainty
    return unc.ufloat(x, np.sqrt(abs(x)))

# # sort out nan values
# nan = np.isnan(data['ZRate'])
# log.info(f"Sort out {sum(nan)} nan values")
# data = data.loc[nan==False]

# # efficiency corrected number of Z bosons in timewindow, w/o deadtime corrections
# data['zDelI_mc'] = data['ZRate'] * data['timewindow'] * data['deadtime']

# data['zDelI_mc'] = data['zDelI_mc'].apply(lambda x: unorm(x))

if data['recZCount'].dtype==object:
    data['recZCount'] = data['recZCount'].apply(lambda x: unc.ufloat_fromstr(x))

data['zLumi_mc'] = data['recZCount'] / xsecI
data['zLumiInst_mc'] = data['zLumi_mc'] / data['timewindow'] * 1000  # convert into /nb
# reference lumi during one measurement
data['dLRec(/nb)'] = data['recLumi'] / data['timewindow'] * 1000  # convert into /nb 

# convert inputs to uncertainties
for key in ('effHLT', 'effID', 'effGlo', 'effSta'):
    data[key] = data[key].apply(lambda x: unc.ufloat_fromstr(x))


do_ratio=True ## activate the ratio plots 

for fill, data_fill in data.groupby("fill"):
    log.info(f"Now at fill {fill}")

    energystr = f"${energy}\,\mathrm{{TeV}} ({args.year})$"
    
    if len(data_fill) == 1:
        log.info("Only one measurement, next fill ")
        continue

    if extLumi:
        # group external reference luminosity
        # # only take range where we measured Zs
        # ref_fill = ref_fill.loc[(ref_fill['time'] >= data_fill['timeDown'].values[0]) & (ref_fill['time'] < data_fill['timeUp'].values[-1])]
        # 

        ref_fill = data_ref.loc[data_ref["fill"] == fill]

        if len(ref_fill) == 0:
            log.warning("Fill is not present in external luminosity source!")

        ref_fill_lumi = []
        for tUp, tDown in data_fill[['timeUp','timeDown']].values:
            ref_fill_lumi.append(ref_fill.loc[(ref_fill['time'] >= tDown) & (ref_fill['time'] < tUp)]['recorded(/nb)'].sum() / ((tUp - tDown) *24 * 3600) )
        
        yRefExt = np.array(ref_fill_lumi)

    yRef = data_fill['dLRec(/nb)'].values
    yRef_lumi_per_fill[fill] = yRef.sum()

    x = data_fill['time'].values
    xUp = data_fill['timeUp'].values
    xDown = data_fill['timeDown'].values

    # convert into hours
    xUp = (xUp - x) * 24 
    xDown = (x - xDown) * 24
    x = x * 24
    x = (x - x[0] + xDown[0])
    
    # x axis range
    xMin = min(x-xDown)
    xMax = max(x+xUp)
    xRange = xMax - xMin
    xMin = xMin - xRange * 0.025
    xMax = xMax + xRange * 0.025
    xRange = xMax - xMin

    if int(max(x)) > 8:
        ticksteps = 2
    else:
        ticksteps = 1
        
    xTicks = np.arange(0, int(xMax)+ticksteps, ticksteps)

    # average pileup as x axis
    xPU = data_fill['pileUp'].values
    xMinPU = min(xPU)
    xMaxPU = max(xPU)   
    xRangePU = xMaxPU - xMinPU 
    xMinPU = xMinPU - xRangePU * 0.025
    xMaxPU = xMaxPU + xRangePU * 0.025
    xRangePU = xMaxPU - xMinPU

    if int(max(xPU)) > 8:
        ticksteps = 2
    else:
        ticksteps = 1
        
    xTicksPU = np.arange(0, int(xMaxPU)+ticksteps, ticksteps)
    
    ### Make plot for efficiency vs time
    plt.clf()
    fig = plt.figure()
    
    ax1 = fig.add_subplot(111)
    
    fig.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.125, hspace=0.0)

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabelEff)
        
    maxY = 0
    minY = 1
    for eff, name, col, mkr, ms in (
        ("effGlo", "$ \epsilon_\mathrm{Glo|Sta}^\mu $", "green", "*", markersize*1.5), 
        ("effID", "$ \epsilon_\mathrm{ID|Glo}^\mu $", "red", "^", markersize*1.2), 
        ("effSta", "$ \epsilon_\mathrm{Sta|Trk}^\mu $", "blue", "s", markersize), 
        ("effHLT", "$ \epsilon_\mathrm{HLT}^\mu $", "black", "s", markersize), 
    ):
        
        y = np.array([y.n for y in data_fill[eff].values])
        yErr = np.array([y.s for y in data_fill[eff].values])
        
        # sort out zeros
        mask = y>0
        yErr = yErr[mask]
        y = y[mask]

        ax1.errorbar(x[mask], y, xerr=(xDown[mask], xUp[mask]), yerr=yErr, 
            label=name,
            marker=mkr, linewidth=0, color=col, ecolor=col, elinewidth=1.0, capsize=1.0, barsabove=True, markersize=ms,
            zorder=1)

        # ax1.plot(x, y, label=name,
        #     marker="o", linewidth=0, color=col, markersize=markersize)

        # yMax = [y.n + y.s for y in data_fill[eff].values]
        # yMin = [y.n - y.s for y in data_fill[eff].values]
                    
        maxY = max(maxY, max(y + yErr))
        minY = min(minY, min(y - yErr))

    leg = ax1.legend(loc="lower left", ncol=4)
    
    yRange = maxY - minY
    ax1.set_ylim([minY-yRange*0.25, 1.02])
    ax1.set_xlim([xMin, xMax])
    ax1.set_xticks(xTicks)

    ax1.set_xticks(xTicks)

    hep.cms.label(ax=ax1, label=args.label, data=True, rlabel=f"Fill {fill} "+f"({args.year}, {energy} TeV)")

    plt.savefig(outDir+f"/fill{fill}_eff.{fmt}")
    plt.close()    

    ### Make plot for lumi vs times
    plt.clf()
    fig = plt.figure()
    if do_ratio:
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
    else:
        ax1 = fig.add_subplot(111)
        
    fig.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.125, hspace=0.0)

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabelLumi)
        
    if args.xsec == "":
        # normalize Z luminosity to reference luminosity    
        # data_fill['zLumiInst_mc'] *= ref_fill['recorded(/pb)'].sum() / data_fill['zLumi_mc'].sum()
        # data_fill['zLumi_mc'] *= ref_fill['recorded(/pb)'].sum() / data_fill['zLumi_mc'].sum()

        data_fill['zLumiInst_mc'] *= data_fill['recLumi'].sum() / data_fill['zLumi_mc'].sum()
        data_fill['zLumi_mc'] *= data_fill['recLumi'].sum() / data_fill['zLumi_mc'].sum()


    y = np.array([yy.n for yy in data_fill['zLumiInst_mc'].values])
    yErr = np.array([y.s for y in data_fill['zLumiInst_mc'].values])

    if args.online:
        lumi_label = "Online luminosity"
    else:
        lumi_label = "Ref. luminosity"
    
    ax1.errorbar(x, yRef, xerr=(xDown, xUp), label=lumi_label, color="red", 
        linestyle='', zorder=0)

    if extLumi:
        ax1.errorbar(x, yRefExt, xerr=(xDown, xUp), label="RAMSES", color="red", linestyle='',
            zorder=0)

    ax1.errorbar(x, y, xerr=(xDown, xUp), yerr=yErr, 
        label="Z boson rate",# measurement",
        fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
        zorder=1)

    leg = ax1.legend(loc="lower left", ncol=2)

    yMin = min(min(y),min(yRef))
    yMax = max(max(y),max(yRef))

    yRange = yMax - yMin 
    ax1.set_ylim([yMin - yRange*0.45, yMax + yRange*0.15])
    ax1.set_xlim([xMin, xMax])
    
    ax1.set_xticks(xTicks)

    if do_ratio:
        ax1.set(xticklabels=[])

        ratios = []
        for zlumi, lumi in data_fill[['zLumi_mc','recLumi']].values:
            ratios.append(zlumi/lumi)

        yRatio = np.array([yy.n for yy in ratios])
        yRatioErr = np.array([y.s for y in ratios])
    
        ax2.set_xlabel(xlabel)
        ax2.set_ylabel("Z Rate / Luminosity")#"Ratio")
        
        ax2.errorbar(x, yRatio, xerr=(xDown, xUp), yerr=yRatioErr, 
            label="Z rate measurement",
            fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)
            
        ax2.plot(np.array([xMin, xMax]), np.array([1.0, 1.0]), color="black",linestyle="-", linewidth=1)
       
        # if not args.noFit:
        #     log.info("Make fit")
            
        #     func = linear
        #     popt, pcov = curve_fit(func, x, yRatio, sigma=yRatioErr, absolute_sigma=True)
        #     perr = np.sqrt(np.diag(pcov))
        #     params = unc.correlated_values(popt, pcov)
        #     log.info(params)     

        #     f = lambda x: func(x, *params)
        #     xMC = np.arange(0,100,0.5)
        #     yMC = np.array([f(x).n for x in xMC])
        #     yErrMC = np.array([f(x).s for x in xMC])

        #     nround = 3
        #     slopes[fill] = params[0].n 
            
        #     # ax2.text(0.01, 0.97, f"$f(x) = ({round(params[0].n,nround)} \\pm {round(params[0].s,nround)}) x + {round(params[1].n,nround)} \\pm {round(params[1].s,nround)}$",  verticalalignment='top', transform=ax2.transAxes,style='italic',fontsize=10)    
        #     ax2.plot(xMC, yMC, color="lime", linestyle="dashed", label="Fit")

        ax2.set_ylim(args.rrange)
        ax2.set_xlim([xMin, xMax])
        ax2.set_xticks(xTicks)
        # leg_lower = ax2.legend(loc="lower right", ncol=2,fontsize=10)



    # align y labels
    ax1.yaxis.set_label_coords(-0.12, 1.0)
    ax2.yaxis.set_label_coords(-0.12, 1.0)

    lumi = round(yRef_lumi_per_fill[fill], 1)
    if lumi>100:
        lumi = int(lumi) 
    # ax1.text(0.96, 0.94, f"$\mathcal{{L}}^\mathrm{'online' if args.online else 'ref'}_\mathrm{{int.}}={lumi}"+"\,\mathrm{pb}^{-1}$", 
    #     transform=ax1.transAxes, verticalalignment="top", horizontalalignment="right")
    hep.cms.label(ax=ax1, label=args.label, data=True, rlabel=f"Fill {fill} "+f"({args.year}, {energy} TeV)")

    plt.savefig(outDir+f"/fill{fill}_lumi.{fmt}")
    plt.close()


    ### Make plot for lumi vs pileup
    plt.clf()
    fig = plt.figure()
    if do_ratio:
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
    else:
        ax1 = fig.add_subplot(111)
        
    fig.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.125, hspace=0.0)

    ax1.set_xlabel("average pileup")
    ax1.set_ylabel(ylabelLumi)

    y = np.array([yy.n for yy in data_fill['zLumiInst_mc'].values])
    yErr = np.array([y.s for y in data_fill['zLumiInst_mc'].values])

    yRef = data_fill['dLRec(/nb)'].values   
    
    ax1.plot(xPU, yRef, #xerr=(xDown, xUp), 
        label="Reference luminosity", color="red", 
        marker = "_",
        linestyle='', zorder=0)

    ax1.errorbar(xPU, y, #xerr=(xDown, xUp), 
        yerr=yErr, 
        label="Z boson rate",
        fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
        zorder=1)

    leg = ax1.legend(loc="lower left", ncol=2)

    yMin = min(min(y),min(yRef))
    yMax = max(max(y),max(yRef))

    yRange = yMax - yMin 
    ax1.set_ylim([yMin - yRange*0.5, yMax + yRange*0.2])
    ax1.set_xlim([xMinPU, xMaxPU])

    if do_ratio:
        ax1.set(xticklabels=[])

        ratios = []
        for zlumi, lumi in data_fill[['zLumi_mc','recLumi']].values:
            ratios.append(zlumi/lumi)

        yRatio = np.array([yy.n for yy in ratios])
        yRatioErr = np.array([y.s for y in ratios])
    
        ax2.set_xlabel("average pileup")
        ax2.set_ylabel("Ratio")
        
        ax2.errorbar(xPU, yRatio, #xerr=(xDown, xUp), 
            yerr=yRatioErr, 
            label="Z rate measurement",
            fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)
            
        ax2.plot(np.array([xMinPU, xMaxPU]), np.array([1.0, 1.0]), color="black",linestyle="-", linewidth=1)

        if not args.noFit:
            log.info("Make fit")
            
            func = linear
            popt, pcov = curve_fit(func, xPU, yRatio, sigma=yRatioErr, absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
            params = unc.correlated_values(popt, pcov)

            if np.isfinite(pcov).all():
                log.info(params)     

                f = lambda x: func(x, *params)
                xMC = np.arange(0,100,0.5)
                yMC = np.array([f(x).n for x in xMC])
                yErrMC = np.array([f(x).s for x in xMC])

                slopes_pileup[fill] = params[0].n 

                nround = 3                
                # ax2.text(0.01, 0.97, f"$f(x) = ({round(params[0].n,nround)} \\pm {round(params[0].s,nround)}) x + {round(params[1].n,nround)} \\pm {round(params[1].s,nround)}$",  verticalalignment='top', transform=ax2.transAxes,style='italic',fontsize=10)    
                ax2.plot(xMC, yMC, color="lime", linestyle="dashed", label="Fit")
            else:
                log.info("Fit did not converge")                 

        ax2.set_ylim(args.rrange)
        ax2.set_xlim([xMinPU, xMaxPU])
        # ax2.set_xticks(xTicksPU)

    # align y labels
    ax1.yaxis.set_label_coords(-0.12, 1.0)
    ax2.yaxis.set_label_coords(-0.12, 0.5)

    lumi = round(yRef_lumi_per_fill[fill], 1)
    if lumi>100:
        lumi = int(lumi) 
    ax1.text(0.96, 0.94, f"$\mathcal{{L}}^{{{'online' if args.online else 'ref'}}}_\mathrm{{int.}}={lumi}"+"\,\mathrm{pb}^{-1}$", 
        transform=ax1.transAxes, verticalalignment="top", horizontalalignment="right")
    hep.cms.label(ax=ax1, label=args.label, data=True, rlabel=f"Fill {fill} "+f"({args.year}, {energy} TeV)")

    plt.savefig(outDir+f"/fill{fill}_lumi_pileup.{fmt}")
    plt.close()
    
    ## make plot eff vs. pileup
    plt.clf()
    fig = plt.figure()

    ax1 = fig.add_subplot(111)

    fig.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.125, hspace=0.0)
    
    ax1.set_xlabel("average pileup")
    ax1.set_ylabel(ylabelEff)
        
    maxY = 0
    minY = 1
    for eff, name, col, mkr, ms in (
        ("effGlo", "$ \epsilon_\mathrm{Glo|Sta}^\mu $", "green", "*", markersize*1.5), 
        ("effID", "$ \epsilon_\mathrm{ID|Glo}^\mu $", "red", "^", markersize*1.2), 
        ("effSta", "$ \epsilon_\mathrm{Sta|Trk}^\mu $", "blue", "s", markersize), 
        ("effHLT", "$ \epsilon_\mathrm{HLT}^\mu $", "black", "s", markersize), 
    ):
        
        y = np.array([y.n for y in data_fill[eff].values])
        yErr = np.array([y.s for y in data_fill[eff].values])

        # sort out zeros
        mask = y>0
        yErr = yErr[mask]
        y = y[mask]

        # y = data_fill[eff].values

        ax1.errorbar(xPU[mask], y, yerr=yErr, 
            label=name,
            marker=mkr, linewidth=0, color=col, ecolor=col, elinewidth=1.0, capsize=1.0, barsabove=True, markersize=ms,
            zorder=1)

        # xerr=(xDown, xUp)
        maxY = max(maxY, max(y + yErr))
        minY = min(minY, min(y - yErr))

    leg = ax1.legend(loc="lower left", ncol=4)
    
    yRange = maxY - minY
    ax1.set_ylim([minY-yRange*0.25, 1.01])
    ax1.set_xlim([xMinPU, xMaxPU])
    #ax1.set_xticks(xTicksPU)

    hep.cms.label(ax=ax1, label=args.label, data=True, rlabel=f"Fill {fill} "+f"({args.year}, {energy} TeV)")

    plt.savefig(outDir+f"/fill{fill}_eff_pileup.{fmt}")
    plt.close()

round_to = lambda x, n=2: x if x == 0 else round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))

# Bar plot: fill x slopes
if not args.noFit:
    for s, suffix in [(slopes, "time"), (slopes_pileup, "pileup")]:
        if len(s)<=0:
            continue

        val = np.array([v for v in s.values()])
        lumis = np.array([yRef_lumi_per_fill[k] for k in s.keys()])

        mean = np.mean(list(val))
        std = np.std(list(val))

        log.info(f"Sample of slopes with mean {mean} and std {std}")

        weighted_mean = np.sum(lumis * val) / np.sum(lumis)
        variance = np.sum(lumis * (val - weighted_mean)**2) / np.sum(lumis)
        weighted_std = np.sqrt(variance)

        log.info(f"Weighted sample of slopes with mean {weighted_mean} and std {weighted_std}")

        nbins = 20
        for xrange in ([-0.02, 0.02], [-0.002, 0.002]):
            # include outliers in first and last bin
            val[val < xrange[0]] = xrange[0]
            val[val > xrange[1]] = xrange[1]

            xlabel = "Slopes of fits"
            # Histo for the fit slopes (for each fill) weighted and unweighted
            for weighted, mu, sigma in ((True, weighted_mean, weighted_std), (False, mean, std)):
                ylabel = "Integrated luminosity [/pb]" if weighted else "Number of fills"

                fig = plt.figure()
                axh = fig.add_subplot(111)
                fig.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.125, hspace=0.0)
                axh.set_xlabel(xlabel)
                axh.set_ylabel(ylabel)
                if weighted:
                    axh.hist(val, weights=lumis, bins=nbins, range=xrange)
                else:
                    axh.hist(val, bins=nbins, range=xrange)

                axh.text(0.97, 0.97, 
                    f"$\\mu$ = {round_to(mu)} \n $\\sigma$ = {round_to(sigma)}",
                    verticalalignment='top', horizontalalignment="right", transform=axh.transAxes
                )

                axh.set_xlim(xrange)
                outname = "hist_slopes"
                if suffix is not None:
                    outname += f"_{suffix}"
                if weighted:
                    outname += "_weighted"

                outname += "_" + str(xrange[1]).replace(".","p")

                plt.savefig(outDir+f"/{outname}.{fmt}")
                plt.close()

        


