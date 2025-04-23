import os,sys
import ROOT
import argparse
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import uncertainties as unc
import pdb
from scipy.stats import norm    # for gauss function
from matplotlib.patches import Patch

import mplhep as hep
hep.style.use(hep.style.CMS)

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from python.corrections import apply_muon_prefire, apply_ECAL_prefire
from common.utils import to_DateTime
from common import plotting

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetCanvasPreferGL(1)
ROOT.gStyle.SetTitleX(.3)

parser = argparse.ArgumentParser()

parser.add_argument("-r","--rates", required=True, nargs='+', help="Nominator csv file with z rates per measurement")
parser.add_argument("-n","--newRates", nargs='+', type=str, default=None, help="Nominator csv file with z rates per measurement")
parser.add_argument("-x","--xsec", type=str,
    help="csv file with z rates per measurement where xsec should be taken from (e.g. from low pileup run)")
parser.add_argument("--label",  default='Work in progress',  type=str, help="specify label ('Work in progress', 'Preliminary', )")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
parser.add_argument("--postfit", action="store_true", default=False, help="Scale by postfit values")

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
    # "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino",],
    "font.size": textsize,
    # 'text.latex.preamble': [r"""\usepackage{bm}"""]
})

mpl.rcParams.update({
    "legend.fontsize" : "medium",
    "axes.labelsize" : "medium",
    "axes.titlesize" : "medium",
    "xtick.labelsize" : "medium",
    "ytick.labelsize" : "medium",
})



# --- uncertainties on PHYSICS luminosity
include_unc_PHYSICS = run2
unc_2016 = 0.0141 # np.sqrt(1.08**2 + 0.88**2 + (0.55-0.25)**2)/100.
unc_2017 = 0.0080 # np.sqrt(0.43**2 + 0.56**2 + (0.86-0.72)**2)/100.
unc_2018 = 0.0102 # np.sqrt(0.51**2 + 0.76**2 + (0.71-0.51)**2)/100.

unc_2016_highPU = 0.012 # 0.0141
unc_2017_highPU = 0.0086 # 0.0080
unc_2018_highPU = 0.0083 # 0.0102

print("relative uncertainties attributed to PHYSICS: ")
print("2016: "+str(unc_2016))
print("2017: "+str(unc_2017))
print("2018: "+str(unc_2018))

# --- if averages should be plotted
plot_averages = run2

########## Data Acquisition ##########


# --- get Z xsec
if args.xsec:
    print("get Z cross section")
    data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])

    if data_xsec['recZCount'].dtype==object:
        data_xsec['recZCount'] = data_xsec['recZCount'].apply(lambda x: unc.ufloat_fromstr(x).n)

    apply_muon_prefire(data_xsec)
    apply_ECAL_prefire(data_xsec)

    print("apply prefire corrections - done")

    print(f"Lumi(2017H) = {data_xsec['recLumi'].sum()}")

    nz = sum(data_xsec['recZCount'])
    lumi = sum(data_xsec['recLumi'])
    if args.postfit:
        nz *= 1.00009373
        lumi *= 1.00821199

    xsec = nz/lumi
    normalize = False
else:
    # normalize everything as no cross section is specified
    xsec = 1      
    normalize = True    

    # # cross section from theory
    # xsec = 772627.88 * 0.995988004755 / 0.3646 * 0.3649 / 1000.
    # normalize = False
    
# --- z luminosity
print("get Z luminosity")
data = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.rates], ignore_index=True, sort=False)

print(f"Before: {data['recLumi'].sum()}")

if args.newRates is not None:
    data_new = pd.concat([pd.read_csv(csv, sep=',',low_memory=False) for csv in args.newRates], ignore_index=True, sort=False)

    common_rows = data.merge(data_new, on=['run', 'measurement'], suffixes=('_old', '_new'))

    # mask = data_new[['run', 'measurement']].isin(data[['run', 'measurement']].to_dict(orient='list')).all(axis=1)

    # data.loc[data['run'] > 315252, 'delLumi'] = data_new[mask]['delLumi']
    # data.loc[data['run'] > 297020, 'recLumi'] = data_new[mask]['recLumi']
    # data.loc[data['run'] > 297020, 'deadtime'] = data_new[mask]['deadtime']
    # data.loc[data['run'] > 297020, 'pileUp'] = data_new[mask]['pileUp']

    # data.loc[data.set_index(['run', 'measurement']).index.isin(common_rows.set_index(['run', 'measurement']).index), 'delLumi']

    data.loc[data.set_index(['run', 'measurement']).index.isin(common_rows.set_index(['run', 'measurement']).index), 'delLumi'] = data_new.set_index(['run', 'measurement']).reindex(data.set_index(['run', 'measurement']).index)['delLumi'].dropna().values
    data.loc[data.set_index(['run', 'measurement']).index.isin(common_rows.set_index(['run', 'measurement']).index), 'recLumi'] = data_new.set_index(['run', 'measurement']).reindex(data.set_index(['run', 'measurement']).index)['recLumi'].dropna().values
    data.loc[data.set_index(['run', 'measurement']).index.isin(common_rows.set_index(['run', 'measurement']).index), 'deadtime'] = data_new.set_index(['run', 'measurement']).reindex(data.set_index(['run', 'measurement']).index)['deadtime'].dropna().values
    data.loc[data.set_index(['run', 'measurement']).index.isin(common_rows.set_index(['run', 'measurement']).index), 'pileUp'] = data_new.set_index(['run', 'measurement']).reindex(data.set_index(['run', 'measurement']).index)['pileUp'].dropna().values

    print(f"After: {data['recLumi'].sum()}")


# --- PHYSICS luminosity
if args.postfit:
    data.loc[((data['run'] > 272006) & (data['run'] < 294645)), 'recLumi'] = data.loc[((data['run'] > 272006) & (data['run'] < 294645)), 'recLumi'] * 0.99036267
    data.loc[((data['run'] > 297045) & (data['run'] < 306463)), 'recLumi'] = data.loc[((data['run'] > 297045) & (data['run'] < 306463)), 'recLumi'] * 1.00686159
    data.loc[((data['run'] > 315251) & (data['run'] < 325176)), 'recLumi'] = data.loc[((data['run'] > 315251) & (data['run'] < 325176)), 'recLumi'] * 0.9964506

lumi_2016 = data.loc[((data['run'] > 272006) & (data['run'] < 294645)), 'recLumi'].sum()/1000.
lumi_2017 = data.loc[((data['run'] > 297045) & (data['run'] < 306463)), 'recLumi'].sum()/1000.
lumi_2018 = data.loc[((data['run'] > 315251) & (data['run'] < 325176)), 'recLumi'].sum()/1000.

print(f"Lumi(2016) = {lumi_2016}")
print(f"Lumi(2017) = {lumi_2017}")
print(f"Lumi(2018) = {lumi_2018}")

# --- Z counts

if data['recZCount'].dtype==object:
    data['recZCount'] = data['recZCount'].apply(lambda x: unc.ufloat_fromstr(x).n)

if args.postfit:
    data.loc[((data['run'] > 272006) & (data['run'] < 294645)), 'recZCount'] = data.loc[((data['run'] > 272006) & (data['run'] < 294645)), 'recZCount'] * 1.0000023
    data.loc[((data['run'] > 297045) & (data['run'] < 306463)), 'recZCount'] = data.loc[((data['run'] > 297045) & (data['run'] < 306463)), 'recZCount'] * 1.00007705
    data.loc[((data['run'] > 315251) & (data['run'] < 325176)), 'recZCount'] = data.loc[((data['run'] > 315251) & (data['run'] < 325176)), 'recZCount'] * 1.00005156

# --->>> prefire corrections
apply_muon_prefire(data)
apply_ECAL_prefire(data)
    
# <<<---

data['zLumi'] = data['recZCount'] / xsec

zlumi_2016 = data.loc[((data['run'] > 272006) & (data['run'] < 294645)), 'zLumi'].sum()/1000.
zlumi_2017 = data.loc[((data['run'] > 297045) & (data['run'] < 306463)), 'zLumi'].sum()/1000.
zlumi_2018 = data.loc[((data['run'] > 315251) & (data['run'] < 325176)), 'zLumi'].sum()/1000.

print(f"Lumi(2016) = {zlumi_2016}")
print(f"Lumi(2017) = {zlumi_2017}")
print(f"Lumi(2018) = {zlumi_2018}")



data['timeDown'] = data['beginTime'].apply(lambda x: to_DateTime(x))
data['timeUp'] = data['endTime'].apply(lambda x: to_DateTime(x))

# bring them in format to sort and plot them
data['timeDown'] = mpl.dates.date2num(data['timeDown'])
data['timeUp'] = mpl.dates.date2num(data['timeUp'])

# center of each time slice
data['time'] = data['timeDown'] + (data['timeUp'] - data['timeDown'])/2

data = data[data['recLumi'] > 0.]
data = data[data['zLumi'] > 0.]

if normalize:
    data['zLumi'] = data['zLumi'] / sum(data['zLumi']) * sum(data['recLumi'])

data['zLumi_to_dLRec'] = data['zLumi'] / data['recLumi']

invalid_runs = {
    275657, 275658, 275659, # Outliers in all those runs of 2016. HFOC was used -> problem there?
    278017, 278018          # More outliers, not clear from where
}

# quick study for invalid runs
if True:
    lumi=0
    lumi_diff=0

    for run in invalid_runs:
        data_invalid = data.loc[data['run'] == run]
        
        lumi+= data_invalid['recLumi'].sum()/1000.
        lumi_diff+= (data_invalid['recLumi'] / data['zLumi_to_dLRec'] * 0.95).sum()/1000.

    print("effected luminosity: {0}/fb".format(lumi))
    print("estimated luminosity difference: {0}/fb".format(lumi-lumi_diff))

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


def make_hist(
    df,
    run_range=None,
    zLumi_name = 'zLumi',
    refLumi_name = 'recLumi',    
    sumN=50,    # make averages of sumN measurements
    # label="Z luminosity / Ref. luminosity",
    # label="$\mathcal{L}^{Z}_{high\ PU} / \mathcal{L}_\mathrm{high\ PU}$",
    label="$(N^{Z}_{high\ PU} / N^{Z}_\mathrm{low\ PU}) / (L_{high\ PU} / L_{low\ PU})$",
    saveas="zcount",
    year=None,
    legend='upper right',
    rangey=[0.89,1.11]
):
    saveas = str(sumN) + "_" + saveas

    if run_range:
        data = df.loc[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]
        if len(df) ==0:
            return
    else:
        # skip 2017 low pileup runs
        data = df.loc[(df["run"] < 306828) | (df["run"] >= 307083)]

    data = data.sort_values(['run','time'])
    
    data['lumiratio'] = data[zLumi_name] / data[refLumi_name]

    lumi_sum = data[refLumi_name].sum()/1000.
    lumi_sum = round(lumi_sum, 1 if lumi_sum <100 else 0)

    # --- sum up each sumN rows
    lumiratio = (data[zLumi_name].groupby(data.index // sumN).sum()) / (data[refLumi_name].groupby(data.index // sumN).sum())

    # --- sum up each sumN rows
    time = data['time'].groupby(data.index // sumN).mean()
    run = data['run'].groupby(data.index // sumN).mean()
    fill = data['fill'].groupby(data.index // sumN).mean()
    # deadtime = data['deadtime'].groupby(data.index // sumN).mean()
    weight = data['weightLumi'].groupby(data.index // sumN).sum()

    # --- make histogram
    # mean and std without outliers
    mean = data['lumiratio'][(data['lumiratio']<2.0) & (data['lumiratio']>0.5)].mean()
    std = data['lumiratio'][(data['lumiratio']<2.0) & (data['lumiratio']>0.5)].std()
    # # mean and std with outliers
    # mean = data['lumiratio'].mean()
    # std = data['lumiratio'].std()

    width = 3*std
    range = (mean - width, mean + width)
    nBins = 60

    # --- make histogram
    # include overflow and underflow in last and first bin
    xx = np.array([min(max(v, mean-width), mean+width) for v in data['lumiratio'].values])
    for weighted in (False, True):
        plt.clf()
        fig = plt.figure()
        fig.subplots_adjust(left=0.15, right=0.99, top=0.93, bottom=0.125)
        ax = fig.add_subplot(111)

        if weighted:
            nEntries, bins, _ = ax.hist(xx, weights=data['weightLumi'].values, bins=nBins, range=range)
            ax.set_ylabel("Integrated luminosity [pb$^{-1}$]", fontsize=textsize)
        else:
            nEntries, bins, _ = ax.hist(xx, bins=nBins, range=range)
            ax.set_ylabel("Number of entries", fontsize=textsize)
        
        if True:
            # # plot a gaussian function with mean and std from distribution for comparison

            hist_integral = sum(nEntries * (bins[1:] - bins[:-1]))
            x = np.linspace(range[0], range[1], 100)
            plt.plot(x, hist_integral*norm.pdf(x,mean,std), color="red", linestyle="solid")

            
        ax.set_xlabel(label, fontsize=textsize)

        hep.cms.label(label=args.label, loc=0, ax=ax, data=True, year=year, lumi=lumi_sum)

        ax.text(0.97, 0.97, "$\\mu$ = {0} \n $\\sigma$ = {1}".format(round(mean,3), round(std,3)), 
            verticalalignment='top', horizontalalignment="right", transform=ax.transAxes)


        ax.set_xlim(range)

        histname = "/hist_weighted_"+saveas if weighted else "/hist_"+saveas
        print("save histogram as {0}".format(outDir+histname))
        plt.xticks(fontsize = labelsize)
        plt.yticks(fontsize = labelsize)

        plt.savefig(outDir+histname+".png")
        plt.savefig(outDir+histname+".pdf")
        plotting.write_index_and_log(outDir, histname, args=args)

        plt.close()

    # --- make scatter
    for xx, xxSum, xlabel, suffix1 in (
        # (data['time'].values, time.values, "Time", "time"),
        # (data['fill'].values, fill.values, "Fill number", "fill"),
        # (data['run'].values, run.values, "Run number", "run"),
        (data['weightLumi'].cumsum().values/1000, weight.cumsum().values/1000., "Integrated luminosity [fb$^{-1}$]", "lumi"),
    ):
        rangex = min(xx)-(max(xx)-min(xx))*0.01, max(xx)+(max(xx)-min(xx))*0.01

        for yy, yySum, ylabel, suffix in (
            (data['lumiratio'].values, lumiratio.values, label, "lumi"),
        ):
            mean = np.mean(yy)
            std = np.std(yy)

            plt.clf()
            fig = plt.figure(figsize=(10.0,4.0))
            ax = fig.add_subplot(111)
            fig.subplots_adjust(left=0.1, right=0.99, top=0.93, bottom=0.15)

            # plot uncertinty bar attributed to PHYSICS luminosity
            if suffix1 in ("lumi", ) and include_unc_PHYSICS:

                starts = np.array([rangex[0],])
                widths = np.array([abs(rangex[1] - rangex[0]),])

                if year == 2016:
                    heights = np.array([unc_2016*2,])
                    # bottoms = np.array([1. - unc_2016,])
                elif year == 2017:
                    heights = np.array([unc_2017*2,])
                    # bottoms = np.array([1. - unc_2017,])
                elif year == 2018:
                    heights = np.array([unc_2018*2,])
                    # bottoms = np.array([1. - unc_2018,])
                else:
                    starts = np.array([rangex[0], lumi_2016, lumi_2016+lumi_2017])
                    lerr = np.array([unc_2016, unc_2017, unc_2018])
                    widths = np.array([abs(rangex[0])+lumi_2016, lumi_2017, rangex[1]-(lumi_2016+lumi_2017)])
                    # bottoms = np.array([1. - unc_2016, 1. - unc_2017, 1. - unc_2018])
                    for i, y in enumerate([2016, 2017, 2018]):
                        x = starts[i]+widths[i]/2
                        ax.text(x, rangey[1]-0.005, y, va='top', ha='center')

                    heights_highPU = np.array([unc_2016_highPU, unc_2017_highPU, unc_2018_highPU])
                    # bottoms_highPU = np.array([1. - unc_2016_highPU, 1. - unc_2017_highPU, 1. - unc_2018_highPU])

                # ax.bar(starts, height=heights, width=widths, bottom=bottoms, align='edge',
                #     color="grey", alpha=0.4, hatch='/', zorder=1, #, alpha=0.6
                #     label="Ref. luminosity uncertainty")

                # ax.stairs(
                #     np.array(bottoms+heights/2.), np.append(starts, starts[-1]+widths[-1]), 
                #     baseline=None, linestyle="-", linewidth=2, color="blue", 
                #     zorder=3, 
                #     label="Ref. luminosity uncertainty"
                # )
                zval = np.array([zlumi_2016/lumi_2016, zlumi_2017/lumi_2017, zlumi_2018/lumi_2018]) 
                # zval = np.array(bottoms)

                if not args.postfit:
                    blue = "#5790fc"
                    ax.stairs(
                        zval*(1-lerr),
                        np.append(starts, starts[-1]+widths[-1]), 
                        baseline=None, linestyle="--", linewidth=2, color="#e76300", 
                        zorder=3, 
                        label="$\Delta (L_{high\ PU}/L_{low\ PU})$",
                    )
                    ax.stairs(
                        zval*(1+lerr), 
                        np.append(starts, starts[-1]+widths[-1]), 
                        baseline=None, linestyle="--", linewidth=2, color="#e76300", 
                        zorder=3, 
                    )

                    # ax.stairs(
                    #     zval-heights_highPU/2., 
                    #     np.append(starts, starts[-1]+widths[-1]), 
                    #     baseline=None, linestyle=":", linewidth=2, color="#e42536", 
                    #     zorder=3, 
                    #     label="$\Delta L_{high\ PU}$",#"Ref. luminosity uncertainty"
                    # )
                    # ax.stairs(
                    #     zval+heights_highPU/2., 
                    #     np.append(starts, starts[-1]+widths[-1]), 
                    #     baseline=None, linestyle=":", linewidth=2, color="#e42536", 
                    #     zorder=3, 
                    # )

                # ax.stairs(
                #     np.array([zlumi_2016/lumi_2016, zlumi_2017/lumi_2017, zlumi_2018/lumi_2018]), 
                #     np.array([rangex[0], lumi_2016, lumi_2016+lumi_2017, rangex[1]]), 
                #     baseline=None, linestyle="-", linewidth=2, color='#e42536', 
                #     zorder=3, label="Average")

                if args.postfit:
                    zerr = np.array([0.008148622083788956, 0.007153578874495117, 0.007134699102224342])
                else:
                    zerr = np.array([0.0072420508145137836, 0.004304985481973237, 0.004528631139759564])

#         "value": np.array([35946.39738941, 38447.87459891, 59142.73452099, 203.4959755, 133740.5024848]),
#         "hesse": np.array([292.9136076, 275.0399035, 421.96561489, 1.57895546, 942.60043282]),

                if not args.postfit:
                    purple = '#7a21dd'
                    ax.stairs(
                        zval* (1+zerr), 
                        np.array([rangex[0], lumi_2016, lumi_2016+lumi_2017, rangex[1]]), 
                        baseline=None, linestyle=":", linewidth=2, color="#ffa90e", 
                        zorder=3,
                        label="$\Delta (N^{Z}_{high\ PU}/N^{Z}_{low\ PU})$"
                        )
                    ax.stairs(
                        zval* (1-zerr), 
                        np.array([rangex[0], lumi_2016, lumi_2016+lumi_2017, rangex[1]]), 
                        baseline=None, linestyle=":", linewidth=2, color="#ffa90e", 
                        zorder=3)

                
                totalerr = np.sqrt(zerr**2 + lerr**2) if not args.postfit else zerr

                orange = '#f89c20'
                ax.stairs(
                    zval* (1+totalerr), 
                    np.array([rangex[0], lumi_2016, lumi_2016+lumi_2017, rangex[1]]), 
                    baseline=None, linestyle="-", linewidth=2, color="#bd1f01", 
                    zorder=3,
                    label="$\Delta Total$" if not args.postfit else "$\Delta L^{comb.}_{high\ PU}$"
                    )
                ax.stairs(
                    zval* (1-totalerr), 
                    np.array([rangex[0], lumi_2016, lumi_2016+lumi_2017, rangex[1]]), 
                    baseline=None, linestyle="-", linewidth=2, color="#bd1f01", 
                    zorder=3)

                ax.plot([starts,starts], (rangey[0]+0.04, rangey[1]), linestyle="--", color="black", zorder=3)            

            ax.scatter(xx, yy, s=data['weightLumi'].values, marker='.', alpha=0.6, color='black',#'', 
                label="Measurement", zorder=2)

            ww = data['weightLumi'].values,
            rms = (sum((ww*yy)**2)/sum(ww))**0.5
            print(f"RMS = {rms}")

            red = '#e42536'# 

            if suffix1 == "lumi" and plot_averages:
                # average lumi bars at centered at half of the lumi in each bar
                xxNew = np.array([xx[0]/2., ])
                xxNew = np.append(xxNew, xx[:-1]+(xx[1:] - xx[:-1])/2.)
                xx = xxNew

                xxErr = np.array([xxSum[0]/2., ])
                xxErr = np.append(xxErr, (xxSum[1:] - xxSum[:-1])/2.)
                                
                xxNew = np.array([xxSum[0]/2., ])
                xxNew = np.append(xxNew, xxSum[:-1]+(xxSum[1:] - xxSum[:-1])/2.)

                ax.errorbar(xxNew, yySum, xerr=(xxErr,xxErr), label="Average", linewidth=2, linestyle="", ecolor="#3f90da", color="#3f90da", zorder=4)

                # ax.errorbar(xxNew, yySum*1.005, xerr=(xxErr,xxErr), linestyle="dashed", ecolor='#e42536', color='#e42536', zorder=3)
                # ax.errorbar(xxNew, yySum*0.995, xerr=(xxErr,xxErr), linestyle="dashed", ecolor='#e42536', color='#e42536', zorder=3)


            hep.cms.label(label=args.label, loc=0, ax=ax, data=True, year=year, lumi=lumi_sum)

            ax.set_xlabel(xlabel, fontsize=textsize)
            ax.set_ylabel(ylabel, fontsize=textsize)
            # ax.set_ylim(rangey)
            ax.set_ylim(rangey)
            ax.set_xlim(rangex)
            if suffix1 in ("lumi", "fill", "run"):
                # plot horizontal line at 1
                ax.plot(np.array(rangex), np.array([1.,1.]), 'k--', linewidth=2)

            if suffix1 in ("fill", ):
                # plot vertical lines

                ax.plot(np.array([6166,6166]), np.array(rangey), 'b-', linewidth=1, zorder=3)
                ax.text(6168, rangey[1], "Start 8b4e scheme", rotation='vertical', verticalalignment='top',
                    fontsize=textsize*0.6, horizontalalignment='left', color="blue")
                ax.plot(np.array([6267,6267]), np.array(rangey), 'b-', linewidth=1, zorder=3)
                ax.text(6269, rangey[0], "Start leveling", rotation='vertical', verticalalignment='bottom',
                    fontsize=textsize*0.6, horizontalalignment='left', color="blue")
                # ax.plot(np.array([6070,6070]), np.array(rangey), 'r-', linewidth=1, zorder=3)
                # ax.plot(np.array([6170,6170]), np.array(rangey), 'r-', linewidth=1, zorder=3)

                def get_fill(x):
                    if len(data.loc[data['run'] > x]) * len(data.loc[data['run'] < x]) > 0:
                        return data.loc[data['run'] > x]['fill'].values[0]
                    else:
                        None 

                # trigger versions
                ax.plot(np.array([get_fill(296070),get_fill(296070)]), np.array(rangey), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(297099),get_fill(297099)]), np.array(rangey), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(297557),get_fill(297557)]), np.array(rangey), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(299368),get_fill(299368)]), np.array(rangey), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(300079),get_fill(300079)]), np.array(rangey), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(302026),get_fill(302026)]), np.array(rangey), 'r--', linewidth=1, zorder=3)
                ax.plot(np.array([get_fill(306416),get_fill(306416)]), np.array(rangey), 'r--', linewidth=1, zorder=3)

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

                    ax.plot(np.array([fill,fill]), np.array(rangey), 'k--', linewidth=1, zorder=3)
                    ax.text(fill+2, rangey[0], fill_label, verticalalignment='bottom', fontsize=textsize, horizontalalignment='left')

            if suffix1 in ("time", ):
                xfmt = mdates.DateFormatter('%Y-%m-%d')
                ax.xaxis.set_major_formatter(xfmt)
            print("save scatter as {0}".format(outDir+"/scatter_"+suffix+"_"+suffix1+"_"+saveas))
            plt.xticks(fontsize = labelsize)
            plt.yticks(fontsize = labelsize)
            if "lower" in legend:
                ncol = 3 if args.postfit else 4
            else:
                ncol = 2

            handles, labels = plt.gca().get_legend_handles_labels()
            if not args.postfit:
                p1 = Patch(color='none')
                p2 = Patch(color='none')
                p3 = Patch(color='none')

                handles = [p1, p2, p3, *handles]
                labels = ["", "", "", *labels]

            ax.legend(handles=handles, labels=labels, loc=legend, ncol=ncol, markerscale=3, scatteryoffsets=[0.5], fontsize=textsize, labelspacing=0.2, handlelength=1.0)#, frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")

            # ax.xaxis.set_label_coords(1, -0.1)
            outname = "scatter_"+suffix+"_"+suffix1+"_"+saveas

            if args.postfit:
                outname += "_postfit"

            plt.savefig(outDir+"/"+outname+".png")
            plt.savefig(outDir+"/"+outname+".pdf")
            plt.close()
            plotting.write_index_and_log(outDir, "/scatter_"+suffix+"_"+suffix1+"_"+saveas, args=args)


# make_hist(data, run_range=(297046,306462),
#     # label="$\mathcal{L}_\mathrm{Z} / \mathcal{L}_\mathrm{C}$", 
#     label="Z luminosity / Ref. luminosity", 
#     saveas="2017_zcountI", title="2017",rangey=[0.89,1.11])
 
make_hist(data, saveas="zcount", year=None, rangey=[0.91,1.05], legend="lower right")
# make_hist(data, label="ZCount(BB) / PHYSICS", saveas="zcountBB")
# make_hist(data, label="ZCount(BE) / PHYSICS", saveas="zcountBE")
# make_hist(data, label="ZCount(EE) / PHYSICS", saveas="zcountEE")
# # make_hist(data, label="ZCount(I) / PHYSICS", saveas="zcountI")

# make_hist(data, run_range=(272007,294645), saveas="2016_zcount", year=2016,rangey=[0.92,1.08])#, rangey=[0.7,1.08])
# make_hist(data, run_range=(272007,278769), saveas="2016preVFP_zcount", year=2016, rangey=[0.92,1.08])
# make_hist(data, run_range=(278769,294645), saveas="2016postVFP_zcount", year=2016, rangey=[0.92,1.08])

# make_hist(data, run_range=(297020,299329), saveas="2017B_zcount", year=2017)#,rangey=[0.85,1.15])
# make_hist(data, run_range=(299337,302029), saveas="2017C_zcount", year=2017)#,rangey=[0.85,1.15])
# make_hist(data, run_range=(302030,303434), saveas="2017D_zcount", year=2017)#,rangey=[0.85,1.15])
# make_hist(data, run_range=(303435,304826), saveas="2017E_zcount", year=2017)#,rangey=[0.85,1.15])
# make_hist(data, run_range=(304911,306462), saveas="2017F_zcount", year=2017)#,rangey=[0.85,1.15])

# make_hist(data, run_range=(297046,306462), saveas="2017_zcount", year=2017,  legend="lower right",rangey=[0.92,1.08])#)
# make_hist(data, run_range=(297046,306462), label="ZCount(I) / PHYSICS", saveas="2017_zcountI", title="2017")#,rangey=[0.85,1.15])
# make_hist(data, run_range=(297046,306462), label="ZCount(BB) / PHYSICS", saveas="2017_zcountBB", title="2017")#,rangey=[0.85,1.15])
# make_hist(data, run_range=(297046,306462), label="ZCount(BE) / PHYSICS", saveas="2017_zcountBE", title="2017")#,rangey=[0.85,1.15])
# make_hist(data, run_range=(297046,306462), label="ZCount(EE) / PHYSICS", saveas="2017_zcountEE", title="2017")#,rangey=[0.85,1.15])
# 
# make_hist(data, run_range=(315252,325175), saveas="2018_zcount", year=2018,rangey=[0.92,1.08])

# make_hist(data, saveas="zcount", year="2022")

