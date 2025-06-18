import pandas as pd
import numpy as np
import argparse
import pdb
import os, sys
import uncertainties as unc
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import chisquare

sys.path.append(os.getcwd())

from common.corrections import apply_muon_prefire, apply_ECAL_prefire

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import linear, pol2

pd.options.mode.chained_assignment = None

from common import parsing, plotting, logging

parser = parsing.parser_plot()
parser.add_argument("-r","--rates", required=True, type=str, help="csv file with z rates per Measurement")
parser.add_argument("-x","--xsec",  type=str, help="csv file where xsec should be taken from (e.g. from low pileup run)")
args = parser.parse_args()
log = logging.setup_logger(__file__, args.verbose)

year = 2023

colors, textsize, labelsize, markersize = plotting.set_matplotlib_style()

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# ------------------------------------------------------------------------------
def make_plots(df,
    yAxis,
    yLabel="sigma",
    xAxis='pileUp',
    yerrors=True,
    run_range=None,
    title="",
    year="2017",
    region='inclusive',
    plot_fills=False,
    plot_all=False,
    nround=0,
    normalized=False,
    make_fit=False,
    resource=""
):
    """
    valid xAxis: 'lumi', 'pileUp', 'measurement', 'time'
    """

    if year == "2022":
        lefttitle = "$\sqrt{s}=13.6\,\mathrm{TeV}$"
    else:
        lefttitle = "$\sqrt{s}=13\,\mathrm{TeV}$"

    if year:
        lefttitle += " $(\mathrm{"+year+"})$"

    if run_range:
        df_selected = df.loc[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]
        if len(df_selected) ==0:
            return
        data = df_selected.copy()
    else:
        data = df.copy()

    if sum(data[yAxis].isnull()) > 0:
        log.info("sort out {0} points with nan".format(sum(data[yAxis].isnull())))
        data = data.loc[~data[yAxis].isnull()]

    # x_step: intervall in which the measurements are collected in one point
    if xAxis == 'lumi':
        xAxis = 'recLumiInst'
        data[xAxis] = data['recLumi'] / data['timewindow']
        xTitle="Inst. luminosity $[\\mathrm{nb}^{-1}\\mathrm{s}^{-1}]$"
        x_step = 0.2
    elif xAxis == 'time':
        data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2
        #data['time'] -= 7200
        xTitle="time"
        x_step = 1
    elif xAxis == 'pileUp':
        xTitle="average number of pileup interactions"
        x_step = 1
    else:
        xTitle=xAxis
        x_step = 0.2

    if yLabel == "sigma":
        yLabel = "$\\sigma^{\\mathrm{fid}}_{\\mathrm{Z}}$\,[pb]"

        if normalized:
            # data[yAxis] /= np.mean(data[yAxis].values)
            # data[yAxis.replace("_mc","")] /= np.mean(data[yAxis.replace("_mc","")].values)
            yLabel = "$\\sigma^{\\mathrm{fid}}_{\\mathrm{Z}} / \\langle \\sigma^{\\mathrm{fid}}_{\\mathrm{Z}} \\rangle$ "

    elif resource != "":
        # plot an efficiency on y axis - load necessery resources

        #  load MC info and format
        if len(year.split(" ")) == 1:
            suffix = year
        elif len(year.split(" ")) == 0:
            log.info("Invalid year")
        elif year.split(" ")[0] == "2016":
            if year == "2016 pre VFP" or year.split(" ")[1] in ("B","C","D","E","F"):
                suffix="2016preVFP"
            elif year == "2016 post VFP" or year.split(" ")[1] in ("G","H"):
                suffix="2016postVFP"
        else:
            suffix = year.split(" ")[0]

        basepath = "/nfs/dust/cms/user/dwalter/CMSSW_10_2_20_UL/src/potato-zcount/plots/"
        genInfoFile = basepath+'GenInfo-V13_02-d20211012-t180358/infoMC_gen_{0}.json'.format(suffix)
        recoInfoFile = basepath+'RecoInfo-V13_02-d20211012-t180324/infoMC_reco_{0}.json'.format(suffix)

        infoReco = pd.read_json(recoInfoFile,orient="index")
        infoTrue = pd.read_json(genInfoFile,orient="index")
        infoReco['bin']= infoReco.index
        infoTrue['bin']= infoTrue.index

        infoReco["True"] = False
        infoTrue["True"] = True

        info = pd.merge(infoTrue, infoReco, how='outer')
        info = info.loc[info['bin'] != 'all']
        info[['lo','hi']] = info['bin'].str.split("To", expand=True)
        info['lo'] = info['lo'].apply(lambda x: int(x.replace("nPU","")))
        info['hi'] = info['hi'].apply(lambda x: int(x.replace("Inf","75")))
        info['center'] = info['lo'] + (info['hi'] - info['lo'])/2.

        info.sort_values('center', inplace=True)

        iReco = info.loc[info['True'] == False]
        iTrue = info.loc[info['True'] == True]

        iReco.reset_index(drop=True)
        iTrue.reset_index(drop=True)

    # data['relstat'] = data[yAxis].apply(lambda x: x.s / x.n)    
    # log.info("sort out {0} points with low statistic".format(len(data['relstat'] > 0.05)))
    # data = data.loc[data["relstat"] < 0.05]

    if plot_fills:
        log.info("make plots for each fill")
        # For each Fill
        for fill, data_fill in data.groupby('fill'):

            subDir = outDir+"/PlotsXSec_"+str(fill)
            if not os.path.isdir(subDir):
                os.mkdir(subDir)
                
            plt.clf()
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.125, hspace=0.0)

            ax1.set_xlabel(xTitle)
            ax1.set_ylabel(yLabel)
            
            ax1.errorbar(data_fill[xAxis].values, 
                data_fill[yAxis].apply(lambda x: x.n).values, 
                data_fill[yAxis].apply(lambda x: x.s).values, 
                label="Measurement",
                marker="o", linewidth=0, color=colors[0], ecolor=colors[0], elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
                zorder=1)

            leg = ax1.legend(loc="upper left", ncol=3)
            
            # ax1.set_ylim([minY-yRange*0.05, maxY + yRange*0.5])
            # ax1.set_xlim([xMin, xMax])

            plt.savefig(subDir+"/"+yAxis+"_vs_"+xAxis+".png")
            plt.close()


    if plot_all:
        log.info("make plot for all measurements")
        if xAxis == 'measurement':
            xPoints = np.arange(0,len(data),1.)
        else:
            xPoints = data[xAxis].values

        plt.clf()
        fig = plt.figure(figsize=(10.0,4.0))
        ax1 = fig.add_subplot(111)
        fig.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.15)

        ax1.set_xlabel(xTitle)
        ax1.set_ylabel(yLabel)

        ax1.errorbar(xPoints, 
            data[yAxis].apply(lambda x: x.n).values, 
            data[yAxis].apply(lambda x: x.s).values, 
            label="Measurement",
            marker="o", linewidth=0, color=colors[0], ecolor=colors[0], elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)

        leg = ax1.legend(loc="upper left", ncol=3)

        # ax1.set_ylim([minY-yRange*0.05, maxY + yRange*0.5])
        # ax1.set_xlim([xMin, xMax])

        outstring = "{0}/{1}_vs_{2}".format(outDir,yAxis,xAxis)
        if run_range:
            outstring += "_run{0}to{1}".format(*run_range)

        plt.savefig(outstring+".png")
        plt.close()


        log.info("Producing all cross sections")
        log.info("the simple average cross section is "+str(sum(data['y'].values)/len(data)))

    # make fit taking into account all points
    #log.info("make fit")
    #func = linear #pol2 #quad
    #popt, pcov = curve_fit(func, data[xAxis], data['y'], sigma=data['y_Err'], absolute_sigma=True)

    #perr = np.sqrt(np.diag(pcov))
    #params = unc.correlated_values(popt, pcov)

    #f = lambda x: func(x, *params)
    #xMC = np.arange(0,100,0.5)
    #yMC = np.array([f(x).n for x in xMC])
    #yErrMC = np.array([f(x).s for x in xMC])

    xx = np.arange(data[xAxis].min()*(1./x_step)//1/(1./x_step), data[xAxis].max()*(1./x_step)//1/(1./x_step), x_step)
    xx = np.append(xx, xx[-1]+x_step)
    xx_centers = []
    yy = []
    yy0 = []

    y0Axis = yAxis.replace("_mc","")
    
    if normalized:
        data[yAxis] = data[yAxis] / (sum(data[yAxis]) / len(data[yAxis].values))
        data[y0Axis] = data[y0Axis] / ( sum(data[y0Axis]) / len(data[y0Axis].values))

    for i in range(0,len(xx)-1):
        dyy = data.loc[(data[xAxis] < xx[i+1]) & (data[xAxis] >= xx[i])]

        if len(dyy)==0:
            continue

        u_y = dyy[yAxis].values
        u_y0 = dyy[y0Axis].values 

        # 1) simple mean
        yy_avg = u_y.mean()
        yy0_avg = u_y0.mean()

        # # 2) or weighted average
        # yy_w = np.array([1./(y.s)**2 for y in u_y])
        # yy_avg = sum(u_y * yy_w) / sum(yy_w)
        # yy0_w = np.array([1./(y.s)**2 for y in u_y0])
        # yy0_avg = sum(u_y0 * yy0_w) / sum(yy0_w)
        
        if yerrors:
            yy_avg_err = yy_avg.s
            yy_avg_mean = yy_avg.n
        else:
            yy_avg_err = 0
            yy_avg_mean = yy_avg

        # sort out data points with large uncertainty
        if yy_avg == 0 or yy_avg_err/yy_avg_mean > 0.02:
            continue

        xx_centers.append((xx[i] + (xx[i+1] - xx[i]) / 2.))
        yy.append(yy_avg)
        yy0.append(yy0_avg)

    xx = np.array(xx_centers)

    if yerrors:
        yy_err = np.array([y.s for y in yy])
        yy = np.array([y.n for y in yy])
        yy0_err = np.array([y.s for y in yy0])
        yy0 = np.array([y.n for y in yy0])

        yMin = min(yy - yy_err)
        yMax = max(yy + yy_err)

    else:
        yy = np.array([y for y in yy])
        yy0 = np.array([y for y in yy0])

        yMin = min(yy)
        yMax = max(yy)

    yRange = abs(yMax - yMin)    

    if make_fit:
        log.info("make fit")
        func = linear #pol2 #quad
        #popt, pcov = curve_fit(func, xx, yy, sigma=yy_err, absolute_sigma=True)
        popt, pcov, infodict, mesg, ier = curve_fit(func, xx, yy, sigma=yy_err, absolute_sigma=True, full_output = True)

        # Compute chi square
        r = yy - func(xx, *popt)
        chisq = np.sum((r/yy_err)**2)
        log.info(chisq)
        ndf = len(yy) - len(popt)
        chisq_ndf = chisq/ndf

        perr = np.sqrt(np.diag(pcov))
        params = unc.correlated_values(popt, pcov)

        log.info(f"Fit parameters: {params}")

        f = lambda x: func(x, *params)
        xMC = np.arange(0,100,0.5)
        yMC = np.array([f(x).n for x in xMC])
        yErrMC = np.array([f(x).s for x in xMC])

    xx_err = np.ones(len(xx_centers))*x_step/2
    
    plt.clf()
    fig = plt.figure(figsize=(10.0,4.0))
    ax1 = fig.add_subplot(111)
    fig.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.16)

    ax1.set_xlabel(xTitle)
    ax1.set_ylabel(yLabel)
    ax1.text(0.74, 0.97, "\\bf{CMS}", verticalalignment='top', transform=ax1.transAxes, weight="bold")
    ax1.text(0.81, 0.97, "\\emph{"+args.label+"}", verticalalignment='top', transform=ax1.transAxes,style='italic')    
    ax1.text(0.74, 0.89, "13.6 TeV (2022)", verticalalignment='top', transform=ax1.transAxes,style='italic') 
    if make_fit:
        ax1.text(0.01, 0.87, " $\chi^2/ndf$ = "+str(chisq_ndf), verticalalignment='top', transform=ax1.transAxes, weight="bold")

        nround = 5
        ax1.text(0.01, 0.97, 
        f"$f(x) = ({round(params[0].n,nround)} \\pm {round(params[0].s,nround)}) x + {round(params[1].n,nround)} \\pm {round(params[1].s,nround)}$", 
        verticalalignment='top', transform=ax1.transAxes,style='italic')    

    xMin = min(xx_centers) - x_step/2
    xMax = max(xx_centers) + x_step/2
    xRange = abs(xMin - xMax)

    # p4 = ax1.errorbar(xx, yy0, xerr=xx_err, yerr=yy0_err, label="Measurements uncorrected",
    #    marker="o", linewidth=0, color="grey", ecolor="grey", elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
    #    zorder=1)
    if yerrors:
        p3 = ax1.errorbar(xx, yy, xerr=xx_err, yerr=yy_err, label="Measurements",
            marker="o", linewidth=0, color="black", ecolor="black", elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)
    else:
        p3 = ax1.errorbar(xx, yy, xerr=xx_err, label="Measurements",
            marker="o", linewidth=0, color="black", ecolor="black", elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)
                    
    ax1.plot(np.array([xMin-xRange*0.02,xMax+xRange*0.02]),np.array([1.,1.]), "k--", linewidth=1)
    
    if make_fit:
        p2 = ax1.plot(xMC, yMC, color=colors[0],linestyle="solid",  label="Fit", linewidth=1)
    
        p1 = ax1.fill_between(xMC, yMC - yErrMC, yMC + yErrMC,
                        color='grey', alpha=0.2, zorder=1) 
        p1 = ax1.fill(np.NaN, np.NaN, color='grey', alpha=0.2, linewidth=0.)    

    # make a separate entry for uncorrected
    # leg_styles = [p3, p4, (p2[0], p1[0])]
    # leg_labels = ['Measurements', 'Uncorrected', 'Fit']
    
        leg_styles = [p3, (p2[0], p1[0])]
        leg_labels = ['Measurements', 'Fit']
    else:
        leg_styles = [p3, ]
        leg_labels = ['Measurements', ]        

    leg = ax1.legend(leg_styles, leg_labels, loc="lower left", ncol=len(leg_labels))

    ax1.set_ylim([yMin-yRange*0.1, yMax+yRange*0.05])
    ax1.set_xlim([xMin-xRange*0.02, xMax+xRange*0.02])

    outstring = "{0}/{1}_vs_{2}".format(outDir,yAxis,xAxis)
    outstring += year.replace(" ","_")

    log.info(f"Save linearity plot `{outstring}`")

    plt.savefig(outstring+".png")
    plt.savefig(outstring+".pdf")
    plt.close()

# ------------------------------------------------------------------------------

plot_lowPU = False
if args.xsec:
    log.info("load low pileup data")
    plot_lowPU = True
    # --- get Z xsec
    data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])

    if data_xsec['recZCount'].dtype==object:
        data_xsec['recZCount'] = data_xsec['recZCount'].apply(lambda x: unc.ufloat_fromstr(x))

    apply_muon_prefire(data_xsec)
    apply_ECAL_prefire(data_xsec)

    log.info("apply prefire corrections - done")

    xsec = sum(data_xsec['recZCount'])/sum(data_xsec['recLumi'])

    

# ------------------------------------------------------------------------------
log.info("load csv file in dataframe")

rates = pd.read_csv(args.rates, sep=',')

invalid_runs = {
    275657, 275658, 275659, # Outliers in all those runs of 2016. HFOC was used -> problem there?
    278017, 278018          # More outliers, not clear from where
}
for run in invalid_runs:
    rates = rates.loc[rates['run'] != run]

# convert to uncertainties
for key in ['recZCount',]:
    rates[key] = rates[key].apply(lambda x: unc.ufloat_fromstr(x))

rates = rates[rates['recLumi'] > 0.]
rates = rates[rates['recZCount'] > 0.]

apply_ECAL_prefire(rates, apply_on="recZCount")
apply_muon_prefire(rates, apply_on="recZCount")

# # sort out outliers
# rates['zLumi'] = rates['recZCount'] / xsec
# rates['zLumi_to_dLRec'] = rates['zLumi'] / rates['recLumi']
# rates = rates.loc[abs(rates['zLumi_to_dLRec']-1) < 0.1]

if args.xsec:
    # include the lowPU data in the dataframe
    rates = pd.concat([rates, data_xsec])

rates['recLumi'] = rates['recLumi'] * 1000    # convert into /nb

rates['xsec_mc'] = rates['recZCount'] / rates['recLumi']

rates['xsec'] = rates['xsec_mc'] / (rates['cIO']**2 * rates['cID'] * rates['cHLT'] * rates['cAcceptanceInner'] * rates['cAcceptanceOuter'])

# only keep necessary data
# rates = rates[["xsec_mc", "xsec", "recLumi", "timewindow", "run", "pileUp", "effHLT", "effID", "effGlo", "effSta"]]

for xAxis in ["pileUp","lumi" ]:
    for key, ylabel, region, mcRes in (
        # ("xsec_mc", "sigma", "I", ""),
        ('effHLT' ,'Muon HLT efficiency', "", "HLT"),
        ('effID' ,'Muon identification efficiency', "", "ID"),
        ('effGlo' ,'Muon global efficiency', "", "Glo"),
        ('effSta' ,'Muon standalone efficiency', "", "Sta"),
        ('meanHLT2' ,'$\mu$(HLT2)', "", ""),
        ('meanHLT1' ,'$\mu$(HLT1)', "", ""),
        ('meanIDFail' ,'$\mu$(IDFail)', "", ""),
        ('meanGloPass' ,'$\mu$(GloPass)', "", ""),
        ('meanGloFail' ,'$\mu$(GloFail)', "", ""),
        ('meanStaPass' ,'$\mu$(StaPass)', "", ""),
        ('meanStaFail' ,'$\mu$(StaFail)', "", ""),
        ('sigmaHLT2' ,'$\sigma$(HLT2)', "", ""),
        ('sigmaHLT1' ,'$\sigma$(HLT1)', "", ""),
        ('sigmaIDFail' ,'$\sigma$(IDFail)', "", ""),
        ('sigmaGloPass' ,'$\sigma$(GloPass)', "", ""),
        ('sigmaGloFail' ,'$\sigma$(GloFail)', "", ""),
        ('sigmaStaPass' ,'$\sigma$(StaPass)', "", ""),
        ('sigmaStaFail' ,'$\sigma$(StaFail)', "", ""),
    ):
        if rates[key].dtype==object:
            rates[key] = rates[key].apply(lambda x: unc.ufloat_fromstr(x))
        
        if key.startswith("xsec"):
            normalized=True
            make_fit=True      
        else:
            normalized=False
            make_fit=False
        
        if key.startswith("sigma") or key.startswith("mean"): 
            yerrors = False
        else:
            yerrors = True
        # # # single eras
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 B", run_range=(272007,275376), normalized=False)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 C", run_range=(275657,276283), normalized=False)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 D", run_range=(276315,276811), normalized=False)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 E", run_range=(276831,277420), normalized=False)
        # # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 F", run_range=(277772,278808), normalized=False)
        # # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 G", run_range=(278820,280385), normalized=False)
        # # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 H", run_range=(280919,284044), normalized=False)
        # # 
        # # # 2017
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, year="2017 B", title="corrected", run_range=(297046,299329), normalized=False)
        # # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, year="2017 C", title="corrected", run_range=(299368,302029), normalized=False)
        # # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, year="2017 D", title="corrected", run_range=(302030,303434), normalized=False)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, year="2017 E", title="corrected", run_range=(303434,304797), normalized=False)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, year="2017 F", title="corrected", run_range=(305040,306462), normalized=False)
        # # 
        # # # ## 2018
        # # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 A", run_range=(315252,316995), normalized=False)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 B", run_range=(317080,319310), normalized=False)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 C", run_range=(319337,320065), normalized=False)
        # # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018 D", run_range=(320673,325175), normalized=False)

        # # # total 2016
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016", run_range=(272007,294645), normalized=True, xAxis=xAxis)
        # total 2016 pre VFP
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 pre VFP", run_range=(272007,278769), normalized=True, xAxis=xAxis)
        # total 2016 post VFP
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2016 post VFP", run_range=(278769,294645), normalized=True, xAxis=xAxis)
        # # total 2017
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2017", run_range=(297046,307083), normalized=True, xAxis=xAxis)
        # # ## total 2018
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, resource=mcRes, title="corrected", year="2018", run_range=(315252,325175), normalized=True, xAxis=xAxis)

        ### Run 3
        # 2022
        make_plots(rates, yAxis=key, yLabel=ylabel, region=region, title="corrected", year="2022", normalized=normalized, make_fit=make_fit, xAxis=xAxis, yerrors=yerrors)
        make_plots(rates, yAxis=key, yLabel=ylabel, region=region, title="corrected", year="2022 BCD", run_range=(355065,359021), normalized=normalized, make_fit=make_fit, xAxis=xAxis, yerrors=yerrors)
        make_plots(rates, yAxis=key, yLabel=ylabel, region=region, title="corrected", year="2022 EFG", run_range=(359022,362760), normalized=normalized, make_fit=make_fit, xAxis=xAxis, yerrors=yerrors)
        # make_plots(rates, yAxis=key, yLabel=ylabel, region=region, title="corrected", year="2022", normalized=normalized, xAxis='pileUp')

