import os
import ROOT
import uncertainties as unc
import numpy as np

from common.logging import child_logger
log = child_logger(__name__)


# ------------------------------------------------------------------------------
def load_input_csv(byLS_data):
    """
    load input by lumisection CSV file, convert it and return the data 
    
    Parameters
    ----------
    byLS_data : str
        path to the file by lumisection CSV file
    """

    import pandas as pd
    
    byLS_file = open(str(byLS_data))
    byLS_lines = byLS_file.readlines()
    byLS_data = pd.read_csv(byLS_data, sep=',', low_memory=False,
        skiprows=lambda x: byLS_lines[x].startswith('#') and not byLS_lines[x].startswith('#run'))
        
    print("INFO:  === formatting csv file...")    # formatting the csv
    byLS_data[['run', 'fill']] = byLS_data['#run:fill'].str.split(':', expand=True).apply(pd.to_numeric)
    byLS_data['ls'] = byLS_data['ls'].str.split(':', expand=True)[0].apply(pd.to_numeric)   

    if 'delivered(/ub)' in byLS_data.columns.tolist():  # convert to /pb
        byLS_data['delivered(/ub)'] = byLS_data['delivered(/ub)'].apply(lambda x: x / 1000000.)
        byLS_data['recorded(/ub)'] = byLS_data['recorded(/ub)'].apply(lambda x: x / 1000000.)
        byLS_data = byLS_data.rename(index=str, columns={'delivered(/ub)': 'delivered(/pb)', 'recorded(/ub)': 'recorded(/pb)'})
    elif 'delivered(/fb)' in byLS_data.columns.tolist():  # convert to /pb
        byLS_data['delivered(/fb)'] = byLS_data['delivered(/fb)'].apply(lambda x: x * 1000.)
        byLS_data['recorded(/fb)'] = byLS_data['recorded(/fb)'].apply(lambda x: x * 1000.)
        byLS_data = byLS_data.rename(index=str, columns={'delivered(/fb)': 'delivered(/pb)', 'recorded(/fb)': 'recorded(/pb)'})

    byLS_data['recorded(/pb)'] = byLS_data['recorded(/pb)'].fillna(0)

    # if there are multiple entries of the same ls (for example from different triggers), 
    #   only keep the one with the highest luminosity.
    byLS_data = byLS_data.sort_values(['fill', 'run', 'ls', 'delivered(/pb)', 'recorded(/pb)'])
    byLS_data = byLS_data.drop_duplicates(['fill', 'run', 'ls'])
    
    return byLS_data

# ------------------------------------------------------------------------------
def unorm(value):
    if value > 0:
        return unc.ufloat(value, np.sqrt(value))
    else: 
        return unc.ufloat(0.0, 0.0)

# ------------------------------------------------------------------------------
def getFileName(directory, run):
    """
    Get root file with histograms
    
    Parameters
    ----------
    directory : str
        directory where the root files are stored
    run : integer
        run number for the rootfile to load
    """
    import glob
    
    # check if run was processed already
    eosFileList = glob.glob(directory + '/*' + str(run) + '*.root')
    if not len(eosFileList) > 0:
        # look one level deeper
        eosFileList = glob.glob(directory + '/*/*' + str(run) + '*.root')
    if not len(eosFileList) > 0:
        print("WARNING: === The file does not (yet) exist for run: " + str(run))
        return None
        
    elif len(eosFileList) > 1:
        print("WARNING: === Multiple files found for run: " + str(run))
        return None
    else:
        return eosFileList[0]

# ------------------------------------------------------------------------------
def load_histogram(
    name,
    fileName,
    lumisections=[0,],
    run=0, 
    prefix="", suffix="", 
    MassBin=50, MassMin=66, MassMax=116, 
    pileup=False
):
    """
    load 2D histograms, project the specified lumisections to the mass axis: 
    - if pileup=False with (mass, lumisections); rebinning possible
    - if pileup=True, (PU, lumisections)
    
    Parameters
    ----------
    name : str
        name of the histogram to load
    fileName : str
        file where the histogram is stored
    lumisections : list
        list of lumisections that are taken from the histogram
    run : integer
        run number 
    prefix/suffix : str
        prefix for naming         
    MassBin/MassMin/MassMax : int
        For rebinning, Number of bins / Lower bound / Upper bound 
    pileup : Boolean
        If the pileup histogram is to be returned
    """

    import ROOT 
    
    file_ = ROOT.TFile(fileName)
    
    h_X_ls = file_.Get("{0}{1}".format(prefix, name)).Clone("{0}{1}{2}".format(prefix, name, suffix))
    h_X_ls.SetDirectory(0)
    h_X = h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}".format(name, run, suffix), lumisections[0], lumisections[0], "e")
    
    for ls in lumisections[1:]:
        h_X.Add(h_X_ls.ProjectionY("h_tmp_{0}_{1}_{2}_{3}".format(name, run, ls, suffix), ls, ls, "e"))
    
    if pileup:
        h_X.SetDirectory(0)
        return h_X
        
    # create new histogram in correct bin range
    hNew = ROOT.TH1D("h_mass_{0}_{1}_{2}_{3}".format(
        name, run, lumisections[0], suffix), "",MassBin, MassMin, MassMax)
    
    for ibin in range(0, h_X.GetNbinsX()+1):
        binCenter = h_X.GetBinCenter(ibin)
        if binCenter < MassMin:
            continue
        elif binCenter > MassMax:
            break
        else:
            content = h_X.GetBinContent(ibin)
            newBin = hNew.FindBin(binCenter)
            hNew.SetBinContent(newBin, hNew.GetBinContent(newBin) + content)
        
    hNew.SetDirectory(0)

    return hNew

# ------------------------------------------------------------------------------
def get_ls_for_next_measurement(
    lumisections, luminosities=None, zcounts=None, 
    lumiPerMeasurement=20, lsPerMeasurement=100, 
    threshold = 0.02
):
    """
    generator that takes the set of lumisections that are process 
    and yields slizes of lists of these lumisections 
    that should be used in the next measurement. 
    It can be run in two modes: 
    - If the parameter `luminosities` is not specified, the number of lsPerMeasurement is taken as criteria
    - If the parameter `luminosities` is specified, the amount of lumiPerMeasurementis taken as criteria
    
    Parameters
    ----------
    lumisections : list
        The list of lumisections that are going to be processed
    luminosities : list, optional
        The list of luminosity in \pb for each lumisection. 
    zcounts : list, optional
        The list of z boson counts for each lumisection.     
    lumiPerMeasurement : float, optional
        The amount of luminosity in \pb required for a measurement
    lsPerMeasurement : int, optional
        The number of lumisections required for a measurement
    threshold : float, optional
        If the luminosity in one lumisection is above this value and the number z counts is zero, the ls is skipped 
    """
        
    if luminosities:
        # make measurement based on number of lumisections
        if len(lumisections) != len(luminosities):
            print("ERROR:  === Same length of lumisections and luminosities is required!")
        
    while len(lumisections) > 0:
        
        # merge data to one measuement if remaining luminosity is too less for two measuements
        if luminosities:
            mergeMeasurements_ = sum(luminosities) < 1.5 * lumiPerMeasurement
        else:
            mergeMeasurements_ = len(lumisections) < 1.5 * lsPerMeasurement        
        
        recLumi_ = 0
        # produce list_good_ls_ with lumisections that are used for one measurement
        list_good_ls_ = []
        while len(lumisections) > 0:
            
            if luminosities and zcounts:
                # consider lumisections where we would expect to have at least any z count 
                #   (for lumi > 0.01 /pb we expect 0.01*500 = 5 Z bosons, the probability to find 0 is < 1%)
                #   (for lumi > 0.02 /pb we expect 0.02*500 = 10 Z bosons, the probability to find 0 is < 0.01%)
                # sort out lumisections without any Z candidate (maybe trigger was off)
                if luminosities[0] > threshold and zcounts[0] == 0:
                    print("WARNING:  === Zero Z boson candidates found {0}/pb while we would expect {1} -> skip lumi section {2}".format(luminosities[0], luminosities[0]*500, lumisections[0]))
                    del lumisections[0]
                    del luminosities[0]
                    del zcounts[0]
                    continue

            list_good_ls_.append(lumisections[0])
            del lumisections[0]
            
            if zcounts:
                del zcounts[0]
                
            if luminosities:
                recLumi_ += luminosities[0]
                del luminosities[0]

                
                # if we have collected enough luminosity
                if not mergeMeasurements_ and recLumi_ >= lumiPerMeasurement:
                    break
            else:
                # if we have collected enough luminosity sections
                if not mergeMeasurements_ and len(list_good_ls_) >= lsPerMeasurement:
                    break
                
        yield list_good_ls_

# ------------------------------------------------------------------------------
def getCorrelation(hPV_data, correlationsFileName, which, region="I"):
    """
    calculate the correlation factor from a resource file by folding the PV histogram with the correlation histogram  

    Parameters
    ----------
    hPV_data : TH1
        1D histogram with number of primary vertices in data
    correlationsFileName : str
        path to the file with correlation factors as function of number of primary vertices
    which: str
        which correlation factor to be read, can be: "IO", "ID", "Acceptance"
    region: str
        eta region for the correlation factor, can be: "I", "BB", "BE", "EE"
    """
    import numpy as np
    import ROOT
    ROOT.gROOT.SetBatch(True) # disable root prompts

    tfile = ROOT.TFile(correlationsFileName,"READ")

    # normalize pileup histogram
    hPV_data.Scale(1./hPV_data.Integral())
    avgPV = hPV_data.GetMean()
    hC = tfile.Get("C{0}_{1}".format(which, region))

    hC.Multiply(hPV_data)

    corr = hC.Integral()
    # if which=="IO":
    #     # load histograms
    #     truth = tfile.Get("truth_{0}".format(region))
    #     Sta = tfile.Get("Sta_{0}".format(region))
    #     Trk = tfile.Get("Trk_{0}".format(region))
    #     StaPass = tfile.Get("StaPass_{0}".format(region))

    #     # fold with pileup histogram
    #     truth.Multiply(hPV_data)
    #     Sta.Multiply(hPV_data)
    #     Trk.Multiply(hPV_data)
    #     StaPass.Multiply(hPV_data)

    #     # calculate correlation factor
    #     nTruth = truth.Integral()
    #     nSta = Sta.Integral()
    #     nTrk = Trk.Integral()
    #     nStaPass = StaPass.Integral()

    #     corr = (nStaPass / nSta) / (nTrk / (2*nTruth))

    # elif which=="ID":
    #     # load histograms
    #     truth = tfile.Get("truth_{0}".format(region))
    #     ID2 = tfile.Get("2ID_{0}".format(region))
    #     ID1 = tfile.Get("1ID_{0}".format(region))

    #     # fold with pileup histogram
    #     truth.Multiply(hPV_data)
    #     ID2.Multiply(hPV_data)
    #     ID1.Multiply(hPV_data)

    #     # calculate correlation factor
    #     nTruth = truth.Integral()
    #     nID1 = ID1.Integral()
    #     nID2 = ID2.Integral()

    #     corr = 4 * nTruth * nID2 / (nID1 + 2*nID2)**2

    tfile.Close()
    print("Correlation coefficienct c{0}_{1} = {2}".format(which, region, corr))
    print("For average primary vertices of <pv> = {0}".format(avgPV))
    return corr

# ------------------------------------------------------------------------------
def writeSummaryCSV(outCSVDir, outName="Mergedcsvfile", writeByLS=True, keys=None):
    """
    Collect "by LS" (and "by measurement") csv files and write one big csv file
    
    Parameters
    ----------
    outCSVDir : str
        Directory where the csv files of each run are located
    writeByLS : boolean, optional
        Boolean if a summary should be created for csv files that contain the Z rates for each lumisection
    """

    import logging as log
    import pandas as pd
    import glob
    
    print("INFO:  === Writing overall CSV file")
    rateFileList = sorted(glob.glob(outCSVDir + '/csvfile??????.csv'))
    df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True, sort=False)
    
    if keys:
        df_merged = df_merged[keys]
    
    with open(outCSVDir + '/' + outName + '_perMeasurement.csv', 'w') as file:
        df_merged.to_csv(file, index=False)

    if writeByLS:
        print("INFO:  ===Writing overall CSV file per LS")
        rateFileList = sorted(glob.glob(outCSVDir + '/csvfile*_*.csv'))
        csvList = []
        # add measurement label to the csv list
        for m in rateFileList:
            csv = pd.read_csv(m)
            measurement = int(m.split("_")[-1][:-4])
            fill = int(m.split("_")[-2].split("csvfile")[-1])
            csv["measurement"] = measurement
            csvList.append(csv)
        df_merged = pd.concat(csvList, ignore_index=True, sort=False)

        # df_merged = pd.concat([pd.read_csv(m) for m in rateFileList], ignore_index=True)

        with open(outCSVDir + '/' + outName + '_perLS.csv', 'w') as file:
            df_merged.to_csv(file, index=False)

# ------------------------------------------------------------------------------
def getEra(run):        
    """
    return the era for a given run number
    
    Parameters
    ----------
    run : int
        Drun number
    """

    if run <= 271658:
        return "Run2016A"
    if run <= 275376:
        return "Run2016B"        
    if run <= 276283:
        return "Run2016C"
    if run <= 276811:
        return "Run2016D"
    if run <= 277420:
        return "Run2016E"
    if run <= 278808:
        return "Run2016F"        
    if run <= 280385:
        return "Run2016G"
    if run <= 284044:
        return "Run2016H"        
    if run <= 297019:
        return "Run2017A"
    if run <= 299329:
        return "Run2017B"        
    if run <= 302029:
        return "Run2017C"
    if run <= 303434:
        return "Run2017D"
    if run <= 304826:
        return "Run2017E"
    if run <= 306462:
        return "Run2017F"        
    if run <= 306826:
        return "Run2017G"
    if run <= 307082:
        return "Run2017H"           
    if run <= 316995:
        return "Run2018A"
    if run <= 319312:
        return "Run2018B"        
    if run <= 320393:
        return "Run2018C"
    if run <= 325273:
        return "Run2018D"     
    if run <= 355769:
        return "Run2022B"   
    if run <= 357482:
        return "Run2022C"
    if run <= 357900:
        return "Run2022D"

    return "Run2022" #default 

def open_workspace(directory, filename, m):
    # loading root workspace
    file_result = f"{directory}/{filename}_{m}.root"

    if not os.path.isfile(file_result):
        log.warning("No result for `{0}`".format(file_result))
        return None, None

    f = ROOT.TFile(file_result,"READ")
    w = f.Get("workspace")
    return f, w

def open_workspace_yield(directory, filename, m, signal_parameters=False):
    # reading fit result from root workspace
    #   if signal_parameters is true, return the mean and standard deviation of the resolution function

    f, w = open_workspace(directory, "workspace_yield_"+filename, m)

    r = f.Get("fitResult")

    # yields are allwed to be negative to obtain reasonable errors, but we cap them here
    if r.covQual() == 3:
        def cap(s):
            va = w.var(s).getVal()
            vs = w.var(s).getError()

            if va < 0:
                return unc.ufloat(0, max(0, vs+va))
            else:
                return unc.ufloat(va, vs)
        Nsig = cap("Nsig")
        Nbkg = cap("Nbkg")
    else:
        log.warning(f"Covariance matrix has status {r.covQual()} and no relyable uncertainties! Use sqrt(N) instead")
        Nsig = unorm(max(0, w.var("Nsig").getVal()))
        Nbkg = unorm(max(0, w.var("Nbkg").getVal()))

    chi2 = w.arg("chi2").getVal()

    if signal_parameters:
        mean = None
        sigma = None
        for p in w.allVars():

            if p.GetName().startswith("sig_mean"):
                mean = p.getVal()
            if p.GetName().startswith("sig_sigma"):
                sigma = p.getVal()

            if mean is not None and sigma is not None:
                break

    f.Close()
    f.Delete()
    w.Delete()
    return Nsig, Nbkg, chi2, mean, sigma


# ------------------------------------------------------------------------------
def chart_to_js(_chart, _name, _data="data.csv"):
    import json, os
    
    chart_json = _chart.to_dict()
    # delete in line data
    del chart_json["datasets"]
    # add external data resource
    chart_json["data"] = {"url": _data}
    with open(_name+'.js','a') as file:
        file.write("""(function(vegaEmbed) {
    var spec = """)
        json.dump(chart_json, file, indent=4)
        file.write(""";
    var embedOpt = {\"renderer\": \"svg\", \"mode\": \"vega-lite\"};
    function showError(el, error){
        el.innerHTML = ('<div class=\"error\" style=\"color:red;\">'
        + '<p>JavaScript Error: ' + error.message + '</p>'
        + \"<p>This usually means there's a typo in your chart specification. \"
        + \"See the javascript console for the full traceback.</p>\"
        + '</div>');
        throw error;
    }""")
        file.write("const el = document.getElementById('vis');\n")
        file.write("vegaEmbed(\"#vis\", spec, embedOpt).catch(error => showError(el, error));\n")
        file.write("})(vegaEmbed);\n")

def chart_to_html(_chart, _name, **kwargs):

    with open(_name+'.js','w') as file:
        file.write("""
<!DOCTYPE html>
<html>
<head>
  <script type="text/javascript" src="https://cdn.jsdelivr.net/npm//vega@5"></script>
  <script type="text/javascript" src="https://cdn.jsdelivr.net/npm//vega-lite@4.8.1"></script>
  <script type="text/javascript" src="https://cdn.jsdelivr.net/npm//vega-embed@6"></script>
</head>
<body>
  <div id="vis"></div>
  <script type="text/javascript">
    """)
    chart_to_js(_chart, _name, **kwargs)
    with open(_name+'.js','a') as file:
        file.write("""
  </script>
</body>
</html>
""")

    os.system("mv {0}.js {0}.html".format(_name))
