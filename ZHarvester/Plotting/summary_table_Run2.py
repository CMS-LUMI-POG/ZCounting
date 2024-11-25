import json
import numpy as np
import argparse
import os
import pdb

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", required=True, type=str, help="input json from make_summary.py")
parser.add_argument("-s","--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

with open(args.input, "r") as file_info:
    info = json.load(file_info)

# labels for table header
labels = [{
    "2016": r"$\delta\NZ_{\mathrm{2016}}$",
    "2016preVFP": r"$\delta\NZ_{\mathrm{2016 preVFP}}$",
    "2016postVFP": r"$\delta\NZ_{\mathrm{2016 postVFP}}$",
    "2016H": r"$\delta\NZ_{\mathrm{2016 H}}$",
    "2017": r"$\delta\NZ_{\mathrm{2017}}$",
    "2017H": r"$\delta\NZ_\lowPU$",
    "2018": r"$\delta\NZ_{\mathrm{2018}}$"
},
{
    "2016": r"$\delta \frac{\NZ_\mathrm{2016}}{\NZ_\lowPU}$ }",
    "2016preVFP": r"$\delta \frac{\NZ_\mathrm{2016 preVFP}}{\NZ_\lowPU}$ }",
    "2016postVFP": r"$\delta \frac{\NZ_\mathrm{2016 postVFP}}{\NZ_\lowPU}$ }",
    "2016H": r"$\delta \frac{\NZ_\mathrm{2016 H}}{\NZ_\lowPU}$ }",
    "2017": r"$\delta \frac{\NZ_\mathrm{2017}}{\NZ_\lowPU}$ }",
    "2018": r"$\delta \frac{\NZ_\mathrm{2018}}{\NZ_\lowPU}$ }"    
}]

# factor = 0 means we get total uncertainty on each period
# factor = 1 means we get uncertainty for the ratio to 2017H
for factor in (0, 1):

    keys = []
    nominal = []
    stat = []           # statistical uncertainty
    bkgUp = []          # background subtraction
    bkgDown = []
    cHLTUp = []             # pileup dependent di-muon HLT MC correction
    cHLTDown = []           # pileup dependent di-muon HLT MC correction
    cIDUp = []              # pileup dependent di-muon ID MC correction
    cIDDown = []            # pileup dependent di-muon ID MC correction
    cIOUp = []              # pileup dependent muon inner/outer track MC correction
    cIODown = []            # pileup dependent muon inner/outer track MC correction
    cAcceptanceUp = []      # pileup dependent MC correction for bias due to kinematic selection of tracks, including the correction to parton level
    cAcceptanceDown = []    # pileup dependent MC correction for bias due to kinematic selection of tracks, including the correction to parton level
    prefireStatUp = []
    prefireStatDown = []
    prefireSystUp = []
    prefireSystDown = []
    prefireECALUp = []
    prefireECALDown = []
    
    base = "/".join(args.input.split("/")[:-2])
    suffix = args.input.split("/")[-1]
    def load(name):
        if not os.path.isfile(base+name+suffix):
            return None
        with open(base+name+suffix, "r") as file_info:
            info = json.load(file_info)
        return info
        
    # ---  variations from alternative fits
    altSigGen = []
    altSig1 = [] 
    altSig2 = []
    altBkg = []
    lumiUp = []
    lumiDown = []
    massUp = []
    massDown = []
    binWidthUp = []
    binWidthDown = []

    altBkgPass = []
    altBkgFail = []

    altSig1eff = [] 
    altSig2eff = []
    altSig1yld = [] 
    altSig2yld = []

    # load alternative infos
    prefix = "summary"
    info_altSigGen    = load(f"/{prefix}_altSigGen/")
    info_altSig1      = load(f"/{prefix}_altSigMCxCB/")
    info_altSig2      = load(f"/{prefix}_altSigMC/")
    info_altBkg       = load(f"/{prefix}_altBkg/")
    info_lumiUp       = load(f"/{prefix}_lumi30/")
    info_lumiDown     = load(f"/{prefix}_lumi15/")
    info_massUp       = load(f"/{prefix}_mass50to130/")
    info_massDown     = load(f"/{prefix}_mass70to110/")
    info_binWidthUp   = load(f"/{prefix}_binWidth1/")
    info_binWidthDown = load(f"/{prefix}_binWidth025/")

    # sortdict = {"2016preVFP":1, "2016postVFP":2, "2016H":3, "2017":4, "2017H":5, "2018":6}
    # sortdict = {"2016preVFP":1, "2016postVFP":2, "2017":4, "2017H":5, "2018":6}
    sortdict = {"2016":1, "2017":4, "2017H":5, "2018":6}

    items = filter(lambda x: x[0] in sortdict.keys(), info.items())

    # for alternative signal - only take the effect on the efficiency
    for key, iEra in sorted(items, key=lambda item: sortdict[item[0]]):

        # for key, iEra in sorted(items, key=lambda item: sortdict[item[0]]):
        # in case of empty info, continue
        if iEra['recZCount'] == 0:
            continue
        # if we compute the ratio with 2017H, we don't do it for 2017H :)
        if factor and key == "2017H":
            continue
            
        keys.append(key)
        # nominal.append(iEra['prefECAL_nominal'])
        # stat.append( 1./iEra['recZCount'] 
        #     + factor*1./info['2017H']['recZCount'] )
        stat.append( (iEra['recZCount_err']/iEra['recZCount'] )**2
            + factor * (info['2017H']['recZCount_err']/info['2017H']['recZCount'])**2 
        )
        # bkgUp.append( (iEra['zRec_bkgUp']-iEra['zRec'])/iEra['zRec'] 
        #     - factor*(info['2017H']['zRec_bkgUp']-info['2017H']['zRec'])/info['2017H']['zRec'] )  # uncertainty on background subtraction
        # bkgDown.append( (iEra['zRec_bkgDown']-iEra['zRec'])/iEra['zRec'] 
        #     - factor*(info['2017H']['zRec_bkgDown']-info['2017H']['zRec'])/info['2017H']['zRec'] )

        for up, down, name in (
            (cHLTUp, cHLTDown, "cHLT"),
            (cIDUp, cIDDown, "cID"),
            (cIOUp, cIODown, "cIO"),
            (cAcceptanceUp, cAcceptanceDown, "cAcceptance"),
        ):
            up.append((iEra['recZCount_{0}Up'.format(name)]-iEra['recZCount'])/iEra['recZCount']
                - factor*(info['2017H']['recZCount_{0}Up'.format(name)]-info['2017H']['recZCount'])/info['2017H']['recZCount'] )
            down.append((iEra['recZCount_{0}Down'.format(name)]-iEra['recZCount'])/iEra['recZCount']
                - factor*(info['2017H']['recZCount_{0}Down'.format(name)]-info['2017H']['recZCount'])/info['2017H']['recZCount'] )

        # --- Muon prefire corrections    
        # correlated between 2017 and 2017H because it's the same statistics used in the measurement of the prefire correction
        if key == "2017": 
            prefireStatUp.append((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
                - factor*((info['2017H']['prefMuon_StatUp']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal']) )
            prefireStatDown.append((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
                - factor*((info['2017H']['prefMuon_StatDown']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal']) )
        else:
            prefireStatUp.append(np.sqrt( ((iEra['prefMuon_StatUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])**2
                + factor*((info['2017H']['prefMuon_StatUp']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal'])**2 ))
            prefireStatDown.append(-np.sqrt( ((iEra['prefMuon_StatDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal'])**2
                + factor*((info['2017H']['prefMuon_StatDown']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal'])**2 ))
        
        prefireSystUp.append((iEra['prefMuon_SystUp']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
            - factor*(info['2017H']['prefMuon_SystUp']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal']
        )
        prefireSystDown.append((iEra['prefMuon_SystDown']-iEra['prefMuon_nominal'])/iEra['prefMuon_nominal']
            - factor*(info['2017H']['prefMuon_SystDown']-info['2017H']['prefMuon_nominal'])/info['2017H']['prefMuon_nominal']
        )
        prefireECALUp.append((iEra['prefECAL_Up']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal']
            - factor*(info['2017H']['prefECAL_Up']-info['2017H']['prefECAL_nominal'])/info['2017H']['prefECAL_nominal']
        )
        prefireECALDown.append((iEra['prefECAL_Down']-iEra['prefECAL_nominal'])/iEra['prefECAL_nominal']
            - factor*(info['2017H']['prefECAL_Down']-info['2017H']['prefECAL_nominal'])/info['2017H']['prefECAL_nominal']
        )

        if info_altBkg:
            altBkgPass.append((info_altBkg[key]['recZCount_altBkgPass']-iEra['recZCount'])/iEra['recZCount']
                - factor*(info_altBkg['2017H']['recZCount_altBkgPass']-info['2017H']['recZCount'])/info['2017H']['recZCount'])
            altBkgFail.append((info_altBkg[key]['recZCount_altBkgFail']-iEra['recZCount'])/iEra['recZCount']
                - factor*(info_altBkg['2017H']['recZCount_altBkgFail']-info['2017H']['recZCount'])/info['2017H']['recZCount'])        

        
        # --- from alternative sources
        for src, arr in (
            (info_altSigGen, altSigGen),
            (info_altSig1, altSig1),
            (info_altSig2, altSig2),
            (info_altBkg, altBkg),
            (info_lumiUp, lumiUp),
            (info_lumiDown, lumiDown),
            (info_massUp, massUp),
            (info_massDown, massDown),
            (info_binWidthUp, binWidthUp),
            (info_binWidthDown, binWidthDown),
        ):
            if src:
                arr.append((src[key]['recZCount']-iEra['recZCount'])/iEra['recZCount']
                - factor*(src['2017H']['recZCount']-info['2017H']['recZCount'])/info['2017H']['recZCount'])

        # split signal shape variations in sources where inner track and where outer track are varied
        for src, arr1, arr2 in (
            (info_altSig1, altSig1eff, altSig1yld),
            (info_altSig2, altSig2eff, altSig2yld),
        ):  
            if src:
                # only take effect on efficiency
                alt = iEra['selZCount'] / src[key]['selZCount'] * src[key]['recZCount']
                altLow = info['2017H']['selZCount'] / src['2017H']['selZCount'] * src['2017H']['recZCount']
                arr1.append((alt - iEra['recZCount'])/iEra['recZCount']
                    - factor*(altLow - info['2017H']['recZCount'])/info['2017H']['recZCount'])

                # only take effect on other parts
                arr2.append((src[key]['selZCount']-iEra['selZCount'])/iEra['selZCount']
                    - factor*(src['2017H']['selZCount']-info['2017H']['selZCount'])/info['2017H']['selZCount'])

    suffix = "_ratio" if factor==1 else ""

    outname = args.input.split("/")[-1].replace("info","result"+suffix)
    outname = outname.replace(".json",".txt")
    # ---- make latex table with variations:
    with open(outDir+"/"+outname, "w") as outfile:

        columns = "l|"
        columns += "".join(["c" for c in range(len(sortdict.keys()) - factor)])
        outfile.write(r"\begin{tabular}{"+columns+"}"+"\n")
        
        if factor:
            outfile.write(" \multirow{3}{*}{ } & \multirow{3}{*}{ " + " & \multirow{3}{*}{ ".join([labels[factor][x] for x in keys]) + r" \\"+"\n")
            outfile.write(" & " + " & ".join([" " for x in keys]) + r" \\"+"\n")
            outfile.write(" & " + " & ".join([" " for x in keys]) + r" \\"+"\n")
        else:
            outfile.write(" & " + " & ".join([labels[factor][x] for x in keys]) + r" \\"+"\n")

        outfile.write(r" \hline "+"\n")

        # outfile.write("Nominal & "+" & ".join([str(abs(int(round(x,0)))) for x in nominal]) + r" \\" +"\n")
        #
        # outfile.write(r" \hline "+"\n")
        
        sysUp = [0 for x in stat]
        sysDown = [0 for x in stat]
        for name, values, connector in (
            ("Correlation \CHLT       ", cHLTUp,          r" & $\pm "),
            ("Correlation \CID        ", cIDUp,           r" & $\pm "),
            ("Correlation \CGlo       ", cIOUp,           r" & $\pm "),
            ("Acceptance              ", cAcceptanceUp,     r" & $\pm "),
            ("Background subtr.       ", bkgUp,           r" & $\pm "),
            ("Muon pref. up (stat.)   ", prefireStatUp,   r" & $"),
            ("Muon pref. down (stat.) ", prefireStatDown, r" & $"),
            ("Muon pref. up (syst.)   ", prefireSystUp,   r" & $"),
            ("Muon pref. down (syst.) ", prefireSystDown, r" & $"),
            ("ECAL pref. up           ", prefireECALUp,   r" & $"),
            ("ECAL pref. down         ", prefireECALDown, r" & $"),
            ("Alt. sig. model (MCxCB) ", altSig1,          r" & $"),
        #    ("Alt. sig. model (MC)    ", altSig2,          r" & $"),
            ("Alt. sig. model (Gen)   ", altSigGen,       r" & $"),
            ("Alt. bkg. model (Pass)  ", altBkgPass,      r" & $"),
            ("Alt. bkg. model (Fail)  ", altBkgFail,      r" & $"),
        #    ("Alt. sig. model (MCxCB) - eff ", altSig1eff,          r" & $"),
        #    ("Alt. sig. model (MC) - eff    ", altSig2eff,          r" & $"),
        #    ("Alt. bkg. model - eff        ", altBkgeff,          r" & $"),
        #    ("Alt. sig. model (MCxCB) - yield ", altSig1yld,          r" & $"),
        #    ("Alt. sig. model (MC) - yield    ", altSig2yld,          r" & $"),
        #    ("Alt. bkg. model - yield        ", altBkgyld,          r" & $"),
            ("Lumi slice up           ", lumiUp,          r" & $"),
            ("Lumi slice down         ", lumiDown,        r" & $"),
            ("Bin width up            ", binWidthUp,      r" & $"),
            ("Bin width down          ", binWidthDown,    r" & $"),
        #    ("Mass range up           ", massUp,        r" & $"),
        #    ("Mass range down         ", massDown,      r" & $"),

        ):
            if len(values) == 0:
                continue

            entries = []
            for i, x in enumerate(values):
                value = round(x*100,2)

                if "\pm" in connector:
                    if value < 0:
                        entries.append("\mp"+str(abs(value))+"$" )
                    else:
                        entries.append("\pm"+str(abs(value))+"$" )
                    sysUp[i] += x**2
                    sysDown[i] += x**2
                elif round(x*100,2) > 0:
                    entries.append("+"+str(value)+"$" )
                    sysUp[i] += x**2
                elif round(x*100,2) < 0:
                    entries.append(str(value)+"$" )
                    sysDown[i] += x**2
                else:
                    entries.append(r"\NA $")
            connector = connector.replace("\pm","")
            outfile.write(name + connector + connector.join(entries) + r" \\"+"\n")
        outfile.write(r" \hline "+"\n")

        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in stat]
        outfile.write("Statistical             & $\pm "+" & $\pm ".join(entries) + r" \\"+"\n")

        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in sysUp]
        outfile.write("Systematic up           & $+"+" & $+".join(entries) + r" \\"+"\n")
        
        entries = [str(round(abs(np.sqrt(x))*100,2))+"$" for x in sysDown]
        outfile.write("Systematic down         & $-"+" & $-".join(entries) + r" \\"+"\n")

        outfile.write(r" \hline "+"\n")
        outfile.write(r" \hline "+"\n")
        
        entries = [str(round(abs(np.sqrt(x + y))*100,2))+"$" for x, y in zip(stat, sysUp)]
        outfile.write("Total up                & $+"+" & $+".join(entries) + r" \\"+"\n")
        
        entries = [str(round(abs(np.sqrt(x + y))*100,2))+"$" for x, y in zip(stat, sysDown)]
        outfile.write("Total down              & $-"+" & $-".join(entries) + r" \\"+"\n")

        outfile.write(r"\end{tabular}"+"\n")
