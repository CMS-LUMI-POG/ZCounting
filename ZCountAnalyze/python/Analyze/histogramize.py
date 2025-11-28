#!/usr/bin/env python3
#
################################################################################
#
# create histograms from ntuples used for extraction of the Z counts
#
################################################################################

import uproot 
import awkward as ak
import numpy as np
import pdb
import vector
import hist
import os

import argparse
parser = argparse.ArgumentParser(prog='./histogramize')
parser.add_argument(
    '-i', '--input', nargs='+',
    help='specify input ntuple root files'
)
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output directory'
)
parser.add_argument(
    '--runs', nargs='*', default=None,
    help='specify runs to be processed'
)
args = parser.parse_args()

### settings
mode = "LUM-21-001" # "LUM-21-001" to create histograms as done in LUM-21-001 or "cmssw" to create the histograms done in the cmssw plugin

# acceptance cuts
ptCut = 25
etaCut = 2.4

# mass binning
mass_lo = 60
mass_hi = 120
mass_nBins = 4*(int)(mass_hi - mass_lo)

# lumisection binning
lumi_nBins = 2500
lumi_lo = 0.5
lumi_hi = 2500.5

# binning for histogram of number of primary vertices (needed for correlation factor)
pv_nBins = 100
pv_lo = 0.5
pv_hi = 100.5

muonMass = 0.105658369

etaBound = 0.9  # below this value is barrel (B) above is endcap (E)

### resources
treename = "zcounting/tree"

filenames = [f+":"+treename for f in args.input]
dirOut = args.output

# New DEBUGs
#print(f"[INFO] Input files: {args.input}")
print(f"[INFO] Output directory: {dirOut}")
print(f"[INFO] Mode: {mode}")
print(f"[INFO] Will process runs: {args.runs if args.runs else 'Auto-detected from files'}")




if not os.path.isdir(dirOut):
    print(f"create output directory {dirOut}")
    os.mkdir(dirOut)


branches_muon = [
    "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", 
    "Muon_nPixelHits", "Muon_nTrackerLayers", "Muon_trackAlgo",
    "Muon_ID", "Muon_triggerBits"
    ]
branches_standalones = [
    "Muon_ptStaReg", "Muon_etaStaReg", "Muon_phiStaReg", "Muon_chargeStaReg", 
    "Muon_ptStaUpd", "Muon_etaStaUpd", "Muon_phiStaUpd", "Muon_chargeStaUpd", 
    "Muon_useUpdated", "Muon_ID",
    "Muon_phi"
    ]
branches_tracks = [
    "Track_pt", "Track_eta", "Track_phi", "Track_charge",
    "Track_nPixelHits", "Track_nTrackerLayers", "Track_trackAlgo" 
    ]
branches_event = ["lumiBlock", "runNumber", "nPV"]
aliases = {

}

runs = args.runs

if runs == None:
    # scan files for uns
    print(f"Scan files for runs")
    runs = set()
    for batch in uproot.iterate(filenames, filter_name="runNumber", library="np"):
        runs.update(set(batch["runNumber"]))
        # New DEBUGs
        print(f"[INFO] Found runs: {sorted(runs)}")

# loop over all runs
print(f"Loop over all runs")
print(runs)

for run in runs:
    print(f"Now at run {run}")
    hists = dict()
    # process the events in batches to avoid using too much memory
    for batch in uproot.iterate(filenames, step_size="100 MB", library="ak", 
        filter_name=branches_event+branches_muon+branches_tracks+branches_standalones, 
        aliases=aliases, cut=f"runNumber == {run}"
    ):
        if len(batch) == 0:
            continue
        print(f"Have a batch with {len(batch)} events")    
        
        events = batch[branches_event]

        # histograms with number of primary vertices vs lumisection
        xx = events["lumiBlock"].to_numpy()
        yy = events["nPV"].to_numpy()
        
        hist = np.histogram2d(xx, yy, bins=[lumi_nBins, pv_nBins], range=((lumi_lo,lumi_hi), (pv_lo, pv_hi)))
        if "h_npv" in hists.keys():
            hists["h_npv"] = (hist[0] + hists["h_npv"][0], hist[1], hist[2])
        else:
            hists["h_npv"] = hist

        # muon objects
        muons = batch[branches_muon]
        # rename
        muons["pt"] = muons["Muon_pt"]
        muons["eta"] = muons["Muon_eta"]
        muons["phi"] = muons["Muon_phi"]
        muons["charge"] = muons["Muon_charge"]

        # index of the muon within the event (needed to disambiguate)
        muons["index"] = ak.local_index(muons["Muon_pt"])

        # collection of muons that pass ID and hlt
        hlt_pass = muons[(muons["pt"] > ptCut) & (abs(muons["eta"]) < etaCut) & (muons["Muon_ID"]>=4) & (muons["Muon_triggerBits"]&1 != 0)]

        # collection of muons that pass ID and fail hlt
        hlt_fail = muons[(muons["pt"] > ptCut) & (abs(muons["eta"]) < etaCut) & (muons["Muon_ID"]>=4) & (muons["Muon_triggerBits"]&1 == 0)]

        # collection of muons that pass global and fail ID
        sel_fail = muons[(muons["pt"] > ptCut) & (abs(muons["eta"]) < etaCut) & (muons["Muon_ID"]==3)]

        if mode == "cmssw":
            # collection of muons that pass standalone and fail global
            glo_fail = muons[(muons["pt"] > ptCut) & (abs(muons["eta"]) < etaCut) & (muons["Muon_ID"]==2)]
        elif mode == "LUM-21-001":
            sta_pass = muons[(muons["pt"] > ptCut) & (abs(muons["eta"]) < etaCut) 
                & (muons["Muon_nPixelHits"] > 0 ) & (muons["Muon_nTrackerLayers"] > 5)
                & (muons["Muon_trackAlgo"] != 13) & (muons["Muon_trackAlgo"] != 14) # veto muon seeded tracks
                ]

            # collection of muons using the standalone track 
            stas = batch[branches_standalones]
            stas["pt"] = stas["Muon_useUpdated"] * stas["Muon_ptStaUpd"] + (1 - stas["Muon_useUpdated"]) * stas["Muon_ptStaReg"] 
            stas["eta"] = stas["Muon_useUpdated"] * stas["Muon_etaStaUpd"] + (1 - stas["Muon_useUpdated"]) * stas["Muon_etaStaReg"] 
            stas["phi"] = stas["Muon_useUpdated"] * stas["Muon_phiStaUpd"] + (1 - stas["Muon_useUpdated"]) * stas["Muon_phiStaReg"] 
            stas["charge"] = stas["Muon_useUpdated"] * stas["Muon_chargeStaUpd"] + (1 - stas["Muon_useUpdated"]) * stas["Muon_chargeStaReg"] 
            stas["index"] = muons["index"]

            glo_pass = stas[(stas["pt"] > ptCut) & (abs(stas["eta"]) < etaCut) & (stas["Muon_ID"]>=3)]
            glo_fail = stas[(stas["pt"] > ptCut) & (abs(stas["eta"]) < etaCut) & (stas["Muon_ID"]==2)]


        # collection of track objects - needed for global efficiency
        trks = batch[branches_tracks]
        # rename
        trks["pt"] = trks["Track_pt"]
        trks["eta"] = trks["Track_eta"]
        trks["phi"] = trks["Track_phi"]
        trks["charge"] = trks["Track_charge"]
        # select good tracks within acceptance
        trks = trks[(trks["pt"] > ptCut) & (abs(trks["eta"]) < etaCut) & (trks["Track_nPixelHits"] > 0) & (trks["Track_nTrackerLayers"] > 5)]

        if mode == "LUM-21-001":
            trks = trks[(trks["Track_trackAlgo"] != 13) & (trks["Track_trackAlgo"] != 14)] # veto muon seeded tracks


        def produce(name, tags, probes=None, disambiguate=False):
            # build pairs from the specified collections, compute masses, and fill histograms

            # 1.) build pairs
            if probes is None:
                # if no probes are specified we do combinations among the tags (i.e. for HLT2 categorie)

                # require pairs with opposite charge
                p_charge = ak.combinations(tags["charge"], 2)        
                lefts, rights = ak.unzip(p_charge)
                mask_p = lefts*rights == -1

                p_pt = ak.combinations(tags["pt"], 2)[mask_p]
                p_eta = ak.combinations(tags["eta"], 2)[mask_p]
                p_phi = ak.combinations(tags["phi"], 2)[mask_p]
            else:
                # require pairs with opposite charge
                p_charge = ak.cartesian([tags["charge"], probes["charge"]])        
                l_charge, r_charge = ak.unzip(p_charge)
                mask_p = l_charge*r_charge == -1


                if disambiguate:
                    # take care that the same object is not used as tag and probe, 
                    #    this can happen e.g. if the charge used for the probe is different then the one used for the tag
                    # s1 = sum(ak.flatten(mask_p))

                    p_idx = ak.cartesian([tags["index"], probes["index"]])
                    l_idx, r_idx = ak.unzip(p_idx)
                    mask_p = mask_p & (l_idx != r_idx)

                    # s2 = sum(ak.flatten(mask_p))

                    # print(f"Reject {s1-s2} pairs due to ambiguity")

                p_pt = ak.cartesian([tags["pt"],probes["pt"]])[mask_p]
                p_eta = ak.cartesian([tags["eta"],probes["eta"]])[mask_p]
                p_phi = ak.cartesian([tags["phi"],probes["phi"]])[mask_p]

            # 2.) calculate mass of each pair
            l_pt, r_pt = ak.unzip(p_pt)
            l_eta, r_eta = ak.unzip(p_eta)
            l_phi, r_phi = ak.unzip(p_phi)

            mu1 = vector.obj(pt=l_pt, phi=l_phi, eta=l_eta, mass=l_pt*0+muonMass)
            mu2 = vector.obj(pt=r_pt, phi=r_phi, eta=r_eta, mass=r_pt*0+muonMass)

            masses = (mu1 + mu2).mass

            mask_mass = (masses > mass_lo) & (masses < mass_hi)

            # 3.) Fill histograms    
            for region, mask_eta in (
                ("BB", (abs(l_eta) < etaBound) & (abs(r_eta) < etaBound)),
                ("BE", (abs(l_eta) < etaBound) != (abs(r_eta) < etaBound)),
                ("EE", (abs(l_eta) >= etaBound) & (abs(r_eta) >= etaBound))
            ):

                yy = masses[mask_eta & mask_mass]   # select tag and probe pairs in eta and mass range
                xx = events["lumiBlock"]

                # pair multiplicities in each event 
                # counts = ak.num(yy)
                
                # bring the event by event array into the same dimension
                if len(xx) > 0:
                    xx = ak.unflatten(xx,counts=1)
                
                # bring the arrays into the same form
                xx, yy = ak.broadcast_arrays(xx,yy)
                
                xx = ak.flatten(xx).to_numpy()
                yy = ak.flatten(yy).to_numpy()
                
                hist = np.histogram2d(xx, yy, bins=[lumi_nBins, mass_nBins], range=((lumi_lo,lumi_hi), (mass_lo, mass_hi)))
                
                histname = f"h_mass_{name}_{region}"
                if histname in hists.keys():
                    hists[histname] = (hist[0] + hists[histname][0], hist[1], hist[2])
                else:
                    hists[histname] = hist

        print("[DEBUG] Producing histograms for category: 2HLT")
        produce("2HLT", hlt_pass)
        print("[DEBUG] Producing histograms for category: 1HLT")
        produce("1HLT", hlt_pass, hlt_fail) 
        print("[DEBUG] Producing histograms for category: SIT_fail")
        produce("SIT_fail", hlt_pass, sel_fail) 

        if mode == "cmssw":
            produce("Glo_fail", hlt_pass, glo_fail) 
            produce("Glo_fail", hlt_pass, trks)
        elif mode == "LUM-21-001":
            # New DEBUGs
            print("[DEBUG] Producing histograms for mode: LUM-21-001")
            produce("Sta_pass", hlt_pass, sta_pass) 
            produce("Sta_fail", hlt_pass, trks)

            produce("Glo_pass", hlt_pass, glo_pass, disambiguate=True) 
            produce("Glo_fail", hlt_pass, glo_fail)            

    # open output root file
    print(f"Write out results in `output_Run{run}.root`")
    output = uproot.recreate(f"{dirOut}/output_Run{run}.root")

    # write the histograms out
    for name, hist in hists.items():
        print(f"[DEBUG] Writing histogram: {name}")
        output[name] = hist



    
