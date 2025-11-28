#!/bin/bash
# Arguments: $1 = run number

RUN=$1
#INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/Muon_Run2022F/Muon/crab_Muon_Run2022F/*/*/*"
#INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/SingleMuon_Run2022C/SingleMuon/crab_SingleMuon_Run2022C/*/*/*"
#INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/Muon_Run2022C/Muon/crab_Muon_Run2022C/*/*/*"
#INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/Muon_Run2022D_v1/Muon/crab_Muon_Run2022D_v1/*/*/*"
#INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/Muon_Run2022E/Muon/crab_Muon_Run2022E/*/*/*"
INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/Muon_Run2022F/Muon/crab_Muon_Run2022F/230215_094428/*/*"
#INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/Muon_Run2022G/Muon/crab_Muon_Run2022G/*/*/*"
#INPUT_PATH="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Data/crab/V04/*/*/*/*/*/*"
OUTPUT_DIR="/eos/cms/store/group/comm_luminosity/ZCounting/2022/Test/Oct25_histos/"
#OUTPUT_DIR="/eos/user/t/tmenezes/ZCounting_Data/Nov27_TestInputHistos/"

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

cd /afs/cern.ch/user/t/tmenezes/work/private/CMSSW_12_4_18/src
eval `scramv1 runtime -sh`

cd "$_CONDOR_SCRATCH_DIR"



if [[ "$RUN" == "NONE" || -z "$RUN" ]]; then
    python3 /afs/cern.ch/user/c/cmszcont/CMSSW_12_4_18/src/ZCounting/ZCountAnalyze/python/Analyze/histogramize.py -i "$INPUT_PATH" -o "$OUTPUT_DIR"
else
    python3 /afs/cern.ch/user/c/cmszcont/CMSSW_12_4_18/src/ZCounting/ZCountAnalyze/python/Analyze/histogramize.py -i "$INPUT_PATH" -o "$OUTPUT_DIR" --run "$RUN"
fi
