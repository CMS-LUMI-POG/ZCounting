import ROOT
import os,sys
from ROOT import TGraphAsymmErrors
from ROOT import TGraphErrors
from ROOT import TColor
from array import array
from ROOT import *
from operator import truediv
import random
import math
import pandas
import os.path
import logging as log
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b","--beginRun",help="first run to analyze [%default]",default=299918)
parser.add_argument("-e","--endRun",help="analyze stops when comes to this run [%default]",default=1000000)
parser.add_argument("-l","--lumiChunk",help="define statistics: measurement less than this to be merged with next measurement [%default]",default=1000000)
parser.add_argument("-s","--sizeChunk",help="define granularity: numbers of LS to be merged for one measurement [%default]",default=50)
parser.add_argument("-v","--verbose",help="increase logging level from INFO to DEBUG",default=False,action="store_true")

args = parser.parse_args()
if args.verbose:
    log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
else:
    log.basicConfig(format="%(levelname)s: %(message)s", level=log.INFO)

#ByLS csv inputs
inFile="/eos/cms/store/user/jsalfeld/2017LumiByLS_notrig_PU.csv"
#inFile="/eos/cms/store/user/jsalfeld/2017LumiByLS_trig.csv"

#Data inputs
eosDir="/eos/cms/store/group/comm_luminosity/ZCounting/DQMFiles2017/cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2017/SingleMuon/"

#MC inputs: to build MC*Gaussian template for efficiency fitting
mcDir="/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/LookupTable/"
mcShapeSubDir="MCFiles/92X_norw_IsoMu27_noIso/"

#Constant settings
currentYear=2017
chunkSize=50
maximumLS=2500
secPerLS=float(23.3)


log.info("Loading C marco...")
ROOT.gROOT.Macro( os.path.expanduser( '~/.rootlogon.C' ) )
ROOT.gROOT.LoadMacro("calculateDataEfficiency_v3.C")
ROOT.gROOT.LoadMacro("calculateZEfficiency.C")
ROOT.gROOT.SetBatch(True)

log.info("Loading input byls csv...")
lumiFile=open(inFile)
lumiLines=lumiFile.readlines()
data=pandas.read_csv(inFile, sep=',',low_memory=False, skiprows=[0,len(lumiLines)-5,len(lumiLines)-4,len(lumiLines)-3,len(lumiLines)-2,len(lumiLines)-1,len(lumiLines)])
log.debug("%s",data.axes)

# TAKE INPUT CSV FILE AND STRUCTURE PER-RUN BASIS, THEN CREATE LIST OF LUMI AND LS`s PER RUN
LSlist=data.groupby('#run:fill')['ls'].apply(list)
recLumiList=data.groupby('#run:fill')['recorded(/pb)'].apply(list)
delLumiList=data.groupby('#run:fill')['delivered(/pb)'].apply(list)
avgpuList=data.groupby('#run:fill')['avgpu'].apply(list)
timeList=data.groupby('#run:fill')['time'].apply(list)

for i in range(0,len(LSlist)):
	LSlist[i]=[int(x.split(':')[0]) for x in LSlist[i]]
fillRunlist=data.drop_duplicates('#run:fill')['#run:fill'].tolist()

log.debug("%s",fillRunlist)
log.debug("length LS list: %i",len(LSlist))
log.debug("length Run list: %i",len(fillRunlist))

log.info("Looping over runs...")
for run_i in range(0,len(fillRunlist)):

    run=int(fillRunlist[run_i].split(':')[0])
    fill=int(fillRunlist[run_i].split(':')[1])

    if run<int(args.beginRun) or run>=int(args.endRun):
        continue;
    if run<299918:#Z Coungig module is enabled since this run 
	continue

    if run<302030:
	era="C"
    elif run>=302030 and run<303434:
	era="D"
    elif run>=303434 and run<304910:
	era="E"
    elif run>=304910:
        era="F"
    else:
        log.error("===RunNum %s cannot be associated with an era",run)
        continue

    log.info("===Running Run %i",run)
    log.info("===Running Fill %i",fill)
        
    LSchunks 	= [LSlist[run_i][x:x+chunkSize] for x in range(0, len(LSlist[run_i]), chunkSize)]
    Del_chunks  = [delLumiList[run_i][x:x+chunkSize] for x in range(0, len(delLumiList[run_i]), chunkSize)]
    Rec_chunks  = [recLumiList[run_i][x:x+chunkSize] for x in range(0, len(recLumiList[run_i]), chunkSize)]
    Avgpu_chunks = [avgpuList[run_i][x:x+chunkSize] for x in range(0, len(avgpuList[run_i]), chunkSize)]
    time_chunks = [timeList[run_i][x:x+chunkSize] for x in range(0, len(timeList[run_i]), chunkSize)]

    log.debug("===Setting up arrays for output csv...")
    fillarray=array('d')
    beginTime=[]
    endTime=[]
    Zrate=array('d')
    instDel=array('d')
    lumiDel=array('d')
    ZyieldDel=array('d')

    ZyieldRec=array('d')
    lumiRec=array('d')
    windowarray=array('d')
    deadTime=array('d')
    beginLS=array('i')
    endLS=array('i')

    HLTeffB=array('d')
    HLTeffE=array('d')
    SITeffB=array('d')
    SITeffE=array('d')
    StaeffB=array('d')
    StaeffE=array('d')

    ZMCeff=array('d')
    ZMCeffBB=array('d')
    ZMCeffBE=array('d')
    ZMCeffEE=array('d')

    ZBBeff=array('d')
    ZBEeff=array('d')
    ZEEeff=array('d')

    nMeasurements=0
    skipStatus=0

    log.info("===Loading input DQMIO.root file...")
    eosFile = ""
    if os.path.isfile(eosDir+"000"+str(run)[:-2]+"xx/DQM_V0001_R000"+str(run)+"__SingleMuon__Run2017"+era+"-PromptReco-v1__DQMIO.root"):
        eosFile = eosDir+"000"+str(run)[:-2]+"xx/DQM_V0001_R000"+str(run)+"__SingleMuon__Run2017"+era+"-PromptReco-v1__DQMIO.root"
    elif os.path.isfile(eosDir+"000"+str(run)[:-2]+"xx/DQM_V0001_R000"+str(run)+"__SingleMuon__Run2017"+era+"-PromptReco-v2__DQMIO.root"):
        eosFile = eosDir+"000"+str(run)[:-2]+"xx/DQM_V0001_R000"+str(run)+"__SingleMuon__Run2017"+era+"-PromptReco-v2__DQMIO.root"
    elif os.path.isfile(eosDir+"000"+str(run)[:-2]+"xx/DQM_V0001_R000"+str(run)+"__SingleMuon__Run2017"+era+"-PromptReco-v3__DQMIO.root"):
        eosFile = eosDir+"000"+str(run)[:-2]+"xx/DQM_V0001_R000"+str(run)+"__SingleMuon__Run2017"+era+"-PromptReco-v3__DQMIO.root"
    else:
        log.warning("===Missing DQM files for Run%i",run)
	continue

    log.info("===Looping over LSchunks...")
    for chunk_i in range(0,len(LSchunks)):
        if skipStatus==2:
            break
        if float(LSchunks[chunk_i][-1]) > maximumLS:
            log.warning("======Losing data after LS %i for Run%i",maximumLS,run)
            skipStatus=2
            while LSchunks[chunk_i][-1] > maximumLS:
                del LSchunks[chunk_i][-1]
                del Del_chunks[chunk_i][-1]
                del Rec_chunks[chunk_i][-1]
                del Avgpu_chunks[chunk_i][-1]
                del time_chunks[chunk_i][-1]

        nMeasurements=nMeasurements+1

        log.info("======Running LSchunk No.%i",chunk_i)
        log.debug("======LS list: %s",LSchunks[chunk_i])

	recLumi_i = sum(Rec_chunks[chunk_i])
	delLumi_i = sum(Del_chunks[chunk_i])	
        deadtime_i = recLumi_i/delLumi_i
        log.debug("======RecLumi: %f",recLumi_i)
        log.debug("======DelLumi: %f",delLumi_i)
        log.debug("======DeadTime: %f",deadtime_i)


	avgPileup_i = sum(Avgpu_chunks[chunk_i])
	avgPileup_i = avgPileup_i/len(Avgpu_chunks[chunk_i])
        log.debug("======avgPU: %f",avgPileup_i)

	datestamp_low=time_chunks[chunk_i][0].split(" ")
	date_low=ROOT.TDatime(currentYear,int(datestamp_low[0].split("/")[0]),int(datestamp_low[0].split("/")[1]),int(datestamp_low[1].split(":")[0]),int(datestamp_low[1].split(":")[1]),int(datestamp_low[1].split(":")[2]))
        datestamp_up=time_chunks[chunk_i][-1].split(" ")
	date_up=ROOT.TDatime(currentYear,int(datestamp_up[0].split("/")[0]),int(datestamp_up[0].split("/")[1]),int(datestamp_up[1].split(":")[0]),int(datestamp_up[1].split(":")[1]),int(datestamp_up[1].split(":")[2]))
        timeWindow_i=(date_up.Convert()-date_low.Convert())+secPerLS
        log.debug("======time_chunks: %s",time_chunks[chunk_i])
        log.debug("======beginTime: %s",date_low.Convert())
        log.debug("======endTime: %s",date_up.Convert())
        log.debug("======timeWindow: %f",timeWindow_i)

        log.debug("Openning DQMIO.root file: %s", eosFile)
        HLTeffB_i=ROOT.calculateDataEfficiency_v3(0,str(eosFile),".",str(run),chunk_i,LSchunks[chunk_i][0],LSchunks[chunk_i][-1],avgPileup_i,"HLT",0,0,0,0,0,recLumi_i)
        HLTeffE_i=ROOT.calculateDataEfficiency_v3(0,str(eosFile),".",str(run),chunk_i,LSchunks[chunk_i][0],LSchunks[chunk_i][-1],avgPileup_i,"HLT",1,0,0,0,0,recLumi_i)
        SITeffB_i=ROOT.calculateDataEfficiency_v3(0,str(eosFile),".",str(run),chunk_i,LSchunks[chunk_i][0],LSchunks[chunk_i][-1],avgPileup_i,"SIT",0,2,1,2,1,recLumi_i,mcDir+mcShapeSubDir+"MuSITEff/MC/probes.root",mcDir)
        SITeffE_i=ROOT.calculateDataEfficiency_v3(0,str(eosFile),".",str(run),chunk_i,LSchunks[chunk_i][0],LSchunks[chunk_i][-1],avgPileup_i,"SIT",1,2,1,2,1,recLumi_i,mcDir+mcShapeSubDir+"MuSITEff/MC/probes.root",mcDir)
        StaeffB_i=ROOT.calculateDataEfficiency_v3(0,str(eosFile),".",str(run),chunk_i,LSchunks[chunk_i][0],LSchunks[chunk_i][-1],avgPileup_i,"Sta",0,2,2,2,2,recLumi_i,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)
        StaeffE_i=ROOT.calculateDataEfficiency_v3(0,str(eosFile),".",str(run),chunk_i,LSchunks[chunk_i][0],LSchunks[chunk_i][-1],avgPileup_i,"Sta",1,2,2,2,2,recLumi_i,mcDir+mcShapeSubDir+"MuStaEff/MC/probes.root",mcDir)
        Zyield_i=ROOT.calculateDataEfficiency_v3(1,str(eosFile),".",str(run),chunk_i,LSchunks[chunk_i][0],LSchunks[chunk_i][-1],avgPileup_i,"HLT",0,0,0,0,0)
        log.debug("======perMuonEff: %f, %f ,%f, %f, %f, %f",HLTeffB_i,HLTeffE_i,SITeffB_i,SITeffE_i,StaeffB_i,StaeffE_i)
        log.debug("======ZRawYield: %f",Zyield_i)

	#ZtoMuMu efficiency considering di-mu correlation from MC
        ZMCEff=ROOT.calculateZEfficiency(0,avgPileup_i,HLTeffB_i,HLTeffE_i,SITeffB_i,SITeffE_i,StaeffB_i,StaeffE_i)
        ZMCEffBB = ROOT.calculateZEfficiency(1,avgPileup_i,HLTeffB_i,HLTeffE_i,SITeffB_i,SITeffE_i,StaeffB_i,StaeffE_i)
        ZMCEffBE = ROOT.calculateZEfficiency(2,avgPileup_i,HLTeffB_i,HLTeffE_i,SITeffB_i,SITeffE_i,StaeffB_i,StaeffE_i)
        ZMCEffEE = ROOT.calculateZEfficiency(3,avgPileup_i,HLTeffB_i,HLTeffE_i,SITeffB_i,SITeffE_i,StaeffB_i,StaeffE_i)

	#ZtoMuMu efficiency purely from data
        ZBBEff=(StaeffB_i*StaeffB_i * SITeffB_i*StaeffB_i * (1-(1-HLTeffB_i)*(1-HLTeffB_i)))
        ZBEEff=(StaeffB_i*StaeffE_i * SITeffB_i*SITeffE_i * (1-(1-HLTeffB_i)*(1-HLTeffE_i)))
        ZEEEff=(StaeffE_i*StaeffE_i * SITeffE_i*StaeffE_i * (1-(1-HLTeffE_i)*(1-HLTeffE_i)))
        log.debug("======ZToMuMuEff: %f",ZMCEff)
        log.debug("======ZToMuMuEff: %f, %f ,%f, %f, %f, %f",ZMCEffBB, ZMCEffBE, ZMCEffEE, ZBBEff, ZBEEff, ZEEEff)

	#End products
        ZXSec  = Zyield_i*(1-0.01)/(ZMCEff*recLumi_i)
        ZRate  = Zyield_i*(1-0.01)/(ZMCEff*timeWindow_i*deadtime_i)
        log.debug("======ZXSec: %f",ZXSec)
        log.debug("======ZRate: %f",ZRate)

	#Variables to write in csv file
        fillarray.append(fill)
        beginTime.append(time_chunks[chunk_i][0])
        endTime.append(time_chunks[chunk_i][-1])
        Zrate.append(ZRate)
        instDel.append(delLumi_i/timeWindow_i)
        lumiDel.append(delLumi_i)
        ZyieldDel.append(Zyield_i/(ZMCEff*deadtime_i))

	#Additional variables to write in efficiency csv file
        ZyieldRec.append(Zyield_i)
        lumiRec.append(recLumi_i)
        windowarray.append(timeWindow_i)
        deadTime.append(deadtime_i)
        beginLS.append(LSchunks[chunk_i][0])
        endLS.append(LSchunks[chunk_i][-1])
	#Efficiency related
    	HLTeffB.append(HLTeffB_i)
    	HLTeffE.append(HLTeffE_i)
        SITeffB.append(SITeffB_i)
        SITeffE.append(SITeffE_i)
        StaeffB.append(StaeffB_i)
        StaeffE.append(StaeffE_i)

        ZMCeff.append(ZMCEff)
        ZMCeffBB.append(ZMCEffBB)
        ZMCeffBE.append(ZMCEffBE)
        ZMCeffEE.append(ZMCEffEE)

        ZBBeff.append(ZBBEff)
        ZBEeff.append(ZBEEff)
        ZEEeff.append(ZEEEff)


    with open('csvfile'+str(run)+'.csv','wb') as file:
        for c in range(0,nMeasurements):
                file.write(str(int(fillarray[c]))+","+str(beginTime[c])+","+str(endTime[c])+","+str(Zrate[c])+","+str(instDel[c])+","+str(lumiDel[c])+","+str(ZyieldDel[c]))
		file.write('\n')

    with open('effcsvfile'+str(run)+'.csv','wb') as file:
        for c in range(0,nMeasurements):
                file.write(str(int(fillarray[c]))+","+str(beginTime[c])+","+str(endTime[c])+","+str(Zrate[c])+","+str(instDel[c])+","+str(lumiDel[c])+","+str(ZyieldDel[c])+","+str(beginLS[c])+","+str(endLS[c])+","+str(lumiRec[c])+","+str(windowarray[c])+","+str(HLTeffB[c])+","+str(HLTeffE[c])+","+str(SITeffB[c])+","+str(SITeffE[c])+","+str(StaeffB[c])+","+str(StaeffE[c])+","+str(ZMCeff[c])+","+str(ZMCeffBB[c])+","+str(ZMCeffBE[c])+","+str(ZMCeffEE[c])+","+str(ZBBeff[c])+","+str(ZBEeff[c])+","+str(ZEEeff[c]))
                file.write('\n')
