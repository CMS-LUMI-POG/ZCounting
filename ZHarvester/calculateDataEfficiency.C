#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TTree.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TPad.h>
#include <TLatex.h>

#include "Utils/ZMMSignals.hh"
#include "Utils/ZMMBackgrounds.hh"
#include "Utils/RooGaussDoubleSidedExp.h"
#endif

// RooFit headers
#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooGaussian.h"

#include "RooPlot.h"
#include "RooGlobalFunc.h"

// constants

// Set up pile-up bounds available get MC template for fitting
const Int_t minPU = 1;
const Int_t maxPU = 60;


// Global settings
//--------------------------------------------------------------------------------------------------
Float_t massLo  = 56.;
Float_t massHi  = 116.;
UInt_t  massBin = (UInt_t)(massHi-massLo);
Float_t massWidth = 1.0;

Float_t npvLo = -0.5;
Float_t npvHi = 74.5;
UInt_t npvBin = (UInt_t)(npvHi-npvLo);

Float_t etaCutTag   = 2.4;
Float_t etaCutProbe = 2.4;
Float_t etaBound    = 0.9;

Float_t ptCutTag = 30.;
Float_t ptCutProbe = 30.;

Float_t energy = 13.0;
Float_t luminosity = 0.0;

TString outputDir="";

char lumitext[100] = "";

// precision to integrate PDF in fit instead of taking value at maximum (0 to deactivate)
// Int_t integrateBinsPrecision = 0;

// Global setters
//--------------------------------------------------------------------------------------------------
// void set_IntegrateBinsPrecision(Float_t precision_){
//     integrateBinsPrecision = precision_;

//     std::cout<<"Set IntegrateBins() precision to "<<precision_<<std::endl;
// }

void set_massRange(Float_t massLo_, Float_t massHi_, UInt_t nBins=0){
    massLo = massLo_;
    massHi = massHi_;
    if(nBins == 0){
        massBin = (UInt_t)(massHi-massLo);
    }
    else{
        massBin = nBins;
    }
    
    massWidth = (massHi-massLo)/massBin;

    std::cout<<"Set mass range to ["<<massLo_<<","<<massHi_<<"] with "<<massBin<<" bins."<<std::endl;
}

void set_npvRange(Float_t npvLo_, Float_t npvHi_){
    npvLo = npvLo_;
    npvHi = npvHi_;
    npvBin = (UInt_t)(npvHi-npvLo);

    std::cout<<"Set npv range to ["<<npvLo_<<","<<npvHi_<<"] with "<<npvBin<<" bins."<<std::endl;
}

void set_ptCut(Float_t pt_){
    ptCutTag = pt_;
    ptCutProbe = pt_;
    std::cout<<"Set pT cut of tag and probe muons to pT > "<<pt_<<" GeV"<<std::endl;
}

void set_etaCut(Float_t eta_){
    etaCutTag = eta_;
    etaCutProbe = eta_;
    std::cout<<"Set eta cut of tag and probe muons to |eta| < "<<eta_<<" "<<std::endl;
}

void set_etaBound(Float_t eta_){
    // set the boundary where barrel and endcap is defined
    etaBound = eta_;
    std::cout<<"Set eta boundary for muons to |eta| < "<<eta_<<" "<<std::endl;
}

void set_output(TString dir_){
    outputDir = dir_;
    std::cout<<"Set output directory to "<<outputDir<<std::endl;
}

void set_lumienergy(Float_t lumi_, Float_t energy_){
    energy = energy_;
    luminosity = lumi_; 

    if (static_cast<int>(energy) == energy) {
        // int value
        sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = %i TeV", luminosity, static_cast<int>(energy));
    }
    else {
        // non-integer value
        sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = %.1f TeV", luminosity, energy);
    }    
}

void set_luminosity(Float_t lumi_){
    set_lumienergy(lumi_, energy);
}

void set_energy(Float_t energy_){
    set_lumienergy(luminosity, energy_);
}

//--------------------------------------------------------------------------------------------------
// generate template for extraction of muon efficiency in barrel (B), endcap (E) or inclusive (I) region
TFile* generateGenTemplate(
    const TString mcfilename
);

//--------------------------------------------------------------------------------------------------
// generate template for extraction of muon efficiency in barrel (B), endcap (E) or inclusive (I) region
TFile* generateTemplate(
    const TString mcfilename,
    const TString effType,
    TH1D *hPV=0
);

//--------------------------------------------------------------------------------------------------
// generate template for extraction of Z yield in {BB, BE, EE, I} region
TFile* generateTemplate_ZYield(
	const TString mcfilename,
	TH1D          *hPV,
    const int     iBin
);

//--------------------------------------------------------------------------------------------------
// generate template for extraction of HLT correlation coefficient
TFile* generateTemplate_cHLT(const TString mcfilename);

//--------------------------------------------------------------------------------------------------
// extract the HLT correlation factor from the MC
double extractCorrelation_HLT(const TString mcfilename, TH1D *hPV, const TString etaRegion);
//--------------------------------------------------------------------------------------------------


std::vector<double> preFit(TH1 *failHist);

//--------------------------------------------------------------------------------------------------
// Signal Model:
// 	    1: Breit-Wigner convolved with Crystal Ball function
// 	    2: MC template convolved with Gaussian
//      3: Breit-Wigner
//      4: MC template
//      5: Breit-Wigner convolved with Gaussian
//      6: MC template convolved with Crystal Ball function

Int_t set_signal_model(
    Int_t model_type,
    CSignalModel *&model,
    RooRealVar &param_mass,
    const Bool_t pass,
    Int_t ibin=0,
    TH1D *hist=0
){

    switch(model_type) {
        case 1:
            model = new CBreitWignerConvCrystalBall(param_mass, pass, ibin);
            return 4;
        case 2:
            model = new CMCTemplateConvGaussian(param_mass, hist, pass, ibin);
            return 2;
        case 3:
            model = new CBreitWigner(param_mass, pass, ibin);
            return 0;
        case 4:
            model = new CMCTemplate(param_mass, hist, pass, ibin);
            return 0;
        case 5:
            model = new CBreitWignerConvGaussian(param_mass, pass, ibin);
            return 2;
        case 6:
            model = new CMCTemplateConvCrystalBall(param_mass, hist, pass, ibin);
            return 4;
    }
    return 0;
}

//--------------------------------------------------------------------------------------------------
// Background Model:
//      0: uniform background distribution
//      1: exponential model
//      2: quadratic model
//      3: exponential + quadratic model
//      4: Das (wide gaussian with exponential tails) model
//      5: Das + decaying exponential model
//      6: RooCMSShape
//      (7: bkg QCD + ttbar template - not working)
Int_t set_background_model(
    Int_t model_type,
    CBackgroundModel *&model,
    RooRealVar &param_mass,
    const Bool_t pass,
    Int_t ibin=0,
    TH1D *hist=0,
    std::vector<double> *params=0
){
    switch(model_type) {
        case 0:
            model = new CUniform(param_mass, pass, ibin);
            return 0;
        case 1:
            model = new CExponential(param_mass, pass, ibin);
            return 1;
        case 2:
            if(params != 0)
                model = new CQuadratic(param_mass, pass, ibin,
                    params->at(0), params->at(1), params->at(2), params->at(3), params->at(4), params->at(5));
            else
                model = new CQuadratic(param_mass, pass, 0, 0.,0.,0.,0.,0.,0.);
            return 3;
        case 3:
            if(params != 0)
                model = new CQuadPlusExp(param_mass, pass, ibin,
                    params->at(0), params->at(1), params->at(2), params->at(3), params->at(4), params->at(5));
            else
                model = new CQuadPlusExp(param_mass, pass, ibin, 0.,0.,0.,0.,0.,0.);
                return 4;
        case 4:
            model = new CDas(param_mass, pass, ibin);
            return 4;
        case 5:
            model = new CDasPlusExp(param_mass, pass, ibin);
            return 6;
        case 6:
            model = new CRooCMSShape(param_mass, pass, ibin);
            return 3;
        case 7:
            model = new CQCD(param_mass, hist, pass, ibin);
            return 1;
    }
    return 0;
}

//--------------------------------------------------------------------------------------------------
TFile* generateGenTemplate(
    const TString mcfilename
){
    const TString histfilename = outputDir+"/histTemplates_Gen.root";


    // We can use the template from the last slize if it exists
    TFile *outfile = TFile::Open(histfilename,"READ");
    if(outfile){
        cout << "Use existing template "<< endl;
        return outfile;
    }

    outfile = TFile::Open(histfilename,"RECREATE");

    cout << "Creating histogram templates... "; cout.flush();

    TFile *infile    = new TFile(mcfilename);
    TTree *eventTree = (TTree*)infile->Get("Truth");

    Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
    Double_t wgt;
    
    eventTree->SetBranchAddress("massGen",       &mass);
    eventTree->SetBranchAddress("ptGen1",        &ptTag);
    eventTree->SetBranchAddress("ptGen2",        &ptProbe);
    eventTree->SetBranchAddress("etaGen1",       &etaTag);
    eventTree->SetBranchAddress("etaGen2",       &etaProbe);
    eventTree->SetBranchAddress("eventWeight",   &wgt);

    TH1D *h_mass = new TH1D("h_mass", "", massBin, massLo, massHi);
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        eventTree->GetEntry(ientry);

        if(mass < massLo)  continue;
        if(mass > massHi)  continue;
        if(ptTag   < ptCutTag)   continue;
        if(ptProbe   < ptCutProbe)   continue;
        if(fabs(etaTag) > etaCutTag) continue;
        if(fabs(etaProbe) > etaCutProbe) continue;

        h_mass->Fill(mass, wgt);
    }

    // set negative bin entries to 0
    for(int i=1; i <= massBin; i++){
        if(h_mass->GetBinContent(i) < 0.)
            h_mass->SetBinContent(i, 0.);
    }


    
    outfile->cd();
    h_mass->Write();
    outfile->Write();

    infile->Close();
    delete infile;

    cout << "Done!" << endl;
    return outfile;
}

//--------------------------------------------------------------------------------------------------
TFile* generateTemplate(
    const TString mcfilename,
    const TString effType,
    TH1D          *hPV  // primary vertex distribution in data, must be normalized to integral 1
){
    const TString histfilename = outputDir+"/histTemplates_"+effType+".root";

    TFile *outfile = 0;

    if(hPV==0){
        // if no PU reweighting is done, we can use the template from the last slize if it exists

        outfile = TFile::Open(histfilename,"READ");
        if(outfile){
            cout << "Use existing template "<< endl;
            return outfile;
        }
    }

    outfile = TFile::Open(histfilename,"RECREATE");

    cout << "Creating histogram templates... "; cout.flush();

    TFile *infile    = new TFile(mcfilename);
    TTree *eventTree = (TTree*)infile->Get(effType);
    TH1D *hPV_scalefactor = 0;

    if(hPV){
        std::cout<<"PV reweighting with <PV> = "<<hPV->GetMean()<<std::endl;

        TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");  // get primary vertex distribution in MC

        hPV_scalefactor = (TH1D*)hPV->Clone("hPV_scalefactor");
        hPV_scalefactor->Divide(hPVtemplate);
    }

    Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
    Double_t wgt;
    Int_t npv;
    Int_t pass;
    Bool_t match1, match2;
    
    eventTree->SetBranchAddress("mass",           &mass);
    eventTree->SetBranchAddress("ptTag",          &ptTag);
    eventTree->SetBranchAddress("ptProbe",        &ptProbe);
    eventTree->SetBranchAddress("etaTag",         &etaTag);
    eventTree->SetBranchAddress("etaProbe",       &etaProbe);
    eventTree->SetBranchAddress("nPV",            &npv);
    eventTree->SetBranchAddress("pass",           &pass);
    eventTree->SetBranchAddress("match1",         &match1);
    eventTree->SetBranchAddress("match2",         &match2);
    eventTree->SetBranchAddress("eventWeight",    &wgt);

    TH1D *h_mass_pass_B = new TH1D("h_mass_pass_B", "", massBin, massLo, massHi);
    TH1D *h_mass_fail_B = new TH1D("h_mass_fail_B", "", massBin, massLo, massHi);
    TH1D *h_mass_pass_E = new TH1D("h_mass_pass_E", "", massBin, massLo, massHi);
    TH1D *h_mass_fail_E = new TH1D("h_mass_fail_E", "", massBin, massLo, massHi);
    TH1D *h_mass_pass_I = new TH1D("h_mass_pass_I", "", massBin, massLo, massHi);
    TH1D *h_mass_fail_I = new TH1D("h_mass_fail_I", "", massBin, massLo, massHi);
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        eventTree->GetEntry(ientry);

        if(!match1 || !match2) continue;
        if(mass < massLo)  continue;
        if(mass > massHi)  continue;
        if(ptTag   < ptCutTag)   continue;
        if(ptProbe   < ptCutProbe)   continue;
        if(fabs(etaTag) > etaCutTag) continue;
        if(fabs(etaProbe) > etaCutProbe) continue;

        if(hPV_scalefactor)
            wgt *= hPV_scalefactor->GetBinContent(hPV_scalefactor->FindBin(std::max(1,npv)));

        if(fabs(etaProbe) < etaBound){
            if(pass) h_mass_pass_B->Fill(mass, wgt);
            else     h_mass_fail_B->Fill(mass, wgt);
        }else{
            if(pass) h_mass_pass_E->Fill(mass, wgt);
            else     h_mass_fail_E->Fill(mass, wgt);
        }
    }

    // set negative bin entries to 0
    for(int i=1; i <= massBin; i++){
        if(h_mass_pass_B->GetBinContent(i) < 0.)
            h_mass_pass_B->SetBinContent(i, 0.);

        if(h_mass_fail_B->GetBinContent(i) < 0.)
            h_mass_fail_B->SetBinContent(i, 0.);

        if(h_mass_pass_E->GetBinContent(i) < 0.)
            h_mass_pass_E->SetBinContent(i, 0.);

        if(h_mass_fail_E->GetBinContent(i) < 0.)
            h_mass_fail_E->SetBinContent(i, 0.);

    }
    h_mass_pass_I->Add(h_mass_pass_B);
    h_mass_fail_I->Add(h_mass_fail_B);
    h_mass_pass_I->Add(h_mass_pass_E);
    h_mass_fail_I->Add(h_mass_fail_E);
    
    // add passing histograms to failing ones to increase statistics
    h_mass_fail_B->Add(h_mass_pass_B);
    h_mass_fail_E->Add(h_mass_pass_E);
    h_mass_fail_I->Add(h_mass_pass_I);

    outfile->cd();
    h_mass_pass_B->Write();
    h_mass_fail_B->Write();
    h_mass_pass_E->Write();
    h_mass_fail_E->Write();
    h_mass_pass_I->Write();
    h_mass_fail_I->Write();
    outfile->Write();

    infile->Close();
    delete infile;

    cout << "Done!" << endl;
    return outfile;
}

//--------------------------------------------------------------------------------------------------
TFile* generateTemplate_ZYield(
    const TString mcfilename,
    TH1D          *hPV,
    const int     iBin
){
    //const TString histfilename = hPV == 0 ? outputDir+"/../histTemplates_HLT.root" : outputDir+"/histTemplates_HLT_"+std::to_string(iBin)+".root";
    const TString histfilename = outputDir+"/histTemplates_HLT.root";


    TFile *outfile = 0;

    if(hPV==0){
        // if no PU reweighting is done, we can use the template from the last slize if it exists

        outfile = TFile::Open(histfilename,"READ");
        if(outfile){
            cout << "Use existing template "<< endl;
            return outfile;
        }
    }

    outfile = TFile::Open(histfilename,"RECREATE");

    cout << "Creating histogram templates... "; cout.flush();

    TFile *infile    = new TFile(mcfilename);
    TTree *eventTree = (TTree*)infile->Get("HLT");    
    TH1D *hPV_scalefactor = 0;

    if(hPV){
        std::cout<<"PV reweighting with <PV> = "<<hPV->GetMean()<<std::endl;

        TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

        hPV_scalefactor = (TH1D*)hPV->Clone("hPV_scalefactor");
        hPV_scalefactor->Divide(hPVtemplate);
    }

    Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
    Double_t wgt;
    Int_t npv;
    Int_t pass;
    Bool_t match1, match2;

    eventTree->SetBranchAddress("mass",           &mass);
    eventTree->SetBranchAddress("ptTag",          &ptTag);
    eventTree->SetBranchAddress("ptProbe",        &ptProbe);
    eventTree->SetBranchAddress("etaTag",         &etaTag);
    eventTree->SetBranchAddress("etaProbe",       &etaProbe);
    eventTree->SetBranchAddress("nPV",            &npv);
    eventTree->SetBranchAddress("pass",           &pass);
    eventTree->SetBranchAddress("match1",       &match1);
    eventTree->SetBranchAddress("match2",     &match2);    
    eventTree->SetBranchAddress("eventWeight",    &wgt);
    
    TH1D *h_mass_0hlt_BB = new TH1D("h_mass_0hlt_BB", "", massBin, massLo, massHi);
    TH1D *h_mass_0hlt_BE = new TH1D("h_mass_0hlt_BE", "", massBin, massLo, massHi);
    TH1D *h_mass_0hlt_EE = new TH1D("h_mass_0hlt_EE", "", massBin, massLo, massHi);
    TH1D *h_mass_0hlt_I  = new TH1D("h_mass_0hlt_I",  "", massBin, massLo, massHi);
    TH1D *h_mass_1hlt_BB = new TH1D("h_mass_1hlt_BB", "", massBin, massLo, massHi);
    TH1D *h_mass_1hlt_BE = new TH1D("h_mass_1hlt_BE", "", massBin, massLo, massHi);
    TH1D *h_mass_1hlt_EE = new TH1D("h_mass_1hlt_EE", "", massBin, massLo, massHi);
    TH1D *h_mass_1hlt_I  = new TH1D("h_mass_1hlt_I",  "", massBin, massLo, massHi);
    TH1D *h_mass_2hlt_BB = new TH1D("h_mass_2hlt_BB", "", massBin, massLo, massHi);
    TH1D *h_mass_2hlt_BE = new TH1D("h_mass_2hlt_BE", "", massBin, massLo, massHi);
    TH1D *h_mass_2hlt_EE = new TH1D("h_mass_2hlt_EE", "", massBin, massLo, massHi);
    TH1D *h_mass_2hlt_I  = new TH1D("h_mass_2hlt_I",  "", massBin, massLo, massHi);

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        eventTree->GetEntry(ientry);
        
        if(!match1 || !match2) continue;
        if(mass < massLo)  continue;
        if(mass > massHi)  continue;
        if(ptTag   < ptCutTag)   continue;
        if(ptProbe   < ptCutProbe)   continue;
        if(fabs(etaTag) > etaCutTag) continue;
        if(fabs(etaProbe) > etaCutProbe) continue;

        if(hPV_scalefactor)
            wgt *= hPV_scalefactor->GetBinContent(hPV_scalefactor->FindBin(std::max(1,npv)));

        if(fabs(etaProbe) < etaBound && fabs(etaTag) < etaBound){
            if(pass==2)         h_mass_2hlt_BB->Fill(mass, wgt);
            else if(pass==1)    h_mass_1hlt_BB->Fill(mass, wgt);
            else                h_mass_0hlt_BB->Fill(mass, wgt);
        }
        else if(fabs(etaProbe) >= etaBound && fabs(etaTag) >= etaBound){
            if(pass==2)         h_mass_2hlt_EE->Fill(mass, wgt);
            else if(pass==1)    h_mass_1hlt_EE->Fill(mass, wgt);
            else                h_mass_0hlt_EE->Fill(mass, wgt);
        }
        else{
            if(pass==2)         h_mass_2hlt_BE->Fill(mass, wgt);
            else if(pass==1)    h_mass_1hlt_BE->Fill(mass, wgt);
            else                h_mass_0hlt_BE->Fill(mass, wgt);
        }
    }

    // set negative bin entries to 0
    for(int i=1; i <= massBin+1; i++){
        if(h_mass_2hlt_BB->GetBinContent(i) < 0.)
            h_mass_2hlt_BB->SetBinContent(i, 0.);

        if(h_mass_2hlt_BE->GetBinContent(i) < 0.)
            h_mass_2hlt_BE->SetBinContent(i, 0.);

        if(h_mass_2hlt_EE->GetBinContent(i) < 0.)
            h_mass_2hlt_EE->SetBinContent(i, 0.);

        if(h_mass_1hlt_BB->GetBinContent(i) < 0.)
            h_mass_1hlt_BB->SetBinContent(i, 0.);

        if(h_mass_1hlt_BE->GetBinContent(i) < 0.)
            h_mass_1hlt_BE->SetBinContent(i, 0.);

        if(h_mass_1hlt_EE->GetBinContent(i) < 0.)
            h_mass_1hlt_EE->SetBinContent(i, 0.);

        if(h_mass_0hlt_BB->GetBinContent(i) < 0.)
            h_mass_0hlt_BB->SetBinContent(i, 0.);

        if(h_mass_0hlt_BE->GetBinContent(i) < 0.)
            h_mass_0hlt_BE->SetBinContent(i, 0.);

        if(h_mass_0hlt_EE->GetBinContent(i) < 0.)
            h_mass_0hlt_EE->SetBinContent(i, 0.);
    }
    // inclusive
    h_mass_0hlt_I->Add(h_mass_0hlt_BB);
    h_mass_0hlt_I->Add(h_mass_0hlt_BE);
    h_mass_0hlt_I->Add(h_mass_0hlt_EE);
    h_mass_1hlt_I->Add(h_mass_1hlt_BB);
    h_mass_1hlt_I->Add(h_mass_1hlt_BE);
    h_mass_1hlt_I->Add(h_mass_1hlt_EE);
    h_mass_2hlt_I->Add(h_mass_2hlt_BB);
    h_mass_2hlt_I->Add(h_mass_2hlt_BE);
    h_mass_2hlt_I->Add(h_mass_2hlt_EE);
    
    outfile->cd();
    
    h_mass_2hlt_BB->Write();
    h_mass_2hlt_BE->Write();
    h_mass_2hlt_EE->Write();
    h_mass_2hlt_I->Write();
    h_mass_1hlt_BB->Write();
    h_mass_1hlt_BE->Write();
    h_mass_1hlt_EE->Write();
    h_mass_1hlt_I->Write();
    h_mass_0hlt_BB->Write();
    h_mass_0hlt_BE->Write();
    h_mass_0hlt_EE->Write();
    h_mass_0hlt_I->Write();
    

    outfile->Write();

    infile->Close();
    delete infile;

    cout << "Done!" << endl;
    return outfile;
}

// generate histogram templates binned as function of number of primary vertices, for the extraction of the correlation coefficient
TFile* generateTemplate_cHLT(
    const TString mcfilename
){
    const TString histfilename = outputDir+"/histTemplates_cHLT.root";
    TFile *outfile = TFile::Open(histfilename,"CREATE");
    if(!outfile){
        cout << "Use existing template "<< endl;
        outfile = TFile::Open(histfilename,"READ");
        return outfile;
    }
    cout << "Creating histogram templates... "; cout.flush();

    TFile *infile    = new TFile(mcfilename);
    TTree *eventTree = (TTree*)infile->Get("HLT");
    TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

    Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
    Double_t wgt;
    Int_t npv;
    Int_t pass;
    Bool_t match1, match2;
    
    eventTree->SetBranchAddress("mass",           &mass);
    eventTree->SetBranchAddress("ptTag",          &ptTag);
    eventTree->SetBranchAddress("ptProbe",        &ptProbe);
    eventTree->SetBranchAddress("etaTag",         &etaTag);
    eventTree->SetBranchAddress("etaProbe",       &etaProbe);
    eventTree->SetBranchAddress("nPV",            &npv);
    eventTree->SetBranchAddress("pass",           &pass);
    eventTree->SetBranchAddress("match1",       &match1);
    eventTree->SetBranchAddress("match2",     &match2);    
    eventTree->SetBranchAddress("eventWeight",    &wgt);
    
    TH1D *h_npv_0hlt_BB = new TH1D("h_npv_0hlt_BB", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_0hlt_BE = new TH1D("h_npv_0hlt_BE", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_0hlt_EE = new TH1D("h_npv_0hlt_EE", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_0hlt_I  = new TH1D("h_npv_0hlt_I",  "", npvBin, npvLo, npvHi);
    TH1D *h_npv_1hlt_BB = new TH1D("h_npv_1hlt_BB", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_1hlt_BE = new TH1D("h_npv_1hlt_BE", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_1hlt_EE = new TH1D("h_npv_1hlt_EE", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_1hlt_I  = new TH1D("h_npv_1hlt_I",  "", npvBin, npvLo, npvHi);
    TH1D *h_npv_2hlt_BB = new TH1D("h_npv_2hlt_BB", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_2hlt_BE = new TH1D("h_npv_2hlt_BE", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_2hlt_EE = new TH1D("h_npv_2hlt_EE", "", npvBin, npvLo, npvHi);
    TH1D *h_npv_2hlt_I  = new TH1D("h_npv_2hlt_I",  "", npvBin, npvLo, npvHi);

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        eventTree->GetEntry(ientry);
        
        if(!match1 || !match2) continue;
        if(mass < massLo)  continue;
        if(mass > massHi)  continue;
        if(ptTag   < ptCutTag)   continue;
        if(ptProbe   < ptCutProbe)   continue;
        if(fabs(etaTag) > etaCutTag) continue;
        if(fabs(etaProbe) > etaCutProbe) continue;

        if(fabs(etaProbe) < etaBound && fabs(etaTag) < etaBound){
            if(pass==2)         h_npv_2hlt_BB->Fill(npv, wgt);
            else if(pass==1)    h_npv_1hlt_BB->Fill(npv, wgt);
            else                h_npv_0hlt_BB->Fill(npv, wgt);
        }
        else if(fabs(etaProbe) >= etaBound && fabs(etaTag) >= etaBound){
            if(pass==2)         h_npv_2hlt_EE->Fill(npv, wgt);
            else if(pass==1)    h_npv_1hlt_EE->Fill(npv, wgt);
            else                h_npv_0hlt_EE->Fill(npv, wgt);
        }
        else{
            if(pass==2)         h_npv_2hlt_BE->Fill(npv, wgt);
            else if(pass==1)    h_npv_1hlt_BE->Fill(npv, wgt);
            else                h_npv_0hlt_BE->Fill(npv, wgt);
        }
    }

    // set negative bin entries to 0
    for(int i=1; i <= npvBin+1; i++){
        if(h_npv_2hlt_BB->GetBinContent(i) < 0.)
            h_npv_2hlt_BB->SetBinContent(i, 0.);

        if(h_npv_2hlt_BE->GetBinContent(i) < 0.)
            h_npv_2hlt_BE->SetBinContent(i, 0.);

        if(h_npv_2hlt_EE->GetBinContent(i) < 0.)
            h_npv_2hlt_EE->SetBinContent(i, 0.);

        if(h_npv_1hlt_BB->GetBinContent(i) < 0.)
            h_npv_1hlt_BB->SetBinContent(i, 0.);

        if(h_npv_1hlt_BE->GetBinContent(i) < 0.)
            h_npv_1hlt_BE->SetBinContent(i, 0.);

        if(h_npv_1hlt_EE->GetBinContent(i) < 0.)
            h_npv_1hlt_EE->SetBinContent(i, 0.);

        if(h_npv_0hlt_BB->GetBinContent(i) < 0.)
            h_npv_0hlt_BB->SetBinContent(i, 0.);

        if(h_npv_0hlt_BE->GetBinContent(i) < 0.)
            h_npv_0hlt_BE->SetBinContent(i, 0.);

        if(h_npv_0hlt_EE->GetBinContent(i) < 0.)
            h_npv_0hlt_EE->SetBinContent(i, 0.);
    }
    // inclusive
    h_npv_0hlt_I->Add(h_npv_0hlt_BB);
    h_npv_0hlt_I->Add(h_npv_0hlt_BE);
    h_npv_0hlt_I->Add(h_npv_0hlt_EE);
    h_npv_1hlt_I->Add(h_npv_1hlt_BB);
    h_npv_1hlt_I->Add(h_npv_1hlt_BE);
    h_npv_1hlt_I->Add(h_npv_1hlt_EE);
    h_npv_2hlt_I->Add(h_npv_2hlt_BB);
    h_npv_2hlt_I->Add(h_npv_2hlt_BE);
    h_npv_2hlt_I->Add(h_npv_2hlt_EE);
    
    outfile->cd();
    
    h_npv_2hlt_BB->Write();
    h_npv_2hlt_BE->Write();
    h_npv_2hlt_EE->Write();
    h_npv_2hlt_I->Write();
    h_npv_1hlt_BB->Write();
    h_npv_1hlt_BE->Write();
    h_npv_1hlt_EE->Write();
    h_npv_1hlt_I->Write();
    h_npv_0hlt_BB->Write();
    h_npv_0hlt_BE->Write();
    h_npv_0hlt_EE->Write();
    h_npv_0hlt_I->Write();
    hPVtemplate->Write();
    

    outfile->Write();

    infile->Close();
    delete infile;

    cout << "Done!" << endl;
    return outfile;
}


//--------------------------------------------------------------------------------------------------
// template<typename T>
Double_t make_plot(
    const Int_t   nfl,
    const RooRealVar &param_mass,
    RooAbsData *data,
    CSignalModel *sigModel,         // only needed for label
    CBackgroundModel *bkgModel,     // only needed for label
    const Int_t nEntries,
    RooAbsPdf* modelPdf,
    RooRealVar* Nsig,
    RooRealVar* Nbkg,
    const Int_t iBin,
    RooRealVar *eff = 0,
    RooRealVar *corr = 0,
    const TString effType="yield",
    const TString etaRegion="I",
    const Int_t passRegion=1,
    const Bool_t plot_ratio=kTRUE,
    const Bool_t logscale=kFALSE
){
    char pname[50];
    TString suffix = "";    
    char ctitle[100];
    char binlabelx[100];
    char binlabely[100];
    const TString xlabel = "tag-and-probe mass [GeV]";
    char ylabel[50];
    char yield[50];
    char nsigstr[100];
    char nbkgstr[100];
    char chi2str[100];
    char effstr[100] = "";
    char cstr[100] = "";

    TString passstr = passRegion ? "Pass" : "Fail";
    TString sigmodstr = "signal"+passstr;
    TString bkgmodstr = "background"+passstr;

    // if(effType=="Sta"){
    //     sigmodstr+= "_1";
    //     bkgmodstr+= "_1";
    // }
    // else if(effType=="Trk"){
    //     sigmodstr+= "_2";
    //     bkgmodstr+= "_2";        
    // }
    // else 
    if(effType == "yield"){
        // passstr = "HLT"+correlationFit;
        sigmodstr += "_"+std::to_string(passRegion);
        bkgmodstr += "_"+std::to_string(passRegion);
    }
    else{
        sigmodstr+= "_0";
        bkgmodstr+= "_0";
    }
    
    if(etaRegion=="B" || etaRegion == "BB") sprintf(binlabelx, "0.0 < |#eta| < %.1f",etaBound);
    else if(etaRegion=="E" || etaRegion == "EE") sprintf(binlabelx, "%.1f < |#eta| < %.1f",etaBound, etaCutTag);
    else sprintf(binlabelx, "|#eta| < %.1f", etaCutTag);


    sprintf(pname,"%s_%s_%i_%d", effType.Data(), etaRegion.Data(), passRegion, iBin);
    sprintf(ctitle,"%s %i (%d) ", effType.Data(), passRegion, iBin);

    
    if(eff != 0){
        if(effType == "yield"){
            sprintf(effstr,"#varepsilon^{#mu}_{HLT} = %.4f_{ -%.4f}^{ +%.4f}",
                eff->getVal(),std::abs(eff->getErrorLo()),eff->getErrorHi());
        }
        else {
            sprintf(effstr,"#varepsilon^{#mu}_{%s} = %.4f_{ -%.4f}^{ +%.4f}", 
                effType.Data(), eff->getVal(), std::abs(eff->getErrorLo()), eff->getErrorHi());
        }
    }
    if(corr != 0){
        if(corr->getErrorLo() == 0 || corr->getErrorHi() == 0)
            sprintf(cstr,"c = %.4f", corr->getVal());
        else
            sprintf(cstr,"c = %.4f_{ -%.4f}^{ +%.4f}",corr->getVal(),std::abs(corr->getErrorLo()),corr->getErrorHi());
    }
    
    const double margin_left = 0.15;
    const double margin_right = 0.03;
    const double margin_bottom_p1 = plot_ratio ? 0.025 : 0.15;
    const double margin_top_p1 = 0.01;

    TCanvas *canvas = 0;
    TPad *pad1 = 0;
    if(plot_ratio){
        canvas = new TCanvas(pname, ctitle, 800,800);
        pad1 = new TPad("pad1", "pad1", 0., 0.35, 1, 1.0);
    }
    else{
        canvas = new TCanvas(pname, ctitle, 800, 540);
        pad1 = new TPad("pad1", "pad1", 0., 0.0, 1, 1.0);       
    }

    canvas->SetTicks();
    pad1->SetLeftMargin(margin_left);
    pad1->SetRightMargin(margin_right);
    pad1->SetTopMargin(margin_top_p1);
    pad1->SetBottomMargin(margin_bottom_p1);
    pad1->SetTickx();
    pad1->SetTicky();
    pad1->Draw();
    pad1->cd();

    const double textsize1 = 28./(pad1->GetWh()*pad1->GetAbsHNDC());

    TLegend *legend = new TLegend(margin_left+0.04, 0.47, margin_left+0.36, 0.79);
    legend->SetNColumns(1);
    legend->SetBorderSize(0);

    TLegendEntry *entry1 = new TLegendEntry();
    entry1->SetTextSize(textsize1);
    entry1->SetLabel("Data");
    entry1->SetOption("PE");
    entry1->SetMarkerStyle(kFullCircle);
    entry1->SetMarkerColor(kBlack);
    entry1->SetTextFont(42);
    entry1->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry1);

    TLegendEntry *entry4 = new TLegendEntry();
    entry4->SetTextSize(textsize1);
    entry4->SetLabel("Full model");
    entry4->SetOption("l");
    entry4->SetLineColor(2);
    entry4->SetLineWidth(2);
    entry4->SetTextFont(42);
    entry4->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry4);

    TLegendEntry *entry2 = new TLegendEntry();
    entry2->SetTextSize(textsize1);
    entry2->SetLabel(bkgModel->model->GetTitle());
    entry2->SetOption("l");
    entry2->SetLineColor(8);
    entry2->SetLineWidth(2);
    entry2->SetTextFont(42);
    entry2->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry2);

    TLegendEntry *entry3 = new TLegendEntry();
    entry3->SetTextSize(textsize1);
    entry3->SetLabel(sigModel->model->GetTitle());
    entry3->SetOption("l");
    entry3->SetLineColor(9);
    entry3->SetLineWidth(2);
    entry3->SetTextFont(42);
    entry3->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry3);

    RooPlot *mframe = param_mass.frame(RooFit::Bins(massBin));
    data->plotOn(mframe, RooFit::MarkerStyle(kFullCircle), RooFit::MarkerSize(0.8), RooFit::DrawOption("ZP"));
    modelPdf->plotOn(mframe, RooFit::Components(bkgmodstr), RooFit::LineColor(8));
    modelPdf->plotOn(mframe, RooFit::Components(sigmodstr), RooFit::LineColor(9));
    modelPdf->plotOn(mframe, RooFit::LineColor(kRed));

    // Construct a histogram with the pulls of the data w.r.t the curve
    RooHist *hpull = mframe->pullHist();

    double a = Nsig->getVal(), aErr = Nsig->getError();//*fitResult);
    double b = Nbkg->getVal(), bErr = Nbkg->getError();//*fitResult);
    double resChi2 = mframe->chiSquare(nfl);

    sprintf(binlabely, "p_{T} > %i GeV",(Int_t)ptCutTag);
    if(massWidth != 1)
        sprintf(ylabel,"Events / %.1f GeV", massWidth);
    else
        sprintf(ylabel,"Events / 1 GeV");
    sprintf(yield,"%u Events",nEntries);
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",a,aErr);
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",resChi2);
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",b,bErr);

    mframe->SetTitle("");

    mframe->GetXaxis()->SetTitleSize(0.);
    mframe->GetXaxis()->SetLabelSize(0.);

    mframe->GetYaxis()->SetTitle(ylabel);
    mframe->GetYaxis()->SetLabelSize(textsize1);
    mframe->GetYaxis()->SetTitleFont(42);
    mframe->GetYaxis()->SetTitleSize(textsize1*1.2);
    mframe->GetYaxis()->SetTitleOffset(1.1);
    
    if(!plot_ratio){
        mframe->GetXaxis()->SetTitle(xlabel);
        mframe->GetXaxis()->SetTitleSize(textsize1*1.2);
        mframe->GetXaxis()->SetTitleOffset(1.);
        mframe->GetXaxis()->SetLabelSize(textsize1);
        mframe->GetXaxis()->SetTitleFont(42);

        suffix += "_noRatio";
    }
    if(logscale){
        pad1->SetLogy(1);
        mframe->SetMinimum(1.);
        mframe->SetMaximum(mframe->GetMaximum()*5);
        
        suffix += "_logscale";
    }
    else{
        pad1->SetLogy(0);
        mframe->SetMinimum(0.);

        if(effType == "Trk"){
            mframe->SetMaximum(mframe->GetMaximum()*2);
        }
        else{
            mframe->SetMaximum(mframe->GetMaximum()*1.2);
        }
    }

    mframe->Draw();

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(11);
    latex->SetTextFont(62);
    latex->SetTextSize(textsize1*1.2);
    latex->DrawLatex(margin_left+0.04, 0.91, "CMS");
    latex->SetTextFont(52);
    latex->SetTextSize(textsize1);
    latex->DrawLatex(margin_left+0.14, 0.91, "Work in progress");
    latex->SetTextFont(42);
    latex->DrawLatex(margin_left+0.04, 0.84, lumitext);

    latex->DrawLatex(margin_left+0.6, 0.41, ctitle);
    latex->DrawLatex(margin_left+0.6, 0.35, yield);
    latex->DrawLatex(margin_left+0.5, 0.76, nsigstr);
    latex->DrawLatex(margin_left+0.5, 0.69, nbkgstr);
    latex->DrawLatex(margin_left+0.5, 0.62, chi2str);
    latex->DrawLatex(margin_left+0.5, 0.55, effstr);
    latex->DrawLatex(margin_left+0.5, 0.48, cstr);

    latex->SetTextAlign(31);
    latex->DrawLatex(1-margin_right-0.04, 0.91, binlabelx);
    latex->DrawLatex(1-margin_right-0.04, 0.84, binlabely);

    legend->Draw("same");
    
    if(plot_ratio){

        canvas->cd();
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
        pad2->SetLeftMargin(margin_left);
        pad2->SetRightMargin(margin_right);
        pad2->SetTopMargin(0.025);
        pad2->SetBottomMargin(0.3);
        pad2->SetTickx();
        pad2->SetTicky();
        pad2->Draw("ALPF");
        pad2->cd();

        const double textsize2 = 28./(pad2->GetWh()*pad2->GetAbsHNDC());

        RooPlot *rframe = param_mass.frame(RooFit::Bins(massBin));

        rframe->SetTitle("");
        rframe->addPlotable(hpull, "PX");

        rframe->GetYaxis()->SetTitle("Pulls");
        rframe->GetYaxis()->SetTitleSize(textsize2*1.2);
        rframe->GetYaxis()->SetTitleOffset(0.6);
        rframe->GetYaxis()->SetLabelSize(textsize2);
        rframe->GetYaxis()->SetTitleFont(42);
        rframe->GetYaxis()->CenterTitle(true);
        rframe->GetYaxis()->SetNdivisions(405);

        rframe->GetXaxis()->SetTitle(xlabel);
        rframe->GetXaxis()->SetTitleSize(textsize2*1.2);
        rframe->GetXaxis()->SetTitleOffset(1.);
        rframe->GetXaxis()->SetLabelSize(textsize2);
        rframe->GetXaxis()->SetTitleFont(42);
        rframe->Draw();

        TLine *line0 = new TLine(massLo, 0., massHi, 0.);
        line0->SetLineStyle(1);
        line0->Draw("same");
    
    }

    canvas->SaveAs(outputDir+"/"+pname+suffix+".png");
    // canvas->SaveAs(outputDir+"/"+pname+suffix+".eps");
    canvas->Close();
    
    return resChi2;
}

//--------------------------------------------------------------------------------------------------
// perform fit on one histogram to extract number of Z 
void getZyield(
    TH1D*          h_yield,         // Histogram with signal and background contribution
    const Int_t    iBin,            // Label of measurement in currect run
    const TString  effType="HLT",   // {Trk, Sta, Glo, Sel, ID, HLT}
    const TString  etaRegion="",    // {B, E, BB, BE, EE or I}
    const Int_t    sigMod=1,
    const Int_t    bkgMod=6,
    const Int_t    passRegion=1,
    const TString  mcfilename="",   
    TH1D*          hPV=0
){
    std::cout<<">>> Do fit in "<< etaRegion 
        <<" for "<< effType <<" efficiency"
        <<" in pass="<< passRegion <<" category "
        <<" for extraction of number of Z"<<std::endl;

    RooRealVar m("m", "mass", massLo, massHi);
    m.setBins(massBin);
    // This is only needed for the convolution of template with an analytic function 
    // set to a high number for good precision
    m.setBins(10000,"cache");  
    // set sideband range for background fit
    m.setRange("rangeLow", massLo, 76);
    m.setRange("rangeHigh", 106, massHi);
    m.setRange("rangeCenter", 81, 101);

    CSignalModel     *sigModel = 0;
    CBackgroundModel *bkgModel = 0;

    Int_t nfl=0;

    TFile *histfile = 0;
    TH1D *h=0;
    if(sigMod%2 == 0) {
        if(sigMod > 10){
            histfile = generateGenTemplate(mcfilename);
            h = (TH1D*)histfile->Get("h_mass");
        }
        else if(effType == "HLT"){
            histfile = generateTemplate_ZYield(mcfilename, hPV, iBin);
            if(passRegion==2)
                h = (TH1D*)histfile->Get("h_mass_2hlt_"+etaRegion);
            else if(passRegion==1)
                h = (TH1D*)histfile->Get("h_mass_1hlt_"+etaRegion);            
            else
                h = (TH1D*)histfile->Get("h_mass_0hlt_"+etaRegion);     
        }
        else{
            histfile = generateTemplate(mcfilename, effType, hPV);
            assert(histfile);
            if(passRegion)
                h = (TH1D*)histfile->Get(Form("h_mass_pass_%s", etaRegion.Data()));
            else
                h = (TH1D*)histfile->Get(Form("h_mass_fail_%s", etaRegion.Data()));
        }
        assert(h);
    }

    nfl += set_signal_model(sigMod%10, sigModel, m, passRegion, 0, h);
    nfl += set_background_model(bkgMod, bkgModel, m, passRegion, 0);
    
    // ---> Add constraint terms
    RooArgSet *constraints = new RooArgSet("constraints");
    RooGaussian *constraint_sig_mean = 0;
    RooGaussian *constraint_sig_sigma = 0;
    RooGaussian *constraint_sig_alpha = 0;
    if(sigMod%10 == 1 || sigMod%10 == 2 || sigMod%10 == 5 || sigMod%10 == 6){

        if(effType == "Trk" || effType == "Glo"){
            std::cout<<">>> set constraints on signal model parameters "<< std::endl; 
            // in case of probes with outer track parameters the fit can be unstable, in particular for failing probes,
            // in these cases, constraints are used to stabalize the fit
            // Construct Gaussian constraints p.d.f on signal model
            constraint_sig_mean = new RooGaussian("constraint_sig_mean","constraint signal mean", *sigModel->mean, RooFit::RooConst(0.), RooFit::RooConst(1.0));
            constraint_sig_sigma = new RooGaussian("constraint_sig_sigma","constraint signal sigma", *sigModel->sigma, RooFit::RooConst(2.5), RooFit::RooConst(1.5));
            
            constraints->add(*constraint_sig_mean);
            constraints->add(*constraint_sig_sigma);

            // Constraint on CB parameter
            if(sigMod%10 == 1 || sigMod%10 == 5){
                constraint_sig_alpha = new RooGaussian("constraint_sig_alpha","constraint signal alpha", *sigModel->alpha, RooFit::RooConst(5.), RooFit::RooConst(3.0));
                constraints->add(*constraint_sig_alpha);
            }
        }
    }

    RooGaussian *constraint_alpha = 0;
    RooGaussian *constraint_beta = 0;
    RooGaussian *constraint_gamma = 0;
    if(bkgMod == 6){
        std::cout<<">>> set constraints on signal model parameters "<< std::endl; 
        // Contraints are used for the RooCMSShape 

        constraint_alpha = new RooGaussian("constraint_alpha","constraint alpha", *bkgModel->alpha, RooFit::RooConst(50.), RooFit::RooConst(20.0));
        constraint_beta = new RooGaussian("constraint_beta","constraint beta", *bkgModel->beta, RooFit::RooConst(0.05), RooFit::RooConst(0.05));
        constraint_gamma = new RooGaussian("constraint_gamma","constraint gamma", *bkgModel->gamma, RooFit::RooConst(0.05), RooFit::RooConst(0.05));

        constraints->add(*constraint_alpha);
        constraints->add(*constraint_beta);
        constraints->add(*constraint_gamma);
    }
    // <---

    RooDataHist *data = new RooDataHist("ZReco","ZReco",RooArgList(m), h_yield);

    const Double_t NsigMax = h_yield->Integral();
    Double_t NsigInit = 0;  // initial signal contribution
    Double_t NbkgInit = 0;  // initial background contribution
    if(effType == "HLT"){
        NsigInit = passRegion ? NsigMax*0.995 : NsigMax*0.95;
        NbkgInit = passRegion ? 0.005 : NsigMax*0.05;
    }
    else if(effType == "Sel" || effType == "ID"){
        NsigInit = passRegion ? NsigMax*0.99 : NsigMax*0.90;
        NbkgInit = passRegion ? 0.01 : NsigMax*0.1;    
    }
    else if(effType == "Glo"){
        NsigInit = passRegion ? NsigMax*0.99 : NsigMax*0.2;
        NbkgInit = passRegion ? 0.01 : NsigMax*0.8;         
    }
    else if(effType == "Sta"){
        NsigInit = passRegion ? NsigMax*0.99 : NsigMax*0.1;
        NbkgInit = passRegion ? 0.01 : NsigMax*0.9;               
    }
    else if(effType == "Trk"){
        NsigInit = passRegion ? NsigMax*0.99 : NsigMax*0.3;
        NbkgInit = passRegion ? 0.01 : NsigMax*0.7;        

        m.setRange("rangeLow", massLo, 76);
        m.setRange("rangeHigh", 104, massHi);
        m.setRange("rangeCenter", 76, 104);
    }
    
    RooRealVar Nsig("Nsig", "Number of signal candidates", NsigInit, -0.25*NsigMax, 1.5*NsigMax);
    RooRealVar Nbkg("Nbkg", "Number of background candidates", NbkgInit, -0.25*NsigMax, 1.5*NsigMax);
    RooAddPdf modelPdf("totalPdf","Z sig+bkg",
        RooArgList(*(sigModel->model),*(bkgModel->model)),RooArgList(Nsig,Nbkg));

    // RooFormulaVar purity("purity","Nsig/(Nsig+Nbkg)",RooArgList(Nsig,Nbkg));

    TFile *fFit = new TFile(Form(
        "%s/workspace_yield_%s_%s_%i_%i.root",outputDir.Data(), effType.Data(), etaRegion.Data(), passRegion, iBin),
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    // w->import(purity);
    w->import(*data);
    w->import(modelPdf);

    RooFitResult *fitResult=0;
    RooFitResult *best_fitResult=0;
    Double_t best_chi2 = 100;
    Int_t best_fit = 0;

    TString minimizer = "Minuit2";

    int attempt = 0;
    do {
        std::cout<<">>> Fit with attempt "<<attempt<<std::endl;
        // different fit attempts:
        // attempt 0: fit background in sideband, then full pdf in full range
        // attempt 1: fit signal in center, then full pdf in full range
        // attempt 2: fit background in sideband, then signal in center, then full pdf in full range

        if(attempt > 0){

            // reset parameters of fit models to initial values
            sigModel->Reset();
            bkgModel->Reset();

            Nsig.setVal(NsigInit);
            Nbkg.setVal(NbkgInit);

            if(attempt==3){
                if(effType == "Trk"){
                    // define regions different
                    m.setRange("rangeLow", massLo, 76);
                    m.setRange("rangeHigh", 104, massHi);
                    m.setRange("rangeCenter", 76, 104);
                }
                else{
                    // set sideband range for background fit
                    m.setRange("rangeLow", massLo, 76);
                    m.setRange("rangeHigh", 106, massHi);
                    m.setRange("rangeCenter", 81, 101);
                }
            }
        }

        if(attempt != 1){

            bkgModel->model->fitTo(*data,
                RooFit::Range("rangeLow,rangeHigh"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Strategy(2), // MINOS STRATEGY
                (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer)
                );
        }

        if(attempt > 0){
            // fit signal shape to central region only

            // bkgModel->Print();
            // bkgModel->freeze_all_parameters(1);

            sigModel->model->fitTo(*data,
                RooFit::Range("rangeCenter"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer)
                );

            // bkgModel->Print();
            // bkgModel->freeze_all_parameters(0);
        }

        fitResult = modelPdf.fitTo(*data,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(2), // MINOS STRATEGY
            (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
            // RooFit::IntegrateBins(integrateBinsPrecision),
            RooFit::Minimizer(minimizer),
            // RooFit::Minos(RooArgSet()),
            RooFit::Save()
            );

        RooChi2Var chi2Var("chi2", "chi 2", modelPdf, *data,
            RooFit::DataError(RooAbsData::Expected) //use Expected contribution from model PDF to calculate uncertainty
        );

        // reduced chi2 value from number of degree of freedom in the fit nDoF = nBins - nParams
        Double_t chi2ndf = chi2Var.getVal() / (data->numEntries() - fitResult->floatParsFinal().getSize());

        if(std::abs(1 - (Nsig.getVal()+Nbkg.getVal())/h_yield->Integral()) > 0.1){
            std::cout<<"WARNING: something went wrong in the fit, we give a bad chi2"<<std::endl;
            chi2ndf = 99;
        }  

        std::cout<<"---------------------------------------" <<std::endl;
        std::cout<<"------ attempt = " << attempt << std::endl;
        std::cout<<"------ chi2/ndf = " << chi2ndf <<std::endl;
        std::cout<<"------ status = " << fitResult->status() <<std::endl;
        std::cout<<"---------------------------------------" <<std::endl;

        if(attempt==0 || (fitResult->status() < best_fitResult->status()) || ((fitResult->status() == best_fitResult->status()) && (chi2ndf < best_chi2))){
            // always take the first attempt,
            // take the other attempt if the fit status os lower (better)
            // of if the fit status is the same, take the other attempt if it has lower chi2
            best_chi2 = chi2ndf;
            best_fit = attempt;
            best_fitResult = fitResult;
            w->saveSnapshot(("snapshot_"+std::to_string(attempt)).c_str(),
                fitResult->floatParsFinal(), kTRUE);
        }

        attempt++;

    } while(attempt < 4 && (best_chi2 > 1.5 || best_fitResult->status() != 0));

    // load best fit values into workspace
    w->loadSnapshot(("snapshot_"+std::to_string(best_fit)).c_str());

    const Double_t chi2 = make_plot(nfl, m, data, sigModel, bkgModel,
        h_yield->Integral(), 
        w->pdf("totalPdf"),
        (RooRealVar*)best_fitResult->floatParsFinal().find("Nsig"), 
        (RooRealVar*)best_fitResult->floatParsFinal().find("Nbkg"),
        iBin, 0, 0,
        effType.Data(), etaRegion.Data(), passRegion);

    RooRealVar chi2_plot("chi2plot","chi2 from plot",chi2);
    RooRealVar chi2_best("chi2","chi2",best_chi2);
    w->import(chi2_plot);
    w->import(chi2_best);

    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2 = " << chi2 <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    delete sigModel;
    delete bkgModel;
    delete data;
    delete histfile;

    if(constraint_sig_mean!=nullptr) delete constraint_sig_mean;
    if(constraint_sig_sigma!=nullptr) delete constraint_sig_sigma;
    if(constraint_sig_alpha!=nullptr) delete constraint_sig_alpha;
    if(constraint_alpha!=nullptr) delete constraint_alpha;
    if(constraint_beta!=nullptr) delete constraint_beta;
    if(constraint_gamma!=nullptr) delete constraint_gamma;

}

//--------------------------------------------------------------------------------------------------
// perform fit on two histograms for tag and probe efficiency
void calculateDataEfficiency(
        TH1D          *passHist,            // histogram with passing probes
        TH1D          *failHist,            // histogram with failing probes
        const Int_t   iBin,                 // Label of measurement number in currect run
        const TString effType,              // "HLT" or "Sel" or "Glo" or "Sta" or "Trk"
        const TString etaRegion,            // Barrel "B", Endcap "E" or Inclusive "I"
        const Int_t   sigpass,              // signal model for PASS sample
        const Int_t   bkgpass,              // background model for PASS sample
        const Int_t   sigfail,              // signal model for FAIL sample
        const Int_t   bkgfail,              // background model for FAIL sample
        TH1D          *hPV=0,
        const TString mcfilename="",        // ROOT file containing MC events to generate templates from
        const TString bkgQCDFilename="",    // ROOT file containing bkg template
        const TString bkgTTFilename=""
){

    std::cout<<">>> Do fit in "<< etaRegion <<" for "<<effType<<std::endl;

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(massBin);
    // This is only needed for the convolution of template with an analytic function 
    // set to a high number for good precision
    m.setBins(10000,"cache");  
    // set sideband range for background fit
    m.setRange("rangeLow", massLo, 76);
    m.setRange("rangeHigh", 106, massHi);
    m.setRange("rangeCenter", 81, 101);

    RooCategory sample("sample","");
    sample.defineType("Pass",1);
    sample.defineType("Fail",2);

    RooDataHist* dataPass     = new RooDataHist("dataPass","dataPass",RooArgSet(m),passHist);
    RooDataHist* dataFail     = new RooDataHist("dataFail","dataFail",RooArgSet(m),failHist);
    RooDataHist* dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("Pass",*dataPass),
        RooFit::Import("Fail",*dataFail));

    TFile *histfile = 0;
    if(sigpass%2 == 0 || sigfail%2 == 0) {
        histfile = generateTemplate(mcfilename, effType, hPV);
        assert(histfile);
    }
    std::vector<double> vBkgPars;
    if(bkgfail== 2 or bkgfail==3){
        vBkgPars = preFit(failHist);
    }
    TFile *histbkgQCDfile = 0;
    // TFile *histbkgTTfile = 0;
    if(bkgfail== 7 or bkgpass==7){
        histbkgQCDfile = new TFile(bkgQCDFilename);
        assert(histbkgQCDfile);
        // if(bkgfail== 7 or bkgpass==7){
        //   histbkgTTfile =  new TFile(bkgTTTemplate);
        //   assert(histbkgTTfile);
        // }
    }
    TH1D *hbkgQCDPass = 0;
    TH1D *hbkgQCDFail = 0;
    // TH1D *hbkgTTPass = 0;
    // TH1D *hbkgTTFail = 0;
    if(bkgpass == 7){
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(
            Form("h_mass_%s_pass_%s", effType.Data(), etaRegion.Data()));
        hbkgQCDPass = hbkg->ProjectionY();
        // smoothening - careful: can introduce bias
        // hbkgQCDPass->GetXaxis()->SetRangeUser(massLo,massHi);
        // hbkgQCDPass->Smooth(10000,"R");
        assert(hbkgQCDPass);
        // if(bkgpass == 7){
        //   hbkgTTPass = (TH1D*)histbkgTTfile->Get(Form("bkg_template_%s%s", effType.Data(), "Pass"));
        //   assert(hbkgTTPass);
        // }
    }
    if(bkgfail == 7){
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(
            Form("h_mass_%s_fail_%s", effType.Data(), etaRegion.Data()));
        hbkgQCDFail = hbkg->ProjectionY();
        // smoothening - careful: can introduce bias
        // hbkgQCDFail->GetXaxis()->SetRangeUser(massLo,massHi);
        // hbkgQCDFail->Smooth(10000,"R");
        assert(hbkgQCDFail);
        // if(bkgfail == 7){
        //   hbkgTTFail = (TH1D*)histbkgTTfile->Get(Form("bkg_template_%s%s", effType.Data(), "Fail"));
        //   assert(hbkgTTFail);
        // }
    }

    CSignalModel     *sigPass = 0;
    CBackgroundModel *bkgPass = 0;
    CSignalModel     *sigFail = 0;
    CBackgroundModel *bkgFail = 0;

    Int_t nflpass=0, nflfail=0;

    TH1D *hPass=0;
    if(sigpass%2 == 0) {
        hPass = (TH1D*)histfile->Get(Form("h_mass_pass_%s", etaRegion.Data()));
        hPass->SetDirectory(0);
    }

    nflpass += set_signal_model(sigpass, sigPass, m, kTRUE, 0, hPass);
    nflpass += set_background_model(bkgpass, bkgPass, m, kTRUE, 0, hbkgQCDPass);

    TH1D *hFail=0;
    if(sigfail%2 == 0) {
        hFail = (TH1D*)histfile->Get(Form("h_mass_fail_%s", etaRegion.Data()));
        hFail->SetDirectory(0);
    }

    nflfail += set_signal_model(sigfail, sigFail, m, kFALSE, 0, hFail);
    nflfail += set_background_model(bkgfail, bkgFail, m, kFALSE, 0, hbkgQCDFail, &vBkgPars);

    // ---> Add constraint terms
    // Construct Gaussian constraints p.d.f on signal model
    RooArgSet *constraints = 0;
    RooGaussian *constraint_pass_mean = 0;
    RooGaussian *constraint_pass_sigma = 0;
    RooGaussian *constraint_fail_mean = 0;
    RooGaussian *constraint_fail_sigma = 0;
    if(sigpass == 2 || sigpass == 6 || sigfail == 2 || sigfail == 6){
        // Construct Gaussian constraints p.d.f on signal model
        constraint_pass_mean = new RooGaussian("constraint_pass_mean","constraint mean", *sigPass->mean, RooFit::RooConst(0.), RooFit::RooConst(1.0));
        constraint_pass_sigma = new RooGaussian("constraint_pass_sigma","constraint sigma", *sigPass->sigma, 
            effType == "Trk" ? RooFit::RooConst(2.5) : RooFit::RooConst(0.5), 
            effType == "Trk" ? RooFit::RooConst(1.0) : RooFit::RooConst(0.5));

        constraint_fail_mean = new RooGaussian("constraint_fail_mean","constraint mean", *sigFail->mean, RooFit::RooConst(0.), RooFit::RooConst(1.0));
        constraint_fail_sigma = new RooGaussian("constraint_fail_sigma","constraint sigma", *sigFail->sigma, 
            effType == "Trk" ? RooFit::RooConst(2.5) : RooFit::RooConst(0.5), 
            effType == "Trk" ? RooFit::RooConst(1.0) : RooFit::RooConst(0.5));

        constraints = new RooArgSet(*constraint_pass_mean, *constraint_pass_sigma, *constraint_fail_mean, *constraint_fail_sigma);
    }
    // <---

    Double_t NsigMax     = passHist->Integral()+failHist->Integral();
    Double_t NbkgFailMax = failHist->Integral();
    Double_t NbkgPassMax = passHist->Integral();
    RooRealVar Nsig("Nsig","Signal Yield", NbkgPassMax*0.99, 0, NsigMax);
    RooRealVar eff("eff","Efficiency",0.98,0.,1.0);
    RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.01*NbkgPassMax, 0.0, NbkgPassMax);
    if(bkgpass==0)
        NbkgPass.setVal(0);
    RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.1*NbkgFailMax, 0.0, NbkgFailMax);

    RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
    RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));

    RooAddPdf modelPass("modelPass","Model for PASS sample",
        RooArgList(*(sigPass->model),*(bkgPass->model)),RooArgList(NsigPass,NbkgPass));
    RooAddPdf modelFail("modelFail","Model for FAIL sample",
        RooArgList(*(sigFail->model),*(bkgFail->model)),RooArgList(NsigFail,NbkgFail));

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(modelPass,"Pass");
    totalPdf.addPdf(modelFail,"Fail");

    TFile *fFit = new TFile(
        Form(outputDir+"/workspace_%s_%s_%i.root", effType.Data(), etaRegion.Data(), iBin),
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(*dataCombined);
    w->import(totalPdf);

    RooMsgService::instance().setSilentMode(kTRUE);

    Int_t strategy = 2; // Minuit strategy
    RooFitResult *fitResult=0;
    RooFitResult *best_fitResult=0;
    Double_t best_chi2 = 100;
    Int_t best_fit = 0;

    TString minimizer = "Minuit2";

    int attempt = 0;  // fit attempt
    do {
        std::cout<<">>> Fit with attempt "<<attempt<<std::endl;

        if(attempt>0){
            // reset parameters of fit models to initial values
            bkgPass->Reset();
            bkgFail->Reset();
            sigPass->Reset();
            sigFail->Reset();

            eff.setVal(0.98);
            Nsig.setVal(NbkgPassMax*0.99);

            if(attempt==3){
                // freeze efficiency so each region is fit seprarately first
                eff.setConstant(kTRUE);
                // define regions different
                m.setRange("rangeLow", massLo, 76);
                m.setRange("rangeHigh", 104, massHi);
                m.setRange("rangeCenter", 76, 104);
            }
        }

        if(attempt!=1){
            bkgFail->model->fitTo(*dataFail,
                RooFit::Range("rangeLow,rangeHigh"),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));
        }

        // bkgFail->Print();
        // bkgFail->freeze_all_parameters(1);
        // NbkgFail.setConstant(1);

        if(attempt>0){
            sigFail->model->fitTo(*dataFail,
                RooFit::Range("rangeCenter"),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));
        }

        modelFail.fitTo(*dataFail,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
            RooFit::Strategy(2), // MINOS STRATEGY
            // RooFit::IntegrateBins(integrateBinsPrecision),
            RooFit::Minimizer(minimizer));

        if(attempt!=1){
            bkgPass->model->fitTo(*dataPass,
                RooFit::Range("rangeLow,rangeHigh"),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));
        }

        if(attempt>0){
            sigPass->model->fitTo(*dataPass,
                RooFit::Range("rangeCenter"),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));
        }

        modelPass.fitTo(*dataPass,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
            RooFit::Strategy(2), // MINOS STRATEGY
            // RooFit::IntegrateBins(integrateBinsPrecision),
            RooFit::Minimizer(minimizer));

        if(attempt > 0){    
            const Double_t n0 = NsigFail.getVal();
            const Double_t n1 = NsigPass.getVal();

            const Double_t iEff = n1/(n1+n0);
            const Double_t iNsig = n1 + n0;

            std::cout<<">>> Set parameters to initial values:"<<std::endl;
            std::cout<<"eff = "<<iEff<<std::endl;
            std::cout<<"Nsig = "<<iNsig<<std::endl;
            std::cout<<"<<< "<<std::endl;

            eff.setVal(iEff);
            Nsig.setVal(iNsig);
        }

        if(attempt==3){
            eff.setConstant(kFALSE);
        }

        // bkgFail->Print();
        // bkgFail->freeze_all_parameters(0);
        // NbkgFail.setConstant(0);

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
            RooFit::Strategy(2), // MINOS STRATEGY
            // RooFit::IntegrateBins(integrateBinsPrecision),
            RooFit::Minimizer(minimizer),
            //RooFit::Minos(RooArgSet()),
            RooFit::Save());

        RooChi2Var chi2Var("chi2", "chi 2", totalPdf, *dataCombined,
            RooFit::DataError(RooAbsData::Expected) //use Expected contribution from model PDF to calculate uncertainty
        );

        // reduced chi2 value from number of degree of freedom in the fit nDoF = nBins - nParams
        Double_t chi2ndf = chi2Var.getVal() / (dataCombined->numEntries() - fitResult->floatParsFinal().getSize());

        if(    std::abs(1 - (NsigPass.getVal()+NbkgPass.getVal())/passHist->Integral()) > 0.1
            || std::abs(1 - (NsigFail.getVal()+NbkgFail.getVal())/failHist->Integral()) > 0.1
        ){
            std::cout<<"WARNING: something went wrong in the fit, we give a bad chi2"<<std::endl;
            chi2ndf = 99;
        }

        std::cout<<"---------------------------------------" <<std::endl;
        std::cout<<"------ attempt = " << attempt << std::endl;
        std::cout<<"------ chi2/ndf = " << chi2ndf <<std::endl;
        std::cout<<"------ status = " << fitResult->status() <<std::endl;
        std::cout<<"---------------------------------------" <<std::endl;
        
        if(attempt==0 || (fitResult->status() < best_fitResult->status()) || ((fitResult->status() == best_fitResult->status()) && (chi2ndf < best_chi2))){
            // always take the first attempt,
            // take the other attempt if the fit status os lower (better)
            // of if the fit status is the same, take the other attempt if it has lower chi2
            best_chi2 = chi2ndf;
            best_fit = attempt;
            best_fitResult = fitResult;
            w->saveSnapshot(("snapshot_"+std::to_string(attempt)).c_str(),
                fitResult->floatParsFinal(), kTRUE);
        }

        attempt++;
        // if(attempt==4 && minimizer == "Minuit2"){
        //     // try all again with different minimizer
        //     minimizer = "Minuit";
        //     attempt = 0;
        // }

    } while(attempt < 4 && (best_chi2 > 1.5 || best_fitResult->status() != 0));

    // load best fit values into workspace
    std::cout<<"------ load snapshot " << best_fit <<std::endl;
    w->loadSnapshot(("snapshot_"+std::to_string(best_fit)).c_str());
    
    RooFormulaVar* NsigFormular = (RooFormulaVar*)w->arg("NsigPass");
    RooRealVar* NsigP = new RooRealVar("NsigP","NsigP",NsigFormular->getVal());
    NsigP->setError(NsigFormular->getPropagatedError(*best_fitResult));
    NsigFormular = (RooFormulaVar*)w->arg("NsigFail");
    RooRealVar* NsigF = new RooRealVar("NsigF","NsigF",NsigFormular->getVal());
    NsigF->setError(NsigFormular->getPropagatedError(*best_fitResult));

    const Double_t chi2pass = make_plot(nflpass, m, dataPass,
        sigPass, bkgPass, 
        passHist->Integral(), 
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("Pass"),
        NsigP, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgPass"),
        iBin, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("eff"), 
        0, effType.Data(), etaRegion.Data(), kTRUE);

    const Double_t chi2fail = make_plot(nflfail, m, dataFail,
        sigFail, bkgFail, 
        failHist->Integral(), 
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("Fail"),
        NsigF, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgFail"),
        iBin, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("eff"), 
        0, effType.Data(), etaRegion.Data(), kFALSE);

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2pass = " << chi2pass <<std::endl;
    std::cout<<"------ chi2fail = " << chi2fail <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    RooRealVar chi2p("chi2pass","chi2pass",chi2pass);
    RooRealVar chi2f("chi2fail","chi2fail",chi2fail);
    RooRealVar chi2("chi2","chi2",best_chi2);

    w->import(chi2p);
    w->import(chi2f);
    w->import(chi2);
    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    delete dataCombined;
    delete dataPass;
    delete dataFail;
    delete sigPass;
    delete bkgPass;
    delete sigFail;
    delete bkgFail;
    delete histfile;
    delete histbkgQCDfile;
    // delete histbkgTTfile;

    if(constraint_pass_mean!=nullptr) delete constraint_pass_mean;
    if(constraint_pass_sigma!=nullptr) delete constraint_pass_sigma;
    if(constraint_fail_mean!=nullptr) delete constraint_fail_mean;
    if(constraint_fail_sigma!=nullptr) delete constraint_fail_sigma;
}

//--------------------------------------------------------------------------------------------------
// perform fit in 2 region (1HLT and 2HLT) to extract the Z yield and efficiency together
void calculateHLTEfficiencyAndYield(
    TH1D          *passHist,            // histogram with events where both muon pass HLT
    TH1D          *failHist,            // histogram with events where one muon pass HLT
	const Int_t   iBin,                 // Label of measurement number in currect run
    const TString etaRegion,            // {BB, BE, EE}
	const Int_t   sigpass,              // signal model for PASS sample
	const Int_t   bkgpass,              // background model for PASS sample
	const Int_t   sigfail,              // signal model for FAIL sample
	const Int_t   bkgfail,              // background model for FAIL sample
    // const Double_t corr = 1.,           // correlation factors for second muon
    TH1D          *hPV=0,
    const TString mcfilename="",        // ROOT file containing MC events to generate templates from
    const TString bkgQCDFilename="",    // ROOT file containing bkg template
    const TString bkgTTFilename=""
){

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(massBin);
    // This is only needed for the convolution of template with an analytic function 
    // set to a high number for good precision
    m.setBins(10000,"cache");
    // set sideband range for background fit
    m.setRange("rangeLow", massLo, 76);
    m.setRange("rangeHigh", 104, massHi);
    m.setRange("rangeCenter", 76, 104);

    RooCategory sample("sample","");
    sample.defineType("Pass",1);
    sample.defineType("Fail",2);

    RooDataHist *dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),passHist);
    RooDataHist *dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),failHist);

    RooDataHist *dataCombined = new RooDataHist(
        "dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("Pass",*dataPass),
        RooFit::Import("Fail",*dataFail));

    TFile *histfile = 0;
    if(sigpass%2 == 0 || sigfail%2 == 0) {
        histfile = generateTemplate_ZYield(mcfilename, hPV, iBin);
        assert(histfile);
    }
    const double corr = extractCorrelation_HLT(mcfilename, hPV, etaRegion);

    std::vector<double> vBkgPars;
    if(bkgfail== 2 or bkgfail==3){
        vBkgPars = preFit(failHist);
    }
    TFile *histbkgQCDfile = 0;
    if(bkgfail== 7 or bkgpass==7){
        histbkgQCDfile = new TFile(bkgQCDFilename);
        assert(histbkgQCDfile);
    }
    TH1D *hbkgQCDPass = 0;
    TH1D *hbkgQCDFail = 0;
    if(bkgpass == 7){
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get("h_mass_HLT_pass_"+ etaRegion);
        hbkgQCDPass = hbkg->ProjectionY();
        assert(hbkgQCDPass);
    }
    if(bkgfail == 7){
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get("h_mass_HLT_fail_"+ etaRegion);
        hbkgQCDFail = hbkg->ProjectionY();
        assert(hbkgQCDFail);
    }

    CSignalModel     *sigPass = 0;
    CBackgroundModel *bkgPass = 0;
    CSignalModel     *sigFail = 0;
    CBackgroundModel *bkgFail = 0;

    Int_t nflpass=0, nflfail=0;

    TH1D *hPass=0;
    if(sigpass%2 == 0) {
        // set signal templates of fail histograms
        hPass = (TH1D*)histfile->Get("h_mass_2hlt_"+ etaRegion);
        hPass->SetDirectory(0);
    }

    nflpass += set_signal_model(sigpass, sigPass, m, kTRUE, 2, hPass);
    nflpass += set_background_model(bkgpass, bkgPass, m, kTRUE, 2, hbkgQCDPass);

    TH1D *hFail=0;
    if(sigfail%2 == 0) {
        // set signal templates of pass histograms
        hFail = (TH1D*)histfile->Get("h_mass_1hlt_"+ etaRegion);
        hFail->SetDirectory(0);
    }

    nflfail += set_signal_model(sigfail, sigFail, m, kTRUE, 1, hFail);
    nflfail += set_background_model(bkgfail, bkgFail, m, kTRUE, 1, hbkgQCDFail);

    // ---> Add constraint terms
    RooArgSet *constraints = 0;
    RooGaussian *constraint_pass_mean = 0;
    RooGaussian *constraint_pass_sigma = 0;
    RooGaussian *constraint_fail_mean = 0;
    RooGaussian *constraint_fail_sigma = 0;
    if(sigpass == 2 || sigpass == 6 || sigfail == 2 || sigfail == 6){
        // Construct Gaussian constraints p.d.f on signal model
        constraint_pass_mean = new RooGaussian("constraint_pass_mean","constraint mean", *sigPass->mean, RooFit::RooConst(0.), RooFit::RooConst(1.0));
        constraint_pass_sigma = new RooGaussian("constraint_pass_sigma","constraint sigma", *sigPass->sigma, RooFit::RooConst(0.5), RooFit::RooConst(0.5));

        constraint_fail_mean = new RooGaussian("constraint_fail_mean","constraint mean", *sigFail->mean, RooFit::RooConst(0.), RooFit::RooConst(1.0));
        constraint_fail_sigma = new RooGaussian("constraint_fail_sigma","constraint sigma", *sigFail->sigma, RooFit::RooConst(0.5), RooFit::RooConst(0.5));

        constraints = new RooArgSet(*constraint_pass_mean, *constraint_pass_sigma, *constraint_fail_mean, *constraint_fail_sigma);
    }
    // <---

    Double_t NsigMax     = passHist->Integral() + failHist->Integral();
    Double_t NbkgFailMax = failHist->Integral();
    Double_t NbkgPassMax = passHist->Integral();

    RooRealVar eff("eff","Efficiency barrel",0.95,0.,1);
    RooRealVar Nsig("Nsig","Signal Yield",0.8*NsigMax,0,1.5*NsigMax);
    RooConstVar c("c", "Correlation factor", corr);

    RooFormulaVar NsigPass("NsigPass","eff*eff*Nsig*c",RooArgList(eff,Nsig,c));
    RooFormulaVar NsigFail("NsigFail","2*eff*(1.0-c*eff)*Nsig",RooArgList(eff,Nsig,c));

    RooRealVar NbkgPass("NbkgPass","Background count in  PASS sample",0.01*NbkgPassMax, 0.0, NbkgPassMax);
    RooRealVar NbkgFail("NbkgFail","Background count in  FAIL sample",0.05*NbkgFailMax, 0.0, NbkgFailMax);

    RooAddPdf modelPass("model2","Model for  PASS sample",
        RooArgList(*(sigPass->model),*(bkgPass->model)),
        RooArgList(NsigPass, NbkgPass));
    RooAddPdf modelFail("model1","Model for  FAIL sample",
        RooArgList(*(sigFail->model),*(bkgFail->model)),
        RooArgList(NsigFail, NbkgFail));

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(modelPass,"Pass");
    totalPdf.addPdf(modelFail,"Fail");

    TFile *fFit = new TFile(
        outputDir+"/workspace_"+etaRegion+"_"+iBin+".root",
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(*dataCombined);
    w->import(totalPdf);

    RooMsgService::instance().setSilentMode(kTRUE);

    Int_t strategy = 2; // Minuit strategy
    RooFitResult *fitResult=0;
    RooFitResult *best_fitResult=0;
    Double_t best_chi2 = 100;
    Int_t best_fit = 0;

    TString minimizer = "Minuit2";

    int attempt = 0;  // fit attempt
    do {
        std::cout<<">>> Fit with attempt "<<attempt<<std::endl;

        if(attempt>0){
            // reset parameters of fit models to initial values
            bkgPass->Reset();
            bkgFail->Reset();
            sigPass->Reset();
            sigFail->Reset();

            eff.setVal(0.95);
            Nsig.setVal(NsigMax);

            if(attempt==3){
                eff.setVal(0.8);
                Nsig.setVal(NsigMax*0.5);
                // freeze efficiency so each region is fit seprarately first
                eff.setConstant(kTRUE);

                m.setRange("rangeLow", massLo, 79);
                m.setRange("rangeHigh", 103, massHi);
                m.setRange("rangeCenter", 79, 103);
            }
            else{
                m.setRange("rangeLow", massLo, 81);
                m.setRange("rangeHigh", 101, massHi);
                m.setRange("rangeCenter", 81, 101);
            }
        }
        
        if(attempt!=1){
                std::cout<<">>> Fit sideband regions in fail"<<std::endl;

            bkgFail->model->fitTo(*dataFail,
                RooFit::Range("rangeLow,rangeHigh"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));
        }

        // bkgFail->Print();
        // bkgFail->freeze_all_parameters(1);

        if(attempt > 0){
            sigFail->model->fitTo(*dataFail,
                RooFit::Range("rangeCenter"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));

        }

        // fit total pdf in Fail
        modelFail.fitTo(*dataFail,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
            RooFit::Strategy(2), // MINOS STRATEGY
            // RooFit::IntegrateBins(integrateBinsPrecision),
            RooFit::Minimizer(minimizer));


        if(attempt!=1){
                std::cout<<">>> Fit sideband regions in fail"<<std::endl;

            bkgPass->model->fitTo(*dataPass,
                RooFit::Range("rangeLow,rangeHigh"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));
        }

        // bkgPass->Print();
        // bkgPass->freeze_all_parameters(1);

        if(attempt > 0){
            sigPass->model->fitTo(*dataPass,
                RooFit::Range("rangeCenter"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
                RooFit::Strategy(2), // MINOS STRATEGY
                // RooFit::IntegrateBins(integrateBinsPrecision),
                RooFit::Minimizer(minimizer));

        }

        // fit total pdf in Pass
        modelPass.fitTo(*dataPass,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
            RooFit::Strategy(2), // MINOS STRATEGY
            // RooFit::IntegrateBins(integrateBinsPrecision),
            RooFit::Minimizer(minimizer));

        if(attempt>0){
            std::cout<<">>> Fit sideband regions in pass"<<std::endl;

            const Double_t n1 = NsigFail.getVal();
            const Double_t n2 = NsigPass.getVal();
            const Double_t iEff = 2*n2/(corr*(n1+2*n2));
            const Double_t iNsig = n2/(corr * iEff * iEff);

            std::cout<<">>> Set parameters to initial values:"<<std::endl;
            std::cout<<"eff = "<<iEff<<std::endl;
            std::cout<<"Nsig = "<<iNsig<<std::endl;
            std::cout<<"<<< "<<std::endl;

            if(attempt==3){
                eff.setConstant(kFALSE);
            }

            eff.setVal(iEff);
            Nsig.setVal(iNsig);
        }

        // bkgFail->Print();
        // bkgPass->Print();
        // bkgFail->freeze_all_parameters(0);
        // bkgPass->freeze_all_parameters(0);

        // fit all regions together
        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            (constraints != nullptr) ? RooFit::ExternalConstraints(*constraints) : RooCmdArg::none(),
            RooFit::Strategy(2), // MINOS STRATEGY
            // RooFit::IntegrateBins(integrateBinsPrecision),
            RooFit::Minimizer(minimizer),
            //RooFit::Minos(RooArgSet()),
            RooFit::Save());

        RooChi2Var chi2Var("chi2", "chi 2", totalPdf, *dataCombined,
            RooFit::DataError(RooAbsData::Expected) //use Expected contribution from model PDF to calculate uncertainty
        );

        // reduced chi2 value from number of degree of freedom in the fit nDoF = nBins - nParams
        const Int_t ndf = dataCombined->numEntries() - fitResult->floatParsFinal().getSize();
        Double_t chi2ndf = chi2Var.getVal() / ndf;

        if(    std::abs(1 - (NsigPass.getVal()+NbkgPass.getVal())/passHist->Integral()) > 0.1
            || std::abs(1 - (NsigFail.getVal()+NbkgFail.getVal())/failHist->Integral()) > 0.1
        ){
            std::cout<<"WARNING: something went wrong in the fit, we give a bad chi2"<<std::endl;
            chi2ndf = 99;
        }        

        std::cout<<"---------------------------------------" <<std::endl;
        std::cout<<"------ attempt = " << attempt << std::endl;
        std::cout<<"------ chi2/"<<ndf<<" = " << chi2ndf <<std::endl;
        std::cout<<"------ status = " << fitResult->status() <<std::endl;
        std::cout<<"------ eff = "<<eff.getVal()<<std::endl;
        std::cout<<"------ Nsig = "<<Nsig.getVal()<<std::endl;          
        std::cout<<"---------------------------------------" <<std::endl;

        if(attempt==0 || (fitResult->status() < best_fitResult->status()) || ((fitResult->status() == best_fitResult->status()) && (chi2ndf < best_chi2))){
            // always take the first attempt,
            // take the other attempt if the fit status os lower (better)
            // of if the fit status is the same, take the other attempt if it has lower chi2
            best_chi2 = chi2ndf;
            best_fit = attempt;
            best_fitResult = fitResult;
            w->saveSnapshot(("snapshot_"+std::to_string(attempt)).c_str(),
                fitResult->floatParsFinal(), kTRUE);
        }

        attempt++;

    } while(attempt < 4 && (best_chi2 > 1.5 || best_fitResult->status() != 0));

    // load best fit values into workspace
    w->loadSnapshot(("snapshot_"+std::to_string(best_fit)).c_str());

    RooFormulaVar* NsigFormular = (RooFormulaVar*)w->arg("NsigPass");
    RooRealVar* NsigP = new RooRealVar("NsigP","NsigP",NsigFormular->getVal());
    NsigP->setError(NsigFormular->getPropagatedError(*best_fitResult));
    NsigFormular = (RooFormulaVar*)w->arg("NsigFail");
    RooRealVar* NsigF = new RooRealVar("NsigF","NsigF",NsigFormular->getVal());
    NsigF->setError(NsigFormular->getPropagatedError(*best_fitResult));

    const Double_t chi2pass = make_plot(nflpass, m, dataPass,
        sigPass, bkgPass,
        passHist->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("Pass"),
        NsigP, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgPass"),
        iBin, (RooRealVar*)best_fitResult->floatParsFinal().find("eff"), 
        (RooRealVar*)&c, "yield", etaRegion.Data(), 2);

    const Double_t chi2fail = make_plot(nflfail, m, dataFail,
        sigFail, bkgFail, 
        failHist->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("Fail"),
        NsigF, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgFail"),
        iBin, (RooRealVar*)best_fitResult->floatParsFinal().find("eff"), 
        (RooRealVar*)&c, "yield", etaRegion.Data(), 1);

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2pass = " << chi2pass <<std::endl;
    std::cout<<"------ chi2fail = " << chi2fail <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    RooRealVar chi2p("chi2pass","chi2pass",chi2pass);
    RooRealVar chi2f("chi2fail","chi2fail",chi2fail);
    RooRealVar chi2("chi2","chi2",best_chi2);

    w->import(chi2p);
    w->import(chi2f);
    w->import(chi2);
    // w->import(c);
    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    delete dataCombined;
    delete dataPass;
    delete dataFail;
    delete sigPass;
    delete bkgPass;
    delete sigFail;
    delete bkgFail;

    if(constraint_pass_mean!=nullptr) delete constraint_pass_mean;
    if(constraint_pass_sigma!=nullptr) delete constraint_pass_sigma;
    if(constraint_fail_mean!=nullptr) delete constraint_fail_mean;
    if(constraint_fail_sigma!=nullptr) delete constraint_fail_sigma;
}

//--------------------------------------------------------------------------------------------------
// perform fit in 3 regions: 2HLT, 1HLT, Sel fail
// extract the Z yield and efficiencies together
void calculateEfficienciesAndYield(
    TH1D          *h_HLT2,            // histogram with events where both muon pass HLT
    TH1D          *h_HLT1,            // histogram with events where one muon pass HLT
    TH1D          *h_SelFail,         // histogram with events where probe fails selection
    const Int_t   iBin,               // Label of measurement number in currect run
    const TString etaRegion,          // {BB, BE, EE}
    const Int_t   sig,                // signal model
    const Int_t   bkg,                // background model
    TH1D          *hPV=0,
    const TString mcfilename="",      // ROOT file containing MC events to generate templates from
    const TString bkgQCDFilename="",  // ROOT file containing bkg template
    const TString bkgTTFilename=""
){
    // --- prepare data

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(massBin);
    // This is only needed for the convolution of template with an analytic function 
    // set to a high number for good precision
    m.setBins(10000,"cache");  
    // set sideband range for background fit
    m.setRange("rangeLow", massLo, 76);
    m.setRange("rangeHigh", 106, massHi);

    RooCategory sample("sample","");
    sample.defineType("HLT1",1);
    sample.defineType("HLT2",2);
    sample.defineType("SelFail",3);

    RooDataHist *dataHLT1 = new RooDataHist("dataHLT1","dataHLT1",RooArgSet(m),h_HLT1);
    RooDataHist *dataHLT2 = new RooDataHist("dataHLT2","dataHLT2",RooArgSet(m),h_HLT2);
    RooDataHist *dataSelFail = new RooDataHist("dataSelFail","dataSelFail",RooArgSet(m),h_SelFail);

    RooDataHist *dataCombined = new RooDataHist(
        "dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("HLT1",*dataHLT1),
        RooFit::Import("HLT2",*dataHLT2),
        RooFit::Import("SelFail",*dataSelFail)
    );
    
    // --- prepare fit models
    const double icHLT = extractCorrelation_HLT(mcfilename, hPV, etaRegion);

    CSignalModel     *sigHLT1 = 0;
    CSignalModel     *sigHLT2 = 0;
    CSignalModel     *sigSelFail = 0;
    CBackgroundModel     *bkgHLT1 = 0;
    CBackgroundModel     *bkgHLT2 = 0;
    CBackgroundModel     *bkgSelFail = 0;

    Int_t nflHLT1, nflHLT2, nflSelFail; 

    TH1D *hHLT1=0;
    TH1D *hHLT2=0;
    TH1D *hSelFail=0;
    if(sig%2==0) {
        TFile *histfileHLT = generateTemplate_ZYield(mcfilename, 0, iBin);
        TFile *histfileSel = generateTemplate(mcfilename, "Sel", 0);
        assert(histfileHLT);
        // set signal templates histograms
        hHLT1 = (TH1D*)histfileHLT->Get("h_mass_1hlt_"+ etaRegion);
        hHLT2 = (TH1D*)histfileHLT->Get("h_mass_2hlt_"+ etaRegion);
        hSelFail = (TH1D*)histfileSel->Get("h_mass_fail_"+ etaRegion);

        hHLT1->SetDirectory(0);
        hHLT2->SetDirectory(0);
        hSelFail->SetDirectory(0);
        
        histfileHLT->Close();
        histfileSel->Close();
    }

    nflHLT1 = set_signal_model(sig, sigHLT1, m, kTRUE, 1, hHLT1);
    nflHLT2 = set_signal_model(sig, sigHLT2, m, kTRUE, 2, hHLT2);
    nflSelFail = set_signal_model(sig, sigSelFail, m, kFALSE, 0, hSelFail);

    nflHLT1 += set_background_model(bkg, bkgHLT1, m, kTRUE, 1);
    nflHLT2 += set_background_model(bkg, bkgHLT2, m, kTRUE, 2);
    nflSelFail += set_background_model(bkg, bkgSelFail, m, kFALSE, 0);

    // --- set fit model
    const Double_t nHLT1Max = h_HLT1->Integral();
    const Double_t nHLT2Max = h_HLT2->Integral();
    const Double_t nSelFailMax = h_SelFail->Integral();

    const Double_t NsigMax = nHLT1Max + nHLT2Max + nSelFailMax;
    
    // initial guesses
    const Double_t nHLT1ExpSig = nHLT1Max*0.99;
    const Double_t nHLT2ExpSig = nHLT2Max*0.995;
    const Double_t nSelFailExpSig = nSelFailMax*0.8;
    
    const Double_t ieffHLT  =  2*nHLT2ExpSig                 / (2*nHLT2ExpSig + nHLT1ExpSig);
    const Double_t ieffSel  = (2*nHLT2ExpSig + nHLT1ExpSig ) / (2*nHLT2ExpSig + nHLT1ExpSig + nSelFailExpSig);

    const Double_t iNsig = (std::pow(nHLT2ExpSig+nHLT1ExpSig/2.,2) * icHLT / nHLT2ExpSig ) / std::pow(ieffSel,2);
    
    // free fit parameters
    RooRealVar Nsig("Nsig","Signal yield",              iNsig,  0., 5*NsigMax);
    RooRealVar effHLT("effHLT","HLT efficiency",        ieffHLT, 0.,1./icHLT);
    RooRealVar effSel("effSel","Selection Efficiency",  ieffSel, 0.,1);

    RooRealVar NbkgHLT1("NbkgHLT1",      "Background count in  HLT1 sample",    nHLT1Max-nHLT1ExpSig,       0., nHLT1Max);
    RooRealVar NbkgHLT2("NbkgHLT2",      "Background count in  HLT2 sample",    nHLT2Max-nHLT2ExpSig,       0., nHLT2Max);    
    RooRealVar NbkgSelFail("NbkgSelFail","Background count in  SelFail sample", nSelFailMax-nSelFailExpSig, 0., nSelFailMax);

    RooConstVar cHLT("cHLT", "Correlation factor between HLT muons", icHLT);

    RooFormulaVar NsigHLT2("NsigHLT2","effHLT*effHLT*cHLT*Nsig*effSel*effSel",        RooArgList(Nsig, effHLT ,cHLT, effSel));
    RooFormulaVar NsigHLT1("NsigHLT1","2*effHLT*(1.0-cHLT*effHLT)*Nsig*effSel*effSel",RooArgList(Nsig, effHLT ,cHLT, effSel));
    RooFormulaVar NsigSelFail("NsigSelFail","2*effHLT*Nsig*effSel*(1-effSel)",        RooArgList(Nsig, effHLT ,cHLT, effSel));

    RooAddPdf mHLT1("mHLT1","Model for HLT1 sample",
        RooArgList(*(sigHLT1->model), *(bkgHLT1->model)),
        RooArgList(NsigHLT1, NbkgHLT1));
    RooAddPdf mHLT2("mHLT2","Model for HLT2 sample",
        RooArgList(*(sigHLT2->model), *(bkgHLT2->model)),
        RooArgList(NsigHLT2, NbkgHLT2));
    RooAddPdf mSelFail("mSelFail","Model for SelFail sample",
        RooArgList(*(sigSelFail->model), *(bkgSelFail->model)),
        RooArgList(NsigSelFail, NbkgSelFail));

    RooSimultaneous totalPdf("totalPdf","totalPdf", sample);
    totalPdf.addPdf(mHLT1,"HLT1");
    totalPdf.addPdf(mHLT2,"HLT2");
    totalPdf.addPdf(mSelFail,"SelFail");

    // --- set workspace for saving results
    
    TFile *fFit = new TFile(
        outputDir+"/workspace_"+etaRegion+"_"+iBin+".root",
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(*dataCombined);
    w->import(totalPdf);
    
    Int_t strategy = 2; // Minuit strategy
    RooFitResult *fitResult=0;
    RooFitResult *best_fitResult=0;
    Double_t best_chi2 = 100;
    Int_t best_fit = 0;

    int attempt = 0;  // fit attempt

    RooMsgService::instance().setSilentMode(kTRUE);

    do {
        std::cout<<">>> Fit with attempt "<<attempt<<std::endl;
        // attempt 0: just fit fail and pass region simultaneously in full range
        // attempt 1: first fit fail and pass pdf separately,
        //    then set proper initial start values
        //    then fit both together,
        // attempt 2: first fit fail pdf in sideband range,
        //    then in complete fail range,
        //    same for pass region
        //    then set proper initial start values
        //    then fit fail and pass region simultaneously in full range
        
        if(attempt>0){
            // reset parameters of fit models to initial values
            sigHLT1->Reset();
            sigHLT2->Reset();
            sigSelFail->Reset();
            bkgHLT1->Reset();
            bkgHLT2->Reset();
            bkgSelFail->Reset();
            
            effHLT.setVal(0.90);
            effSel.setVal(0.90);
            Nsig.setVal(NsigMax);

            if(attempt==3){
                effHLT.setVal(0.80);
                effSel.setVal(0.80);
                Nsig.setVal(NsigMax*0.5);
                // freeze efficiency so each region is fit seprarately first
                effHLT.setConstant(kTRUE);
                effSel.setConstant(kTRUE);

                m.setRange("rangeLow", massLo, 79);
                m.setRange("rangeHigh", 103, massHi);
                m.setRange("rangeCenter", 79, 103);
            }
            else{
                m.setRange("rangeLow", massLo, 81);
                m.setRange("rangeHigh", 101, massHi);
                m.setRange("rangeCenter", 81, 101);
            }

            if(attempt==2){
                std::cout<<">>> Fit sideband regions in fail"<<std::endl;

                fitResult = bkgHLT1->model->fitTo(*dataHLT1,
                    RooFit::Range("rangeLow,rangeHigh"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());

                fitResult = sigHLT1->model->fitTo(*dataHLT1,
                    RooFit::Range("rangeCenter"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());

            }
            std::cout<<">>> Fit full range pdf in fail "<<std::endl;
            // fit total pdf in Fail
            fitResult = mHLT1.fitTo(*dataHLT1,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                // RooFit::Minos(RooArgSet(Nsig)),
                RooFit::Save());

            const Double_t n1 = NsigHLT1.getVal();

            if(attempt==2){
                std::cout<<">>> Fit sideband regions in pass"<<std::endl;

                fitResult = bkgHLT2->model->fitTo(*dataHLT2,
                    RooFit::Range("rangeLow,rangeHigh"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());

                fitResult = sigHLT2->model->fitTo(*dataHLT2,
                    RooFit::Range("rangeCenter"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());
            }

            std::cout<<">>> Fit full range pdf in pass "<<std::endl;
            // fit total pdf in Fail
            fitResult = mHLT2.fitTo(*dataHLT2,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                // RooFit::Minos(RooArgSet(Nsig)),
                RooFit::Save());

            const Double_t n2 = NsigHLT2.getVal();

            if(attempt==2){
                std::cout<<">>> Fit sideband regions in pass"<<std::endl;

                fitResult = bkgSelFail->model->fitTo(*dataSelFail,
                    RooFit::Range("rangeLow,rangeHigh"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());

                fitResult = sigSelFail->model->fitTo(*dataSelFail,
                    RooFit::Range("rangeCenter"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());
            }

            std::cout<<">>> Fit full range pdf in pass "<<std::endl;
            // fit total pdf in Fail
            fitResult = mSelFail.fitTo(*dataSelFail,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                // RooFit::Minos(RooArgSet(Nsig)),
                RooFit::Save());

            const Double_t nSelFail = NsigSelFail.getVal();
            
            const Double_t iEffHLT  =  2*n2        / (2*n2 + n1);
            const Double_t iEffSel  = (2*n2 + n1 ) / (2*n2 + n1 + nSelFail);
            const Double_t iNsig = (std::pow(n2+n1/2.,2) * icHLT / n2 ) / std::pow(ieffSel,2);


            std::cout<<">>> Set parameters to initial values:"<<std::endl;
            std::cout<<"eff(HLT) = "<<iEffHLT<<std::endl;
            std::cout<<"eff(Sel) = "<<iEffSel<<std::endl;
            std::cout<<"Nsig = "<<iNsig<<std::endl;
            std::cout<<"<<< "<<std::endl;

            if(attempt==3){
                effHLT.setConstant(kFALSE);
                effSel.setConstant(kFALSE);
            }

            effHLT.setVal(iEffHLT);
            effSel.setVal(iEffSel);
            Nsig.setVal(iNsig);
        }

        // first fit each region separately
        mHLT1.fitTo(*dataHLT1,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(strategy) // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(Nsig)),
        );
        mHLT2.fitTo(*dataHLT2,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(strategy) // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(Nsig)),
        );
        mSelFail.fitTo(*dataSelFail,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(strategy) // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(Nsig)),
        );

        RooFitResult *fitResult=0;

        // now fit the full model
        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(0), // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(1), // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(Nsig, eff)),
            RooFit::Save());

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(strategy), // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(Nsig, eff)),
            RooFit::Save());


        RooChi2Var chi2Var("chi2", "chi 2", totalPdf, *dataCombined,
            RooFit::DataError(RooAbsData::Expected) //use Expected contribution from model PDF to calculate uncertainty
        );

        // reduced chi2 value from number of degree of freedom in the fit nDoF = nBins - nParams
        const Int_t ndf = dataCombined->numEntries() - fitResult->floatParsFinal().getSize();
        Double_t chi2ndf = chi2Var.getVal() / ndf;

        std::cout<<"---------------------------------------" <<std::endl;
        std::cout<<"------ attempt = " << attempt << std::endl;
        std::cout<<"------ chi2/"<<ndf<<" = " << chi2ndf <<std::endl;
        std::cout<<"------ eff(HLT) = "<<effHLT.getVal()<<std::endl;
        std::cout<<"------ eff(Sel) = "<<effSel.getVal()<<std::endl;
        std::cout<<"------ Nsig = "<<Nsig.getVal()<<std::endl;          
        std::cout<<"---------------------------------------" <<std::endl;

        if(NsigHLT2.getVal() > h_HLT2->Integral() || NsigHLT1.getVal() > h_HLT1->Integral()){
            std::cout<<"WARNING: something went wrong in the fit, we give a bad chi2"<<std::endl;
            chi2ndf = 99;
        }        

        if(attempt==0 || (chi2ndf < best_chi2 && chi2ndf > 0) ){
            best_chi2 = chi2ndf;
            best_fit = attempt;
            best_fitResult = fitResult;
            w->saveSnapshot(("snapshot_"+std::to_string(attempt)).c_str(),
                fitResult->floatParsFinal(), kTRUE);
        }

        attempt++;

    } while(attempt < 4 && best_chi2 > 1.5);

    // load best fit values into workspace
    w->loadSnapshot(("snapshot_"+std::to_string(best_fit)).c_str());

    RooFormulaVar* NsigFormular = (RooFormulaVar*)w->arg("NsigHLT2");
    RooRealVar* NsHLT2 = new RooRealVar("NsigP","NsigP",NsigFormular->getVal());
    NsHLT2->setError(NsigFormular->getPropagatedError(*best_fitResult));
    
    NsigFormular = (RooFormulaVar*)w->arg("NsigHLT1");
    RooRealVar* NsHLT1 = new RooRealVar("NsigF","NsigF",NsigFormular->getVal());
    NsHLT1->setError(NsigFormular->getPropagatedError(*best_fitResult));

    NsigFormular = (RooFormulaVar*)w->arg("NsigSelFail");
    RooRealVar* NsSelFail = new RooRealVar("NsigSelF","NsigSelF",NsigFormular->getVal());
    NsSelFail->setError(NsigFormular->getPropagatedError(*best_fitResult));

    const Double_t chi2hlt2 = make_plot(nflHLT2, m, dataHLT2,
        sigHLT2, bkgHLT2,
        h_HLT2->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("HLT2"),
        NsHLT2, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgHLT2"),
        iBin, (RooRealVar*)best_fitResult->floatParsFinal().find("effHLT"), 
        (RooRealVar*)&cHLT, "yield", etaRegion.Data(), 2);

    const Double_t chi2hlt1 = make_plot(nflHLT1, m, dataHLT1,
        sigHLT1, bkgHLT1, 
        h_HLT1->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("HLT1"),
        NsHLT1, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgHLT1"),
        iBin, (RooRealVar*)best_fitResult->floatParsFinal().find("effHLT"), 
        (RooRealVar*)&cHLT, "yield", etaRegion.Data(), 1);

    const Double_t chi2fail = make_plot(nflSelFail, m, dataSelFail,
        sigSelFail, bkgSelFail, 
        h_SelFail->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("SelFail"),
        NsSelFail, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgSelFail"),
        iBin, (RooRealVar*)best_fitResult->floatParsFinal().find("effSel"), 
        0, "Sel", etaRegion.Data(), 0);

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2hlt2 = " << chi2hlt2 <<std::endl;
    std::cout<<"------ chi2hlt1 = " << chi2hlt1 <<std::endl;
    std::cout<<"------ chi2fail = " << chi2fail <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    RooRealVar vChi2hlt2("chi2hlt2","chi2hlt2",chi2hlt2);
    RooRealVar vChi2hlt1("chi2hlt1","chi2hlt1",chi2hlt1);
    RooRealVar vChi2fail("chi2fail","chi2fail",chi2fail);
    RooRealVar vChi2("chi2","chi2",best_chi2);

    w->import(vChi2hlt2);
    w->import(vChi2hlt1);
    w->import(vChi2fail);
    w->import(vChi2);
    // w->import(cHLT);
    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    delete dataCombined;
    delete dataHLT1;
    delete dataHLT2;
    delete dataSelFail;
    delete sigHLT1;
    delete sigHLT2;
    delete sigSelFail;
    delete bkgHLT1;
    delete bkgHLT2;
    delete bkgSelFail;
}


//--------------------------------------------------------------------------------------------------
// perform fit in 3 region (0HLT, 1HLT and 2HLT) to extract the Z yield, efficiency and correlation together
//  the 0 HLT histogram has to be taken from an unbiased trigger selection e.g. MET dataset
//  Therefore this method is used to estimate the correlation coefficient
void calculateHLTCorrelation(
    TH1          *h_HLT0,            // histogram with events where no muon pass HLT
    TH1          *h_HLT1,            // histogram with events where one muon pass HLT
    TH1          *h_HLT2,            // histogram with events where one muon pass HLT
	const Int_t   iBin,                 // Label of measurement number in currect run
    const TString etaRegion,            // {BB, BE, EE}
	const Int_t   sig,                 // signal model
	const Int_t   bkg,                 // background model
    TH1D          *hPV=0,
    const TString mcfilename=""        // ROOT file containing MC events to generate templates from
){

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(massBin);
    // This is only needed for the convolution of template with an analytic function 
    // set to a high number for good precision
    m.setBins(10000,"cache");  
    // set sideband range for background fit
    m.setRange("rangeLow", massLo, 76);
    m.setRange("rangeHigh", 106, massHi);

    RooCategory sample("sample","");
    sample.defineType("HLT0",1);
    sample.defineType("HLT1",2);
    sample.defineType("HLT2",3);

    RooDataHist *dataHLT0 = new RooDataHist("dataHLT0","dataHLT0",RooArgSet(m),h_HLT0);
    RooDataHist *dataHLT1 = new RooDataHist("dataHLT1","dataHLT1",RooArgSet(m),h_HLT1);
    RooDataHist *dataHLT2 = new RooDataHist("dataHLT2","dataHLT2",RooArgSet(m),h_HLT2);

    RooDataHist *dataCombined = new RooDataHist(
        "dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("HLT0",*dataHLT0),
        RooFit::Import("HLT1",*dataHLT1),
        RooFit::Import("HLT2",*dataHLT2));


    const double corr = extractCorrelation_HLT(mcfilename, hPV, etaRegion);
    std::vector<double> vBkgPars;
    if(bkg == 2 or bkg == 3){
        vBkgPars = preFit(h_HLT0);
    }

    CSignalModel     *sig0 = 0;
    CBackgroundModel *bkg0 = 0;
    CSignalModel     *sig1 = 0;
    CBackgroundModel *bkg1 = 0;
    CSignalModel     *sig2 = 0;
    CBackgroundModel *bkg2 = 0;

    Int_t nfl0, nfl1, nfl2;

    TH1D *h0=0;
    TH1D *h1=0;
    TH1D *h2=0;
    if(sig%2 == 0) {
        TFile *histfile = generateTemplate_ZYield(mcfilename, 0, iBin);
        assert(histfile);
        // set signal templates of fail histograms
        h0 = (TH1D*)histfile->Get("h_mass_0hlt_"+ etaRegion);
        h0->SetDirectory(0);
        // set signal templates of pass histograms
        h2 = (TH1D*)histfile->Get("h_mass_2hlt_"+ etaRegion);
        h2->SetDirectory(0);
        // set signal templates of pass histograms
        h1 = (TH1D*)histfile->Get("h_mass_1hlt_"+ etaRegion);
        h1->SetDirectory(0);
        
        histfile->Close();
    }

    nfl0 = set_signal_model(sig, sig0, m, kTRUE, 0, h0);
    nfl0 += set_background_model(bkg, bkg0, m, kTRUE, 0, 0);
    nfl1 = set_signal_model(sig, sig1, m, kTRUE, 1, h1);
    nfl1 += set_background_model(bkg, bkg1, m, kTRUE, 1, 0);
    nfl2 = set_signal_model(sig, sig2, m, kTRUE, 2, h2);
    nfl2 += set_background_model(bkg, bkg2, m, kTRUE, 2, 0);

    Double_t NsigMax  = h_HLT0->Integral() + h_HLT1->Integral() + h_HLT2->Integral();
    Double_t Nbkg0Max = h_HLT0->Integral();
    Double_t Nbkg1Max = h_HLT1->Integral();
    Double_t Nbkg2Max = h_HLT2->Integral();

    RooRealVar eff("eff","Efficiency",0.90,0.01,0.99);
    RooRealVar Nsig("Nsig","Signal Yield ",0.8*NsigMax, 0, NsigMax);
    RooRealVar c("c", "Correlation factor ", corr, 0.9, 1.1);

    RooFormulaVar Nsig0("Nsig0","max(0,1-2*eff+c*eff*eff)*Nsig",RooArgList(eff,Nsig,c));
    RooFormulaVar Nsig1("Nsig1","2*eff*max(0,1.0-c*eff)*Nsig",RooArgList(eff,Nsig,c));
    RooFormulaVar Nsig2("Nsig2","eff*eff*Nsig*c",RooArgList(eff,Nsig,c));

    RooRealVar Nbkg0("NbkgHLT0","Background count in HLT 0 sample",0.75*Nbkg0Max, 0.0, Nbkg0Max);
    RooRealVar Nbkg1("NbkgHLT1","Background count in HLT 1 sample",0.1*Nbkg1Max, 0.0, Nbkg1Max);
    RooRealVar Nbkg2("NbkgHLT2","Background count in HLT 2 sample",0.01*Nbkg2Max, 0.0, Nbkg2Max);

    // // turn out background for closure test
    // RooRealVar Nbkg0("NbkgHLT0","Background count in HLT 0 sample",0.0, 0.0, 0.0);
    // RooRealVar Nbkg1("NbkgHLT1","Background count in HLT 1 sample",0.0, 0.0, 0.0);
    // RooRealVar Nbkg2("NbkgHLT2","Background count in HLT 2 sample",0.0, 0.0, 0.0);

    RooAddPdf model0("model0","Model for HLT 0 sample",
        RooArgList(*(sig0->model),*(bkg0->model)),
        RooArgList(Nsig0, Nbkg0));
    RooAddPdf model1("model1","Model for HLT 1 sample",
        RooArgList(*(sig1->model),*(bkg1->model)),
        RooArgList(Nsig1, Nbkg1));
    RooAddPdf model2("model2","Model for HLT 2 sample",
        RooArgList(*(sig2->model),*(bkg2->model)),
        RooArgList(Nsig2, Nbkg2));

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(model0,"HLT0");
    totalPdf.addPdf(model1,"HLT1");
    totalPdf.addPdf(model2,"HLT2");

    TFile *fFit = new TFile(
        outputDir+"/workspace_"+etaRegion+"_"+iBin+".root",
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(*dataCombined);
    w->import(totalPdf);

    Int_t strategy = 2; // Minuit strategy
    RooFitResult *fitResult=0;
    RooFitResult *best_fitResult=0;
    Double_t best_chi2 = 100;
    Int_t best_fit = 0;

    int attempt = 0;  // fit attempt

    RooMsgService::instance().setSilentMode(kTRUE);

    do {
        std::cout<<">>> Fit with attempt "<<attempt<<std::endl;
        // attempt 0: just fit fail and pass region simultaneously in full range
        // attempt 1: first fit fail and pass pdf separately,
        //    then set proper initial start values
        //    then fit both together,
        // attempt 2: first fit fail pdf in sideband range,
        //    then in complete fail range,
        //    same for pass region
        //    then set proper initial start values
        //    then fit fail and pass region simultaneously in full range

        if(attempt>0){
            // reset parameters of fit models to initial values
            bkg0->Reset();
            bkg1->Reset();
            bkg2->Reset();
            sig0->Reset();
            sig1->Reset();
            sig2->Reset();

            c.setVal(1.015);
            eff.setVal(0.85);
            Nsig.setVal(NsigMax);

            if(attempt==3){
                c.setVal(1.03);
                eff.setVal(0.80);
                Nsig.setVal(NsigMax*0.5);
                // freeze efficiency and correlation so each region is fit seprarately first
                eff.setConstant(kTRUE);
                c.setConstant(kTRUE);

                m.setRange("rangeLow", massLo, 78);
                m.setRange("rangeHigh", 104, massHi);
                m.setRange("rangeCenter", 79, 103);
            }
            else{
                m.setRange("rangeLow", massLo, 81);
                m.setRange("rangeHigh", 101, massHi);
                m.setRange("rangeCenter", 81, 101);
            }

            if(attempt==2){
                std::cout<<">>> Fit sideband regions in fail"<<std::endl;

                fitResult = bkg0->model->fitTo(*dataHLT0,
                    RooFit::Range("rangeLow,rangeHigh"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());


                fitResult = sig0->model->fitTo(*dataHLT0,
                    RooFit::Range("rangeCenter"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());
            }
            std::cout<<">>> Fit full range pdf"<<std::endl;
            // fit total pdf in Fail
            fitResult = model0.fitTo(*dataHLT0,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                // RooFit::Minos(RooArgSet(Nsig)),
                RooFit::Save());

            const Double_t n0 = Nsig0.getVal();

            if(attempt==2){
                std::cout<<">>> Fit sideband regions in fail"<<std::endl;

                fitResult = bkg1->model->fitTo(*dataHLT1,
                    RooFit::Range("rangeLow,rangeHigh"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());

                fitResult = sig1->model->fitTo(*dataHLT1,
                    RooFit::Range("rangeCenter"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());
            }
            std::cout<<">>> Fit full range pdf"<<std::endl;
            // fit total pdf in Fail
            fitResult = model1.fitTo(*dataHLT1,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                // RooFit::Minos(RooArgSet(Nsig)),
                RooFit::Save());

            const Double_t n1 = Nsig1.getVal();

            if(attempt==2){
                std::cout<<">>> Fit sideband regions in fail"<<std::endl;

                fitResult = bkg2->model->fitTo(*dataHLT2,
                    RooFit::Range("rangeLow,rangeHigh"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());

                fitResult = sig2->model->fitTo(*dataHLT2,
                    RooFit::Range("rangeCenter"),
                    RooFit::PrintEvalErrors(-1),
                    RooFit::PrintLevel(-1),
                    RooFit::Warnings(0),
                    RooFit::Strategy(strategy), // MINUIT STRATEGY
                    // RooFit::Minos(RooArgSet(Nsig)),
                    RooFit::Save());
            }
            std::cout<<">>> Fit full range pdf"<<std::endl;
            // fit total pdf in Fail
            fitResult = model2.fitTo(*dataHLT2,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                // RooFit::Minos(RooArgSet(Nsig)),
                RooFit::Save());

            const Double_t n2 = Nsig2.getVal();

            const Double_t iNsig    = n0 + n1 + n2;
            const Double_t iEff     = (n1+2*n2)/(2*iNsig);
            const Double_t iC       = (4*n2*iNsig)/((n1+2*n2)*(n1+2*n2));

            std::cout<<">>> Set parameters to initial values:"<<std::endl;
            std::cout<<"Nsig = "<<iNsig<<std::endl;
            std::cout<<"eff = "<<iEff<<std::endl;
            std::cout<<"C = "<<iC<<std::endl;
            std::cout<<"<<< "<<std::endl;

            if(attempt==3){
                eff.setConstant(kFALSE);
                c.setConstant(kFALSE);
            }
            c.setVal(iC);
            eff.setVal(iEff);
            Nsig.setVal(iNsig);
        }
        // fit all regions together
        std::cout<<"--- fit all regions together -- "<<std::endl;
        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(0), // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(1), // MINUIT STRATEGY
            // RooFit::Minos(RooArgSet(Nsig, eff)),
            RooFit::Save());

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(strategy), // MINUIT STRATEGY
            // RooFit::Minos(kTRUE),
            RooFit::Save());

        RooChi2Var chi2Var("chi2", "chi 2", totalPdf, *dataCombined,
            RooFit::DataError(RooAbsData::Expected) //use Expected contribution from model PDF to calculate uncertainty
        );

        // reduced chi2 value from number of degree of freedom in the fit nDoF = nBins - nParams
        Double_t chi2ndf = chi2Var.getVal() / (dataCombined->numEntries() - fitResult->floatParsFinal().getSize());

        std::cout<<"---------------------------------------" <<std::endl;
        std::cout<<"------ attempt = " << attempt << std::endl;
        std::cout<<"------ chi2/ndf = " << chi2ndf <<std::endl;
        std::cout<<"---------------------------------------" <<std::endl;

        if(    std::abs(1 - (Nsig0.getVal()+Nbkg0.getVal())/h_HLT0->Integral()) > 0.1
            || std::abs(1 - (Nsig1.getVal()+Nbkg1.getVal())/h_HLT1->Integral()) > 0.1
            || std::abs(1 - (Nsig2.getVal()+Nbkg2.getVal())/h_HLT2->Integral()) > 0.1
        ){
            std::cout<<"WARNING: something went wrong in the fit, we give a bad chi2"<<std::endl;
            chi2ndf += 99;
        }

        if(attempt==0 || (chi2ndf < best_chi2 && chi2ndf > 0) ){
            best_chi2 = chi2ndf;
            best_fit = attempt;
            best_fitResult = fitResult;
            w->saveSnapshot(("snapshot_"+std::to_string(attempt)).c_str(),
                fitResult->floatParsFinal(), kTRUE);
        }

        attempt++;

    } while(attempt < 4); // && best_chi2 > 2)

    // load best fit values into workspace
    w->loadSnapshot(("snapshot_"+std::to_string(best_fit)).c_str());

    RooFormulaVar* NsigFormular = (RooFormulaVar*)w->arg("Nsig0");
    RooRealVar* NsigHLT0 = new RooRealVar("NsigHLT0","NsigHLT0",NsigFormular->getVal());
    NsigHLT0->setError(NsigFormular->getPropagatedError(*best_fitResult));
    NsigFormular = (RooFormulaVar*)w->arg("Nsig1");
    RooRealVar* NsigHLT1 = new RooRealVar("NsigHLT1","NsigHLT1",NsigFormular->getVal());
    NsigHLT1->setError(NsigFormular->getPropagatedError(*best_fitResult));
    NsigFormular = (RooFormulaVar*)w->arg("Nsig2");
    RooRealVar* NsigHLT2 = new RooRealVar("NsigHLT2","NsigHLT2",NsigFormular->getVal());
    NsigHLT2->setError(NsigFormular->getPropagatedError(*best_fitResult));    

    const Double_t chi2hlt0 = make_plot(nfl0, m, dataHLT0,
        sig0, bkg0, 
        h_HLT0->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("HLT0"),
        NsigHLT0, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgHLT0"),
        iBin, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("eff"), 
        (RooRealVar*)best_fitResult->floatParsFinal().find("c"), 
        "yield", etaRegion.Data(), 0);

    const Double_t chi2hlt1 = make_plot(nfl1, m, dataHLT1,
        sig1, bkg1, 
        h_HLT1->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("HLT1"),
        NsigHLT1, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgHLT1"),
        iBin, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("eff"), 
        (RooRealVar*)best_fitResult->floatParsFinal().find("c"), 
        "yield", etaRegion.Data(), 1);

    const Double_t chi2hlt2 = make_plot(nfl2, m, dataHLT2,
        sig2, bkg2, 
        h_HLT2->Integral(),
        ((RooSimultaneous*)w->pdf("totalPdf"))->getPdf("HLT2"),
        NsigHLT2, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("NbkgHLT2"),
        iBin, 
        (RooRealVar*)best_fitResult->floatParsFinal().find("eff"), 
        (RooRealVar*)best_fitResult->floatParsFinal().find("c"), 
        "yield", etaRegion.Data(), 2);

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2hlt0 = " << chi2hlt0 <<std::endl;
    std::cout<<"------ chi2hlt1 = " << chi2hlt1 <<std::endl;
    std::cout<<"------ chi2hlt2 = " << chi2hlt2 <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    RooRealVar chi2_0("chi2hlt0","chi2hlt0",chi2hlt0);
    RooRealVar chi2_1("chi2hlt1","chi2hlt1",chi2hlt1);
    RooRealVar chi2_2("chi2hlt2","chi2hlt2",chi2hlt2);
    RooRealVar chi2("chi2","chi2",best_chi2);

    w->import(chi2_0);
    w->import(chi2_1);
    w->import(chi2_2);
    w->import(chi2);
    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    delete dataCombined;
    delete dataHLT0;
    delete dataHLT1;
    delete dataHLT2;
    delete sig0;
    delete bkg0;
    delete sig1;
    delete bkg1;
    delete sig2;
    delete bkg2;
}


//--------------------------------------------------------------------------------------------------
double extractCorrelation_HLT(const TString mcfilename, TH1D *hPV, const TString etaRegion){
    
    if(mcfilename == "")
        return 1.0;
    
    TFile* histfile = generateTemplate_cHLT(mcfilename);
    
    TH1D* h0 = (TH1D*)histfile->Get("h_npv_0hlt_"+etaRegion);        
    TH1D* h1 = (TH1D*)histfile->Get("h_npv_1hlt_"+etaRegion);    
    TH1D* h2 = (TH1D*)histfile->Get("h_npv_2hlt_"+etaRegion);
    TH1D* hPV_MC = (TH1D*)histfile->Get("hPV");
    
    if(hPV==0){
        std::cout<<"WARNING: data pileup histogram is 0, impossible to calculate HLT correlation! return 1."<<std::endl;
        return 1.;
    }
    
    TH1D* hRatio = (TH1D*)hPV->Clone("hPV_ratio");
    hRatio->Divide(hPV_MC);

    double n0 = 0;
    double n1 = 0;
    double n2 = 0;
    
    for(int i=0; i<= hPV_MC->GetNbinsX()+1; i++){
        
        const double wgt = hRatio->GetBinContent(i);
        
        n0 += wgt*h0->GetBinContent(i);
        n1 += wgt*h1->GetBinContent(i);
        n2 += wgt*h2->GetBinContent(i);
    }

    if(n1+n2 == 0){
        std::cout<<"WARNING: no HLT events found in data, impossible to calculate HLT correlation! return 1"<<std::endl;
        return 1.;
    }

    const double n = n0+n1+n2;
    
    const double corr = 4*n*n2 / std::pow(n1 + 2*n2, 2);
    
    std::cout<<"EXTRACT c("<<etaRegion<<") = "<<corr<<std::endl;
    
    return corr;
    
}


//--------------------------------------------------------------------------------------------------
std::vector<double> preFit(TH1* failHist){
    std::cout<<"PREFIT"<<std::endl;
    // do not fit between 81 to 101 GeV

    //  std::vector<float> v = {1.,0.,1.,0.,1.,0.};return v;
    TH1D *h = new TH1D("h", "", massBin, massLo, massHi);
    TF1 *fq = new TF1("fq", "[0]+[1]*x+[2]*x*x", massLo, massHi);   // quadratic
    TF1 *fl = new TF1("fl", "[0]+[1]*x", massLo, massHi);   // linear
    
    // define sideband region
    for(int i = 0; i < (int)((81-massLo)/massWidth); i++){
        h->SetBinContent(i+1, failHist->GetBinContent(i+1));
        h->SetBinError(i+1, failHist->GetBinError(i+1));
    }
    for(int i = (int)((massHi-101)/massWidth); i < massBin; i++){
        h->SetBinContent(i+1, failHist->GetBinContent(i+1));
        h->SetBinError(i+1, failHist->GetBinError(i+1));
    }

    h->Fit("fq");
    std::cout<<fq->GetParameter(0)<<","<<fq->GetParError(0)<<" ,"<<fq->GetParameter(1)<<","<<fq->GetParError(1)<<","<<fq->GetParameter(2)<<","<<fq->GetParError(2)<<std::endl;
    std::vector<double> v;

    if(fq->GetParameter(2) > 0.){
        h->Fit("fl");
        std::cout<<"SWITCH TO LINEAR!"<<std::endl;
        std::cout<<fl->GetParameter(0)<<","<<fl->GetParError(0)<<" ,"<<fl->GetParameter(1)<<","<<fl->GetParError(1)<<std::endl;
        v = {fl->GetParameter(0), fl->GetParError(0), fl->GetParameter(1), fl->GetParError(1), 0.,0.};
    }
    else{
        v = {fq->GetParameter(0), fq->GetParError(0), fq->GetParameter(1), fq->GetParError(1), fq->GetParameter(2), fq->GetParError(2)};
    }

    return v;
}
