//------------------------------------------------------------------------------
//
// BLUE: A ROOT class implementing the Best Linear Unbiased Estimate method.
//
// Copyright (C) 2012-2019, Richard.Nisius@mpp.mpg.de
// All rights reserved
//
// This file is part of BLUE - Version 2.2.0.
//
// BLUE is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// For the licensing terms see the file COPYING or http://www.gnu.org/licenses.
//
//------------------------------------------------------------------------------
#include "TROOT.h"
#include "Blue.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include <iostream>                 // standard I/O cout/cin
#include <math.h>                   // sqrt


void LUMI_COMBINATION_2018_LEGACY_LUMI(Int_t Flag = 0){

    gROOT->SetBatch(1);

  //----------------------------------------------------------------------------
  // 
  //
  //
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  2;
  static const Int_t NumUnc =  2;
  static const Int_t MaxObs =  1;
  Int_t NumObs =  1;

  // The names
  TString NamEst[NumEst] = {"   P_18", "   Z_18"};
  TString NamUnc[NumUnc] = {" Uncorr P", "  Corr P1"};
  TString NamObs[MaxObs] = {"L(2018)"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {0, 0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... Combination of 2018 PHYSICS luminosity ");
    printf(" with lumi from low PU(2017) and Z counts for 2016, 2017 and 2018");
    printf(" Flag = %2i \n", Flag);
    printf(" Blind combination \n");
  }else if(Flag == 1){
    printf("... Combination of 2018 PHYSICS luminosity with Z counts for 2017 and 2018");
    printf(" Flag = %2i \n", Flag);
    printf(" Unblind combination \n");
  }else{
    printf("... : Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Character to set the file names
  char Buffer[100];

  // Estimates 0-6
  // First the F0

  // set some values

  // measured physics luminosities
  // Double_t lP16 = 36.31;
  // Double_t lP17 = 41.48;
  Double_t lP18 = 59.83;

  // using updated HF normtags
  // Double_t lZ16 = lP16 * 0.978;
  // Double_t lZ17 = lP17 * 0.995;
  Double_t lZ18 = lP18 * 0.985;

  if(Flag == 0){
    // lZ16 = lP16;
    // lZ17 = lP17;
    lZ18 = lP18;
  }

  // Some toy value for the ratio of Z count such that equal to ratio
  // Correlated relative uncertainty on Z counts
  // (uncorrelated) luminosity uncertainty of low PU sample: 0.2%
  // limited statistics of low PU sample: 0.27%
  // experimental uncertainty on ratio: 0.5%
  // --> 0.6%

  // uncertainties in absolute numbers
  // 0 Uncorrelated lumi
  // 1 full Correlated PHYSICS lumi
  // 2 correlated PHYSICS lumi 17/18
  // 3 correlated PHYSICS lumi 17/17H
  // 4-7 correlated Z
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    // Cor:    0             1         
    //         Uncorr        Corr1    
    lP18,      0.0051*lP18,  0.0071*lP18, 
    lZ18,      0.0096*lZ18,  0.0051*lZ18,
  };

  static const Int_t LenCor = NumEst * NumEst;

  // Uncertainty source 1: Uncorrelated PHYSICS luminosity
  Double_t Cor00[LenCor] = {
  // P(17) P(18) L(17) L(18)
    +1.00, 0.00, // P(18)
    +0.00, 1.00  // L(18)
  };
  // Uncertainty source 2: correlated PHYSICS luminosity 16,17,18,17H -> correlated with Z luminosity
  Double_t Cor01[LenCor] = {
      +1.00, 1.00,
      +1.00, 1.00
  };

  //-- Local Structures for Blue output
  // TMatrices
  TMatrixD* LocRho    = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRhoRes = new TMatrixD(NumObs,NumObs);
  TMatrixD* LocWeight = new TMatrixD(NumEst,NumObs);
  //-- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.3f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%5.3f";
  const TString ForRho = ForWei;
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 0){
      myBlue->FillCor(k,&Cor00[0]);
    }
    else{
      myBlue->FillCor(k,&Cor01[0]);
    }
  }

    printf("... LUMI_COMBINATION_3: Fix the input \n");

    myBlue->FixInp();

    printf("... LUMI_COMBINATION_3: Print estimators \n");
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->PrintEst(i);
    }
    sprintf(Buffer,"LUMI_COMBINATION_3 %i",Flag);
    myBlue->PrintCompatEst(Buffer);
    printf("... LUMI_COMBINATION_3: The correlations");
    printf(" of the estimates in %%\n");
    myBlue->GetRho(LocRho);
    LocRho->operator*=(100);
    myBlue->PrintMatrix(LocRho,"%+4.0f");

    myBlue->Solve();
    printf("... LUMI_COMBINATION_3: The Luminosity 17/18 Combination");
    printf(" Flag = %2i.\n",Flag);
    myBlue->PrintResult();
    // myBlue->LatexResult("LUMI_COMBINATION_3");

    printf("... LUMI_COMBINATION_3: The correlations");
    printf(" of the observables in %%\n");
    myBlue->GetRhoRes(LocRhoRes);
    LocRhoRes->operator*=(100);
    myBlue->PrintMatrix(LocRhoRes,"%+4.0f");

    printf("... LUMI_COMBINATION_3: A difference is observed in the weights");
    printf(" that is under discussion with the authors \n");
    printf("... LUMI_COMBINATION_3: The observed weights\n");
    myBlue->GetWeight(LocWeight);
    myBlue->PrintMatrix(LocWeight," %+5.3f");
    // myBlue->PrintCompatObs();

    // --- Make combination plot

    // initil estimates
    Double_t* RetEstVal = new Double_t(NumEst);
    myBlue->GetEstVal(RetEstVal);


    // initil uncertainties
    Double_t* RetEstUnc = new Double_t(NumEst);
    myBlue->GetEstUnc(RetEstUnc);


    // result
    static const Int_t LenXRes = NumObs * (NumUnc+1);
    Double_t* RetResult = new Double_t(LenXRes);
    myBlue->GetResult(RetResult);

    // total uncertainty on result
    Double_t* RetUncert = new Double_t(NumObs);
    myBlue->GetUncert(RetUncert);

    // compute the sum and uncertainties of all
    double RetEstVal_PHYSICS[2];
    double RetEstVal_Z[2];
    double RetEstUnc_PHYSICS[2];
    double RetEstUnc_Z[2];
    double RetEstVal_COMB[2];
    double RetEstUnc_COMB[2];

    RetEstVal_COMB[1] = 0.;
    RetEstVal_PHYSICS[1] = 0.;
    RetEstVal_Z[1] = 0.;
    RetEstUnc_COMB[1] = 0.;
    RetEstUnc_PHYSICS[1] = 0.;
    RetEstUnc_Z[1] = 0.;

    // propagation of uncertainties
    for(int i=0; i < 1; i++){
        for(int j=0; j < 1; j++){
            RetEstUnc_PHYSICS[1] += (*LocRho)[i][j]/100. * RetEstUnc[i] * RetEstUnc[j];
            RetEstUnc_Z[1] += (*LocRho)[i+1][j+1]/100. * RetEstUnc[i+1] * RetEstUnc[j+1];
            RetEstUnc_COMB[1] += (*LocRhoRes)[i][j]/100. * RetUncert[i] * RetUncert[j];
        }
        RetEstVal_PHYSICS[1] += RetEstVal[i];
        RetEstVal_Z[1] += RetEstVal[i+1];
        RetEstVal_COMB[1] += RetResult[i * (NumUnc+1)];
    }

    RetEstUnc_PHYSICS[1] = std::sqrt(RetEstUnc_PHYSICS[1]);
    RetEstUnc_Z[1]       = std::sqrt(RetEstUnc_Z[1]);
    RetEstUnc_COMB[1]    = std::sqrt(RetEstUnc_COMB[1]);

    // put in values and normalize to luminosity COMBINATION

    RetEstVal_PHYSICS[0] = RetEstVal[0] / RetResult[0];
    RetEstVal_Z[0] = RetEstVal[3] / RetResult[0];

    RetEstUnc_PHYSICS[0] = RetEstUnc[0] / RetResult[0];
    RetEstUnc_Z[0] = RetEstUnc[3] / RetResult[0];

    RetEstVal_COMB[0] = RetResult[0] / RetResult[0];
    RetEstUnc_COMB[0] = RetUncert[0] / RetResult[0];

    std::cout<<"--- RunII ---"<<std::endl;
    std::cout<<"PHYSICS    : "<< RetEstVal_PHYSICS[1] <<"+/-"<<RetEstUnc_PHYSICS[1] <<std::endl;
    std::cout<<"Z          : "<< RetEstVal_Z[1] <<"+/-"<<RetEstUnc_Z[1] <<std::endl;
    std::cout<<"Combination: "<< RetEstVal_COMB[1] <<"+/-"<<RetEstUnc_COMB[1] <<std::endl;

    // and sums
    RetEstUnc_PHYSICS[1] /= RetEstVal_COMB[1];
    RetEstUnc_Z[1]       /= RetEstVal_COMB[1];
    RetEstUnc_COMB[1]    /= RetEstVal_COMB[1];
    RetEstVal_Z[1]       /= RetEstVal_COMB[1];
    RetEstVal_PHYSICS[1] /= RetEstVal_COMB[1];
    RetEstVal_COMB[1]    /= RetEstVal_COMB[1];

    std::cout<<"--- 2018 (normalized) ---"<<std::endl;
    std::cout<<"PHYSICS    : "<< RetEstVal_PHYSICS[0] <<"+/-"<<RetEstUnc_PHYSICS[0] <<std::endl;
    std::cout<<"Z          : "<< RetEstVal_Z[0]       <<"+/-"<<RetEstUnc_Z[0]       <<std::endl;
    std::cout<<"Combination: "<< RetEstVal_COMB[0]    <<"+/-"<<RetEstUnc_COMB[0]    <<std::endl;

    // // put into graphs
    // double y[]  = {3., 2., 1., 0.};
    // double ey[] = {0.2, 0.2, 0.2, 0.2};
    // double ey2[] = {0., 0., 0., 0.};
    // std::cout<<"--- Set Graphs >>>"<<std::endl;
    // TGraphErrors *gPHYSICS = new TGraphErrors(4, RetEstVal_PHYSICS, y, RetEstUnc_PHYSICS, ey);
    // TGraphErrors *gZ       = new TGraphErrors(4, RetEstVal_Z,       y, RetEstUnc_Z,       ey2);
    // TGraphErrors *gCOMB    = new TGraphErrors(4, RetEstVal_COMB,    y, RetEstUnc_COMB,    ey);
    // std::cout<<"--- Set Graphs <<<"<<std::endl;
    //
    // Double_t margin_left = 0.2;
    // Double_t margin_right = 0.01;
    // Double_t margin_top = 0.06;
    // Double_t margin_bottom = 0.15;
    //
    // std::cout<<"--- Set Canvas ---"<<std::endl;
    // TCanvas *canvas = new TCanvas("canvasGCC", "", 600, 600);
    // TPad *pad1 = new TPad("pad1GCC", "pad1GCC", 0., 0.0, 1, 1.0);
    // canvas->SetTicks();
    // pad1->SetLeftMargin(margin_left);
    // pad1->SetRightMargin(margin_right);
    // pad1->SetTopMargin(margin_top);
    // pad1->SetBottomMargin(margin_bottom);
    // pad1->SetTickx();
    // pad1->SetTicky();
    // pad1->Draw();
    // pad1->cd();
    //
    //
    //
    // Double_t textsize1 = 28./(pad1->GetWh()*pad1->GetAbsHNDC());
    //
    // gCOMB->SetTitle(" ");
    // TAxis *yAxis = gCOMB->GetYaxis();
    // yAxis->SetLabelSize(textsize1);
    // yAxis->SetTitleFont(42);
    // yAxis->SetTitleSize(textsize1*1.2);
    // yAxis->SetTitle("Global correlation coeff.");
    // TAxis *xAxis = gCOMB->GetXaxis();
    // xAxis->SetLabelSize(textsize1);
    // xAxis->SetTitleFont(42);
    // xAxis->SetTitleSize(textsize1*1.2);
    // xAxis->SetTitle("Relative uncertainty");
    //
    // std::cout<<"--- Set Ticks ---"<<std::endl;
    // // yAxis->SetNdivisions(-504);
    // // yAxis->ChangeLabel(0,-1,-1,-1,-1,-1,"RunII");
    // // yAxis->ChangeLabel(1,-1,-1,-1,-1,-1,"2018");
    // // yAxis->ChangeLabel(2,-1,-1,-1,-1,-1,"2017");
    // // yAxis->ChangeLabel(3,-1,-1,-1,-1,-1,"2016");
    // // for(int i=0;i<4;i++){
    // //     yAxis->ChangeLabel(1,-1,-1,-1,-1,-1,year[i]);
    // // }
    //
    // // yWidth = globalMax - globalMin
    // // gGccScan_avg.SetMinimum(globalMin-yWidth*0.01)
    // // gGccScan_avg.SetMaximum(globalMax+yWidth*0.2)
    // gCOMB->SetMinimum(-0.5);
    // gCOMB->SetMaximum(3.5);
    // gCOMB->SetFillColor(2);
    // gCOMB->SetFillStyle(3001);
    //
    // gPHYSICS->SetMarkerStyle(20);
    // gPHYSICS->SetMarkerColor(1);
    // gPHYSICS->SetMarkerSize(2);
    //
    // gZ->SetLineColor(1);
    // gZ->SetFillColor(0);
    //
    // std::cout<<"--- Draw 1 ---"<<std::endl;
    // gCOMB->Draw("A2");
    // // std::cout<<"--- Draw 2 ---"<<std::endl;
    // // gPHYSICS->Draw("P same");
    // // std::cout<<"--- Draw 3 ---"<<std::endl;
    // // gZ->Draw("2 same");
    // // std::cout<<"--- Draw 4 ---"<<std::endl;
    // // gCOMB->Draw("A same");
    //
    // // cms(x=xCMS, y=yCMS, textsize=textsize1)
    // //
    // // ltex2 = ROOT.TLatex()
    // // ltex2.SetNDC()
    // // ltex2.SetTextAlign(11)
    // // ltex2.SetTextFont(42)
    // // ltex2.SetTextSize(textsize1)
    // // ltex2.DrawLatex(margin_left, 1.02-margin_top, observables[_name])
    // //
    // // legend = ROOT.TLegend(0.3 ,0.68, 0.8, 0.9)
    // // legend.SetTextSize(textsize1*0.8)
    // //
    // // legend.AddEntry(gGccScan_avg,"#rho_{avg} (Nominal fit)","p")
    // // legend.AddEntry(gGccScan_max,"#rho_{max} (Nominal fit)","p")
    // // legend.AddEntry(gGccScanStat_avg,"#rho_{avg} (Stat. only fit)","p")
    // // legend.AddEntry(gGccScanStat_max,"#rho_{max} (Stat. only fit)","p")
    // //
    // // legend.Draw()
    // std::cout<<"--- Save ---"<<std::endl;
    //
    // // canvas->Draw();
    //
    // canvas->SaveAs("Combination.png");
    // // canvas->SaveAs("Combination.eps");
    //
    // std::cout<<"--- Close ---"<<std::endl;
    // canvas->Close();





  // Delete Object
  std::cout<<"--- Delete ---"<<std::endl;
  delete myBlue; myBlue = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;
  LocWeight->Delete(); LocWeight = NULL;
  // gPHYSICS->Delete(); gPHYSICS = NULL;
  // gZ->Delete(); gZ = NULL;
  // gCOMB->Delete(); gCOMB = NULL;
  // canvas->Delete(); canvas = NULL;
  std::cout<<"--- Return ---"<<std::endl;
  return;
}
