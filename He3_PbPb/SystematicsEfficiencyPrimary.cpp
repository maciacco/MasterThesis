// Spectra.cpp
// This macro computes the systematic uncertainty due to efficiency and TFF

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using utils::TTList;
using namespace he3;

const int N_TRIALS = 10000;

void SystematicsEfficiencyPrimary(const char *outFileName = "SystematicsEfficiencyPrimary_try", const char *outFileOption = "recreate", const char *ratioFile = "SpectraHe3_kINT7", const char *efficiencyFile = "EfficiencyHe3_kINT7", const char *primaryFile = "PrimaryHe3_kINT7")
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(44);

  TFile *inFileRatio = TFile::Open(Form("%s/%s.root", kOutDir, ratioFile));
  TFile *inFileEff = TFile::Open(Form("%s/%s.root", kOutDir, efficiencyFile));
  TFile *inFilePrim = TFile::Open(Form("%s/%s.root", kOutDir, primaryFile));
  if (!inFileRatio || !inFileEff || !inFilePrim)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  TH1D *fRatio[kNCentClasses];
  TH1D *fRatioDistribution[kNCentClasses];
  TH1D *fEff;
  TH1D *fPrim;

  for (int iCent = 0; iCent < kNCentClasses; ++iCent)
  {
    fRatioDistribution[iCent] = new TH1D(Form("fRatioDistribution_%.0f_%.0f",kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]),";R(^{3}#bar{He}/^{3}He);Entries",400,0.8,1.2);
    // std::cout << "read: " << Form("%s_%d_%d/fATOFrawYield_%.0f_%.0f", cutSettings, binCounting, bkg_shape, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]) << std::endl;
    fRatio[iCent] = new TH1D(*(TH1D *)inFileRatio->Get(Form("1.0_89_0.1_2.5_1_1_1/fASpectra_%.0f_%.0f",kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1])));

    TH1D *fSpectra[2];
    for (int iTrial = 0; iTrial < N_TRIALS; ++iTrial){
      fRatio[iCent]->Clear();
      // compute corrected spectra
      for (int iMatt = 0; iMatt < 2; ++iMatt)
      {
        outFile.cd();
        double cent_bin_lim_min = kCentBinsLimitsHe3[iCent][0], cent_bin_lim_max = kCentBinsLimitsHe3[iCent][1];

        // extract nsigma from standard gaussian
        double nsigma = gRandom->Gaus();

        //sec->Fit(&fitFuncSec,"R");
        fSpectra[iMatt] = (TH1D *)inFileRatio->Get(Form("1.0_89_0.1_2.5_1_1_1/f%sSpectra_%.0f_%.0f",kAntimatterMatter[iMatt],kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
        fEff = (TH1D*)inFileEff->Get(Form("f%sEff_TPC_%.0f_%.0f",kAntimatterMatter[iMatt],kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
        fPrim = (TH1D *)inFilePrim->Get(Form("f%sPrimFrac_%.0f_%.0f",kAntimatterMatter[iMatt],kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    
        int pTbinMax = 11;
        if (iCent < kNCentClasses - 1)
          pTbinMax = 13;
        // std::cout<<"entering pt loop..."<<std::endl;
        for (int iPtBin = 3; iPtBin < pTbinMax + 1; ++iPtBin)
        {
          double spec_tmp = fSpectra[iMatt]->GetBinContent(iPtBin);
          double eff_tmp = fEff->GetBinContent(iPtBin);
          double eff_tmp_err = fEff->GetBinError(iPtBin);
          double prim_tmp = fPrim->GetBinContent(iPtBin);
          double prim_tmp_err = fPrim->GetBinError(iPtBin);
          double tmp_shift = gRandom->Gaus(1,TMath::Sqrt(eff_tmp_err*eff_tmp_err/eff_tmp/eff_tmp/* +prim_tmp_err*prim_tmp_err/prim_tmp/prim_tmp */));
          fSpectra[iMatt]->SetBinContent(iPtBin, spec_tmp*tmp_shift);
          fSpectra[iMatt]->SetBinError(iPtBin, fSpectra[iMatt]->GetBinError(iPtBin)*tmp_shift);
        }
      }
      if (iTrial % 100 == 0){
        //fSpectra[0]->Write();
        //fSpectra[1]->Write();
      }

      // compute ratios
      int pTbinMax = 11;
      if (iCent < kNCentClasses - 1)
        pTbinMax = 13;
      for (int iPtBin = 3; iPtBin < pTbinMax + 1; ++iPtBin)
      {
        double antiSpec = fSpectra[0]->GetBinContent(iPtBin);
        double spec = fSpectra[1]->GetBinContent(iPtBin);
        double antiSpecErr = fSpectra[0]->GetBinError(iPtBin);
        double specErr = fSpectra[1]->GetBinError(iPtBin);
        if (spec > 1.e-8 && antiSpec > 1.e-8)
        {
          fRatio[iCent]->SetBinContent(iPtBin, antiSpec / spec);
          fRatio[iCent]->SetBinError(iPtBin, antiSpec / spec * TMath::Sqrt(antiSpecErr * antiSpecErr / antiSpec / antiSpec + specErr * specErr / spec / spec));
        }
      }
      fRatio[iCent]->Fit("pol0","Q");
      double chi2 = fRatio[iCent]->GetFunction("pol0")->GetChisquare();
      double ndf = fRatio[iCent]->GetFunction("pol0")->GetNDF();
      if (chi2/ndf < 10.){
        fRatioDistribution[iCent]->Fill(fRatio[iCent]->GetFunction("pol0")->GetParameter(0));
      }
      fSpectra[0]->Clear();
      fSpectra[1]->Clear();
    }
    fRatioDistribution[iCent]->Write();
  }
  outFile.Close();
}