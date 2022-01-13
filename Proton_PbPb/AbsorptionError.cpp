// Spectra.cpp
// This macro computes the absorption systematic uncertainty

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
using namespace proton;

const int N_TRIALS = 10000;

double uncertaintyPt(int iMatt, double pt){
  if (iMatt == 1)
    return 0.00557*TMath::Power(pt,-0.32839);
  return 0.01788*TMath::Power(pt,-0.58252);
  //return 0;
};

void AbsorptionError(const char *outFileName = "AbsError", const char *outFileOption = "recreate", const char *ratioFile = "SpectraProton")
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(44);

  TFile *inFileRatio = TFile::Open(Form("%s/%s.root", kOutDir, ratioFile));
  if (!inFileRatio)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  TH1D *fRatio[kNCentClasses];
  TH1D *fAbsErrorDist[kNCentClasses];
  for (int iCent = 0; iCent < kNCentClasses; ++iCent)
  {
    // std::cout << "read: " << Form("%s_%d_%d/fATOFrawYield_%.0f_%.0f", cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]) << std::endl;
    fRatio[iCent] = new TH1D(*(TH1D *)inFileRatio->Get(Form("fASpectra_%.0f_%.0f",kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1])));
    fAbsErrorDist[iCent] = new TH1D(Form("fFitPar_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), 2500, 0., 5.);

    TH1D *fSpectra[2];
    for (int iTrial = 0; iTrial < N_TRIALS; ++iTrial){
      fRatio[iCent]->Clear();
      // compute corrected spectra
      for (int iMatt = 0; iMatt < 2; ++iMatt)
      {
        outFile.cd();
        double cent_bin_lim_min = kCentBinsLimitsProton[iCent][0], cent_bin_lim_max = kCentBinsLimitsProton[iCent][1];

        // extract nsigma from standard gaussian
        double nsigma = gRandom->Gaus();

        //sec->Fit(&fitFuncSec,"R");
        fSpectra[iMatt] = (TH1D *)inFileRatio->Get(Form("f%sSpectra_%.0f_%.0f",kAntimatterMatter[iMatt],kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
        int pTbinMax = 24;
        // std::cout<<"entering pt loop..."<<std::endl;
        for (int iPtBin = 5; iPtBin < pTbinMax + 1; ++iPtBin)
        {
          double pt = fSpectra[iMatt]->GetBinCenter(iPtBin);
          double spec_tmp = fSpectra[iMatt]->GetBinContent(iPtBin);
          double correction = 1+nsigma*uncertaintyPt(iMatt,pt);
          //if (iTrial % 100 == 0){std::cout << "Correction: " << correction <<std::endl;}
          fSpectra[iMatt]->SetBinContent(iPtBin, spec_tmp*correction);
          fSpectra[iMatt]->SetBinError(iPtBin, fSpectra[iMatt]->GetBinError(iPtBin)*correction);
        }
      }
      if (iTrial % 100 == 0){
        //fSpectra[0]->Write();
        //fSpectra[1]->Write();
      }

      // compute ratios
      int pTbinMax = 24;
      for (int iPtBin = 5; iPtBin < pTbinMax + 1; ++iPtBin)
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
        fAbsErrorDist[iCent]->Fill(fRatio[iCent]->GetFunction("pol0")->GetParameter(0));
      }
      fSpectra[0]->Clear();
      fSpectra[1]->Clear();
    }
    fAbsErrorDist[iCent]->Write();
  }
  outFile.Close();
}