// Systematics.cpp
// This macro computes systematic errors on pol0 fit parameter

#include <stdlib.h>
#include <string.h>

#include <TFile.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <Riostream.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveStats.h>

#include "../utils/Config.h"

using namespace triton;

// #define USE_COUNTER


double he3CorrectionPt(int iMatt, double pt){
  if (iMatt == 1)
    return 1.010*TMath::Power(pt, -0.00000239);
  return 1.034*TMath::Power(pt, -0.007623);
};

double tCorrectionPt(int iMatt, double pt){
  if (iMatt == 1)
    return 1.03776*TMath::Power(pt, -0.00576917);
  return 1.03887*TMath::Power(pt, -0.00675656);
}

void Systematics(const int points = kNPoints, const bool cutVar = true, const bool binCountingVar = true, const bool expVar = true, const bool sigmoidVar = false, const char *outFileName = "SystematicsAll")
{
  gStyle->SetOptStat(110001110);
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.85);
  TStopwatch swatch;
  swatch.Start(true);

  TFile *specFile = TFile::Open(Form("%s/SpectraHe3Syst.root", kOutDir));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");
  TDirectory *cdHist = outFile->mkdir("hist");

  for (int iC = 0; iC < kNCentClasses; ++iC)
  {
    TDirectory *cdFits = outFile->mkdir(Form("fits_%.0f_%.0f", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]));
    TH1D fFitPar(Form("fFitPar_%.0f_%.0f", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), 2000, 0.5, 1.5);
    TH1D fProb(Form("fProb_%.0f_%.0f", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), 1000., 0., 1.);

    double pTbins[kNPtBins + 1] = {1.4f,1.6f,2.f,2.4f,3.f};

    TH1D fRatio("fRatio", "fRatio", kNPtBins, pTbins);
    int nUsedPtBins = 4;

    // fill and fit ratio histogram nPoint times
    int iP = 0;
    while (iP < kNPoints)
    {
      for (int iPtBin = 2; iPtBin < nUsedPtBins + 1; ++iPtBin)
      {
        char hname[100];
        double cutDCAzRnd = 1.;
        if (cutVar)
          cutDCAzRnd = (int)(gRandom->Rndm() * (kNCutDCAz)) * 0.2 + kCutDCAz[0];
        double cutTPCclsRnd = 69;
        if (cutVar)
          cutTPCclsRnd = (int)(gRandom->Rndm() * (kNCutTPCClusters)) * 2 + kCutTPCClusters[0];
        double cutDCAxyRnd = 0.1;
        if (cutVar){
          cutDCAxyRnd = (int)(gRandom->Rndm() * (kNCutDCAxy)) * 0.01 + kCutDCAxy[0];
        }
        double cutChi2TPCRnd = 2.50;
        if (cutVar)
          cutChi2TPCRnd = (int)(gRandom->Rndm() * (kNCutChi2TPC)) * 0.50 + kCutChi2TPC[0];

        int binCountingFlagRnd = 1;
        if (binCountingVar)
          (gRandom->Rndm() > .5) ? binCountingFlagRnd = 1 : binCountingFlagRnd = 0;

        // EXPO AND POL1 SAMPLES ARE MERGED
        int expFlag = 0;
        // if (expVar)
        //   expFlag = 0;
        (gRandom->Rndm() > .5) ? expFlag = 1 : expFlag = 0;

        int sigmoidFlagRnd = 0;
        if (sigmoidVar)
        {
          if (fRatio.GetBinCenter(iPtBin) < 6.3)
          {
            (gRandom->Rndm() > .5) ? sigmoidFlagRnd = 1 : sigmoidFlagRnd = 0;
          }
        }

        if (cutTPCclsRnd < 99.5)
          snprintf(hname, 24, "%1.1f_0%.0f_%1.2f_%1.2f_%d_%d_%d", cutDCAzRnd, cutTPCclsRnd, cutDCAxyRnd, cutChi2TPCRnd, binCountingFlagRnd, expFlag, sigmoidFlagRnd);
        else
          snprintf(hname, 24, "%1.1f_%.0f_%1.2f_%1.2f_%d_%d_%d", cutDCAzRnd, cutTPCclsRnd, cutDCAxyRnd, cutChi2TPCRnd, binCountingFlagRnd, expFlag, sigmoidFlagRnd);

        std::cout<<"hname="<<hname<<std::endl;
        TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", hname, kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]));
        if (!h) continue;
        double correction = he3CorrectionPt(true,fRatio.GetBinCenter(iPtBin))/he3CorrectionPt(false,fRatio.GetBinCenter(iPtBin))*tCorrectionPt(false,fRatio.GetBinCenter(iPtBin))*tCorrectionPt(true,fRatio.GetBinCenter(iPtBin));
        fRatio.SetBinContent(iPtBin, h->GetBinContent(iPtBin)*correction);
        fRatio.SetBinError(iPtBin, h->GetBinError(iPtBin)*correction);
      }
      std::cout << "p=" << iP << std::endl;
      TF1 fitFunc("fitFunc", "pol0");
      auto fit = fRatio.Fit(&fitFunc, "QS");

      int ndf = 2;
      if (fit->Status() == 0 && fit->Prob()>0.025 && fit->Prob()<0.975 && fit->Ndf() == ndf) //fit->Prob()>0.9 &&
      { // check chi2
        //if (fitFunc.GetParameter(0)<1.08) continue;
        fFitPar.Fill(fitFunc.GetParameter(0));
        fProb.Fill(fitFunc.GetProb());
        cdFits->cd();
        // fRatio.Write();
        fRatio.Reset();
        ++iP;
      }
    }

    cdHist->cd();
    fFitPar.GetYaxis()->SetTitle("Entries");
    fFitPar.GetXaxis()->SetTitle("R (^{3}#bar{H}/^{3}H)");
    fFitPar.SetDrawOption("histo");
    fFitPar.Rebin(8);
    fFitPar.SetFillStyle(3345);
    fFitPar.SetLineWidth(2);
    fFitPar.SetLineColor(kBlue);
    fFitPar.SetFillColor(kBlue);
    fFitPar.Write();

    TCanvas cFitPar("cFitPar", "cFitPar");
    cFitPar.cd();
    if (iC<2)fFitPar.Rebin(2);
    cFitPar.SetTicks(1, 1);
    fFitPar.GetXaxis()->SetRangeUser(fFitPar.GetMean()-5*fFitPar.GetRMS(),fFitPar.GetMean()+5*fFitPar.GetRMS());
    fFitPar.Draw("");
    cFitPar.Print(Form("%s/Systematics_%s.pdf", kPlotDir, fFitPar.GetName()));

    fProb.GetYaxis()->SetTitle("Entries");
    fProb.GetXaxis()->SetTitle("#it{P}( #chi^{2} > #chi^{2}_{#it{obs}} )");
    fProb.Write();
  }
  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;
}
