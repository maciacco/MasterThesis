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

#include "../utils/Config.h"

using namespace he3;

// #define USE_COUNTER

void Systematics(const int points = kNPoints, const bool cutVar = true, const bool binCountingVar = true, const bool expVar = true, const bool sigmoidVar = true, const char *outFileName = "SystematicsAll")
{
  gStyle->SetOptStat(110001110);
  TStopwatch swatch;
  swatch.Start(true);

  TFile *specFile = TFile::Open(Form("%s/SpectraHe3Syst.root", kOutDir));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");
  TDirectory *cdHist = outFile->mkdir("hist");

  for (int iC = 0; iC < kNCentClasses; ++iC)
  {
    TDirectory *cdFits = outFile->mkdir(Form("fits_%.0f_%.0f", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]));
    TH1D fFitPar(Form("fFitPar_%.0f_%.0f", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), 2000, 0.8, 1.0);
    TH1D fProb(Form("fProb_%.0f_%.0f", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]), 1000., 0., 1.0);

    double pTbins[kNPtBins + 1] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 6.5f, 7.f, 8.f, 10.f};

    TH1D fRatio("fRatio", "fRatio", kNPtBins, pTbins);
    int nUsedPtBins = 13;
    if (iC == kNCentClasses - 1)
    {
      int nPtBins = 13;
      double pTbinsNew[] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 7.f, 8.f, 10.f};
      fRatio.SetBins(nPtBins, pTbinsNew);
      nUsedPtBins = 9;
    }

    // fill and fit ratio histogram nPoint times
    int iP = 0;
    while (iP < kNPoints)
    {
      for (int iPtBin = 3; iPtBin < nUsedPtBins + 3; ++iPtBin)
      {
        char hname[100];
        double cutDCAzRnd = 1.;
        if (cutVar)
          cutDCAzRnd = gRandom->Rndm() * 11 * 0.1 + kCutDCAz[0];
        double cutTPCclsRnd = 89;
        if (cutVar)
          cutTPCclsRnd = gRandom->Rndm() * 30 + kCutTPCClusters[0];

        int binCountingFlagRnd = 1;
        if (binCountingVar)
          (gRandom->Rndm() > .4999f) ? binCountingFlagRnd = 1 : binCountingFlagRnd = 0;

        // EXPO AND POL1 SAMPLES ARE SPLIT
        int expFlag = 1;
        if (expVar)
          expFlag = 0;
        //(gRandom->Rndm() > .4999f) ? expFlag = 1 : expFlag = 0;

        int sigmoidFlagRnd = 1;
        if (sigmoidVar)
        {
          if (fRatio.GetBinCenter(iPtBin) < 6.3)
          {
            (gRandom->Rndm() > .4999f) ? sigmoidFlagRnd = 1 : sigmoidFlagRnd = 0;
          }
        }

        if (cutTPCclsRnd < 99.5)
          snprintf(hname, 14, "%1.1f_0%.0f_%d_%d_%d", cutDCAzRnd, cutTPCclsRnd, binCountingFlagRnd, expFlag, sigmoidFlagRnd);
        else
          snprintf(hname, 14, "%1.1f_%.0f_%d_%d_%d", cutDCAzRnd, cutTPCclsRnd, binCountingFlagRnd, expFlag, sigmoidFlagRnd);

        // std::cout<<"hname="<<hname<<std::endl;
        TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", hname, kCentBinsLimitsHe3[iC][0], kCentBinsLimitsHe3[iC][1]));
        fRatio.SetBinContent(iPtBin, h->GetBinContent(iPtBin));
        fRatio.SetBinError(iPtBin, h->GetBinError(iPtBin));
      }
      std::cout << "p=" << iP << std::endl;
      TF1 fitFunc("fitFunc", "pol0");
      auto fit = fRatio.Fit(&fitFunc, "QS");

      int ndf = 10;
      if (iC == 2)
        ndf = 8;
      if (fit->Status() == 0 && fit->Prob() > 0.025 && fit->Prob() < 0.975 && fit->Ndf() == ndf)
      { // check chi2
        fFitPar.Fill(fitFunc.GetParameter(0));
        fProb.Fill(fitFunc.GetProb());
        cdFits->cd();
        fRatio.Write();
        fRatio.Reset();
        ++iP;
      }
    }

    cdHist->cd();
    fFitPar.GetYaxis()->SetTitle("Entries");
    fFitPar.GetXaxis()->SetTitle("R (^{3}#bar{He}/^{3}He)");
    fFitPar.SetDrawOption("histo");
    fFitPar.Rebin(8);
    fFitPar.SetFillStyle(3345);
    fFitPar.SetLineWidth(2);
    fFitPar.SetLineColor(kBlue);
    fFitPar.SetFillColor(kBlue);
    fFitPar.Write();

    TCanvas cFitPar("cFitPar", "cFitPar");
    cFitPar.cd();
    cFitPar.SetTicks(1, 1);
    fFitPar.GetXaxis()->SetRangeUser(0.92,0.99);
    fFitPar.Draw("");
    cFitPar.Print(Form("%s/Systematics_%s.png", kPlotDir, fFitPar.GetName()));

    fProb.GetYaxis()->SetTitle("Entries");
    fProb.GetXaxis()->SetTitle("#it{P}( #chi^{2} > #chi^{2}_{#it{obs}} )");
    fProb.Write();
  }
  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;
}
