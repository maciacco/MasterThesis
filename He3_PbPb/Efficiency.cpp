// Efficiency.cpp
// This macro computes the efficiency x acceptance correction

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "Utils.h"
#include "Config.h"

using utils::Eff;
using utils::EffErr;

void Efficiency(const float cutDCAz = 1.f, const int cutTPCcls = 89, const char *inFileNameMC = "TreeOutMC", const char *outFileNameEff = "EfficiencyHe3")
{
  // make signal extraction plots directory
  system(Form("mkdir %s/efficiency", kPlotDir));

  TFile inFile(Form("%s/%s.root", kResDir, inFileNameMC));
  TFile outFile(Form("%s/%s.root", kOutDir, outFileNameEff), "RECREATE");

  gStyle->SetOptStat(0);

  // define pt bins
  double pTbins[kNPtBins + 1] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 6.5f, 7.f, 8.f, 10.f};

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    // make plot subdirectory
    system(Form("mkdir %s/efficiency/%s_%1.1f_%d", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls));

    // get histograms from file
    TH2F *fTotal = (TH2F *)inFile.Get(TString::Format("%.1f_%d_/f%sTotal", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));
    TH2F *fITS_TPC = (TH2F *)inFile.Get(TString::Format("%.1f_%d_/f%sITS_TPC", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    { // loop over centrality
      TH1D *fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      TH1D *fITS_TPC_Pt = fITS_TPC->ProjectionY(TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(kNPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), pTbins);
      fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), pTbins);
      TH1D fEffPt(TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), kNPtBins, pTbins);

      if (iCent == kNCentClasses - 1)
      {
        int nPtBins = 13;
        double pTbinsNew[] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 7.f, 8.f, 10.f};
        fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(nPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), pTbinsNew);
        fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(nPtBins, TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), pTbinsNew);
        fEffPt = *(TH1D *)fEffPt.Rebin(nPtBins, TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), pTbinsNew);
      }

      for (int iPtBin = 1; iPtBin < fEffPt.GetNbinsX() + 1; ++iPtBin)
      {
        fEffPt.SetBinContent(iPtBin, Eff(fITS_TPC_Pt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
        fEffPt.SetBinError(iPtBin, EffErr(&fEffPt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
      }
      fEffPt.SetMarkerStyle(20);
      fEffPt.SetMarkerSize(0.8);
      fEffPt.GetXaxis()->SetTitle(kAxisTitlePt);
      fEffPt.GetYaxis()->SetTitle("#epsilon #times Acc");
      fEffPt.SetOption("PE");
      outFile.cd();
      fEffPt.Write();

      // save plot image
      TCanvas canv;
      fEffPt.Draw("");
      canv.Print(Form("%s/efficiency/%s_%1.1f_%d/cent_%.0f_%.0f.png", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    }
  }
  outFile.Close();
}