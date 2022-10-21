// Efficiency.cpp
// This macro computes the efficiency x acceptance correction

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using utils::Eff;
using utils::EffErr;
using namespace triton;

void Efficiency(const float cutDCAz = 1.f, const int cutTPCcls = 89, const float cutDCAxy = 0.1f, const float cutChi2TPC = 2.5f, const char *inFileNameMC = "TreeOutMC", const char *outFileNameEff = "EfficiencyHe3")
{
  // make signal extraction plots directory
  system(Form("mkdir %s/efficiency", kPlotDir));

  TFile inFile(Form("%s/%s.root", kResDir, inFileNameMC));
  TFile outFile(Form("%s/%s.root", kOutDir, outFileNameEff), "RECREATE");

  gStyle->SetOptStat(0);

  // define pt bins
  double pTbins[kNPtBins + 1] = {1.4f,1.6f,2.f,2.4f,3.f};

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    // make plot subdirectory
    system(Form("mkdir %s/efficiency/%s_%1.1f_%d_%1.2f_%1.2f", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC));

    // get histograms from file
    TH2F *fTotal = (TH2F *)inFile.Get(TString::Format("%.1f_%d_%.2f_%1.2f/f%sTotal", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, kAntimatterMatter[iMatt]));
    TH2F *fITS_TPC = (TH2F *)inFile.Get(TString::Format("%.1f_%d_%.2f_%1.2f/f%sITS_TPC", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, kAntimatterMatter[iMatt]));

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    { // loop over centrality
      TH1D *fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      TH1D *fITS_TPC_Pt = fITS_TPC->ProjectionY(TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(kNPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), pTbins);
      fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), pTbins);
      TH1D fEffPt(TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), kNPtBins, pTbins);

      for (int iPtBin = 2; iPtBin < fEffPt.GetNbinsX() + 1; ++iPtBin)
      {
        fEffPt.SetBinContent(iPtBin, Eff(fITS_TPC_Pt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
        fEffPt.SetBinError(iPtBin, EffErr(&fEffPt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
      }
      fEffPt.SetMarkerStyle(20);
      fEffPt.SetMarkerSize(0.8);
      fEffPt.GetXaxis()->SetTitle(kAxisTitlePt);
      fEffPt.GetYaxis()->SetTitle("#epsilon #times A");
      fEffPt.SetOption("PE");
      outFile.cd();
      fEffPt.Write();

      // save plot image
      TCanvas canv;
      fEffPt.Draw("");
      canv.Print(Form("%s/efficiency/%s_%1.1f_%d_%1.1f/cent_%.0f_%.0f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls, cutDCAxy, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    }
  }
  outFile.Close();
}