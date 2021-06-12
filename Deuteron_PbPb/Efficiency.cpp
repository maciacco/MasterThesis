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

using namespace utils;
using namespace deuteron;

void Efficiency(const char *cutSettings = "", const char *inFileNameMC = "mc", const char *outFileNameEff = "EfficiencyDeuteron")
{
  // make signal extraction plots directory
  system(Form("mkdir %s/efficiency", kPlotDir));

  TFile inFile(Form("%s/%s.root", kDataDir, inFileNameMC));
  TFile outFile(Form("%s/%s.root", kOutDir, outFileNameEff), "RECREATE");

  gStyle->SetOptStat(0);

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    // make plot subdirectory
    system(Form("mkdir %s/efficiency/%s_%s_", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList
    std::string listName = Form("mpuccio_deuterons_%s", cutSettings);
    TTList *list = (TTList *)inFile.Get(listName.data());

    // get histograms from file
    TH2F *fTotal = (TH2F *)list->Get(TString::Format("f%sTotal", kAntimatterMatter[iMatt]).Data());
    TH2F *fITS_TPC_TOF = (TH2F *)list->Get(TString::Format("f%sITS_TPC_TOF", kAntimatterMatter[iMatt]).Data());

    for (int iCent = 0; iCent < kNCentClasses+1; ++iCent)
    { // loop over centrality
      int cent_bin_min = kCentBinsDeuteron[iCent][0];
      int cent_bin_max = kCentBinsDeuteron[iCent][1];
      double cent_bin_lim_min = kCentBinsLimitsDeuteron[iCent][0];
      double cent_bin_lim_max = kCentBinsLimitsDeuteron[iCent][1];

      TH1D *fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
      TH1D *fITS_TPC_TOF_Pt = fITS_TPC_TOF->ProjectionY(TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
      fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(kNPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), kPtBins);
      fITS_TPC_TOF_Pt = (TH1D *)fITS_TPC_TOF_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), kPtBins);
      TH1D fEffPt(TString::Format("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], cent_bin_lim_min, cent_bin_lim_max), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]), kNPtBins, kPtBins);

      for (int iPtBin = 1; iPtBin < fEffPt.GetNbinsX() + 1; ++iPtBin)
      {
        fEffPt.SetBinContent(iPtBin, Eff(fITS_TPC_TOF_Pt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
        fEffPt.SetBinError(iPtBin, EffErr(&fEffPt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
      }
      fEffPt.SetMarkerStyle(20);
      fEffPt.SetMarkerSize(0.8);
      fEffPt.GetYaxis()->SetRangeUser(0., 1.);
      fEffPt.GetXaxis()->SetRangeUser(0., 8.);
      fEffPt.GetXaxis()->SetTitle("#it{p}_{T}");
      fEffPt.GetYaxis()->SetTitle("#epsilon #times Acc");
      fEffPt.SetOption("PE");
      outFile.cd();
      fEffPt.Write();

      // save plot image
      TCanvas canv;
      fEffPt.Draw("");
      canv.Print(Form("%s/efficiency/%s_%s_/cent_%.0f_%.0f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, cent_bin_lim_min, cent_bin_lim_max));
    }
  }
  outFile.Close();
}