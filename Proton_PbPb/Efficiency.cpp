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
using namespace proton;

void Efficiency(const char *cutSettings = "", const char *inFileNameMC = "AnalysisResults_LHC21l5_full_largeDCA_cutChi2", const char *outFileNameEff = "EfficiencyProton_OldMethod")
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
    std::string listName = Form("nuclei_proton_mcTrue_%s", cutSettings);
    TTList *list = (TTList *)inFile.Get(listName.data());

    // get histograms from file
    TH2F *fTotal = (TH2F *)list->Get(TString::Format("f%sTotal", kAntimatterMatter[iMatt]).Data());
    TH2F *fITS_TPC = (TH2F *)list->Get(TString::Format("f%sITS_TPC", kAntimatterMatter[iMatt]).Data());
    TH2F *fITS_TPC_TOF = (TH2F *)list->Get(TString::Format("f%sITS_TPC_TOF", kAntimatterMatter[iMatt]).Data());

    for (int iCent = 0; iCent < kNCentClasses + 1; ++iCent) // SET FIRST CENTRALITY BIN TO 1 EXCEPT FOR LHC16h7c_g4_2
    {                                                       // loop over centrality
      int cent_bin_min = kCentBinsProton[iCent][0];
      int cent_bin_max = kCentBinsProton[iCent][1];
      double cent_bin_lim_min = kCentBinsLimitsProton[iCent][0];
      double cent_bin_lim_max = kCentBinsLimitsProton[iCent][1];

      TH1D *fTotal_Pt;
      TH1D *fITS_TPC_TOF_Pt;
      TH1D *fITS_TPC_Pt;
      
      if (iCent < -1)
      {
        /* fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
        fITS_TPC_TOF_Pt = fITS_TPC_TOF->ProjectionY(TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max); */
      }
      else
      {
        fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
        fITS_TPC_TOF_Pt = fITS_TPC_TOF->ProjectionY(TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
        fITS_TPC_Pt = fITS_TPC->ProjectionY(TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
      }
      fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(kNPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), kPtBins);
      fITS_TPC_TOF_Pt = (TH1D *)fITS_TPC_TOF_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), kPtBins);
      fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), kPtBins);
      TH1D fEffPt(TString::Format("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], cent_bin_lim_min, cent_bin_lim_max), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);
      TH1D fEffPt_TPC(TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], cent_bin_lim_min, cent_bin_lim_max), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);

      fEffPt.Divide(fITS_TPC_TOF_Pt, fTotal_Pt, 1, 1, "B");
      fEffPt.SetMarkerStyle(20);
      fEffPt.SetMarkerSize(0.8);
      fEffPt.GetYaxis()->SetRangeUser(0., 1.);
      fEffPt.GetXaxis()->SetRangeUser(1., 2.);
      fEffPt.GetXaxis()->SetTitle("#it{p}_{T}");
      fEffPt.GetYaxis()->SetTitle("#epsilon #times A");
      fEffPt.SetOption("PE");
      outFile.cd();
      fEffPt.Write();

      fEffPt_TPC.Divide(fITS_TPC_Pt, fTotal_Pt, 1, 1, "B");
      fEffPt_TPC.SetMarkerStyle(20);
      fEffPt_TPC.SetMarkerSize(0.8);
      fEffPt_TPC.GetYaxis()->SetRangeUser(0., 1.);
      fEffPt_TPC.GetXaxis()->SetRangeUser(1., 2.);
      fEffPt_TPC.GetXaxis()->SetTitle("#it{p}_{T}");
      fEffPt_TPC.GetYaxis()->SetTitle("#epsilon #times A");
      fEffPt_TPC.SetOption("PE");
      outFile.cd();
      fEffPt_TPC.Write();

      // save plot image
      TCanvas canv;
      fEffPt.Draw("");
      canv.Print(Form("%s/efficiency/%s_%s_/cent_%.0f_%.0f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, cent_bin_lim_min, cent_bin_lim_max));
    }
  }
  outFile.Close();
}