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

void Efficiency(const char *cutSettings = "", const char *inFileNameMC = "mc", const char *outFileNameEff = "EfficiencyProton")
{
  // make signal extraction plots directory
  system(Form("mkdir %s/efficiency", kPlotDir));

  /* TFile inFile1(Form("%s/%s_finePtBinning1.root", kDataDir, inFileNameMC));
  TFile inFile2(Form("%s/%s_finePtBinning2.root", kDataDir, inFileNameMC));
  TFile inFile3(Form("%s/%s_finePtBinning3.root", kDataDir, inFileNameMC)); */
  TFile inFile1(Form("%s/%s.root", kDataDir, inFileNameMC));
  TFile outFile(Form("%s/%s.root", kOutDir, outFileNameEff), "RECREATE");

  gStyle->SetOptStat(0);

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    // make plot subdirectory
    system(Form("mkdir %s/efficiency/%s_%s_", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList
    std::string listName = Form("nuclei_proton_%s", cutSettings);
    TTList *list1 = (TTList *)inFile1.Get(listName.data());
   /*  TTList *list2 = (TTList *)inFile2.Get(listName.data());
    TTList *list3 = (TTList *)inFile3.Get(listName.data()); */

    // get histograms from file
    TH2F *fTotal1 = (TH2F *)list1->Get(TString::Format("f%sTotal", kAntimatterMatter[iMatt]).Data());
    TH2F *fITS_TPC_TOF1 = (TH2F *)list1->Get(TString::Format("f%sITS_TPC_TOF", kAntimatterMatter[iMatt]).Data());
    /* TH2F *fTotal2 = (TH2F *)list2->Get(TString::Format("f%sTotal", kAntimatterMatter[iMatt]).Data());
    TH2F *fITS_TPC_TOF2 = (TH2F *)list2->Get(TString::Format("f%sITS_TPC_TOF", kAntimatterMatter[iMatt]).Data());
    TH2F *fTotal3 = (TH2F *)list3->Get(TString::Format("f%sTotal", kAntimatterMatter[iMatt]).Data());
    TH2F *fITS_TPC_TOF3 = (TH2F *)list3->Get(TString::Format("f%sITS_TPC_TOF", kAntimatterMatter[iMatt]).Data()); */

    TH2F *fTotal = (TH2F *)fTotal1->Clone(fTotal1->GetName());
    //fTotal->Add(fTotal2);
    //fTotal->Add(fTotal3);
    TH2F *fITS_TPC_TOF = (TH2F *)fITS_TPC_TOF1->Clone(fITS_TPC_TOF1->GetName());
    //fITS_TPC_TOF->Add(fITS_TPC_TOF2);
    //fITS_TPC_TOF->Add(fITS_TPC_TOF3);

    for (int iCent = 0; iCent < kNCentClasses+1; ++iCent) // SET FIRST CENTRALITY BIN TO 1 EXCEPT FOR LHC16h7c_g4_2
    { // loop over centrality
      int cent_bin_min = kCentBinsProton[iCent][0];
      int cent_bin_max = kCentBinsProton[iCent][1];
      double cent_bin_lim_min = kCentBinsLimitsProton[iCent][0];
      double cent_bin_lim_max = kCentBinsLimitsProton[iCent][1];

      TH1D *fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
      TH1D *fITS_TPC_TOF_Pt = fITS_TPC_TOF->ProjectionY(TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
      fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(kNPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), kPtBins);
      fITS_TPC_TOF_Pt = (TH1D *)fITS_TPC_TOF_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), kPtBins);
      TH1D fEffPt(TString::Format("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], cent_bin_lim_min, cent_bin_lim_max), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);

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