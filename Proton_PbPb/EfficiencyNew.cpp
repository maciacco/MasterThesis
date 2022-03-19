// EfficiencyNew.cpp
// This macro computes the efficiency x acceptance correction

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TDirectory.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using namespace utils;
using namespace proton;

void EfficiencyNew(const char *cutSettings = "", const char *inFileNameMC = "mc_20g7_20210929", const char *outFileNameEff = "EfficiencyProton_LongMCTracks", const char *signalOut = "SignalProtonMC", const char *primOut = "PrimaryProtonMC", const char *saving_mode="recreate")
{
  // make signal extraction plots directory
  system(Form("mkdir %s/efficiency", kPlotDir));

  TFile inFilePrimary(Form("%s/%s.root", kOutDir, primOut));
  TFile inFileSignal(Form("%s/%s.root", kOutDir, signalOut));
  TFile inFileMC(Form("%s/%s.root", kDataDir, inFileNameMC));
  TFile outFile(Form("%s/%s.root", kOutDir, outFileNameEff), saving_mode);

  gStyle->SetOptStat(0);

  if (outFile.GetDirectory(Form("%s_", cutSettings)))
    return;
  TDirectory *dirOutFile = outFile.mkdir(Form("%s_", cutSettings));
  dirOutFile->cd();
  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    // make plot subdirectory
    system(Form("mkdir %s/efficiency/%s_%s_", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList
    std::string listName_mcTrue = Form("nuclei_proton_mcTrue_%s", cutSettings);
    TTList *list_mcTrue = (TTList *)inFileMC.Get(listName_mcTrue.data());
    if (!list_mcTrue) return;

    TH2F *fTotal = (TH2F *)list_mcTrue->Get(TString::Format("f%sTotal", kAntimatterMatter[iMatt]).Data());

    TH1D *fTotal_Pt;
    TH1D *fITS_TPC_TOF_Pt;
    TH1D *fITS_TPC_Pt;
    //TF1 *sec_f;
    TH1D *fSec;
    TH1D *fSec_TPC;
    
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    /* for (int iCent = 0; iCent < kNCentClasses; ++iCent) */ // SET FIRST CENTRALITY BIN TO 1 EXCEPT FOR LHC16h7c_g4_2
    {                                                       // loop over centrality
      int cent_bin_min = kCentBinsProton[iCent][0];
      int cent_bin_max = kCentBinsProton[iCent][1];
      double cent_bin_lim_min = kCentBinsLimitsProton[iCent][0];
      double cent_bin_lim_max = kCentBinsLimitsProton[iCent][1];

      fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
      //fITS_TPC_TOF_Pt = fITS_TPC_TOF3->ProjectionY(TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), cent_bin_min, cent_bin_max);
      //sec_f = (TF1 *)inFilePrimary.Get(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fSec = (TH1D*)inFilePrimary.Get(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fSec_TPC = (TH1D*)inFilePrimary.Get(Form("f%sPrimFracTPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fITS_TPC_TOF_Pt = (TH1D*)inFileSignal.Get(Form("%s_%d_%d/f%sTOFrawYield_%.0f_%.0f", cutSettings, 1, 1, kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fITS_TPC_Pt = (TH1D*)inFileSignal.Get(Form("%s_%d_%d/f%sTPCrawYield_%.0f_%.0f", cutSettings, 1, 1, kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      if (!fSec || !fITS_TPC_TOF_Pt || !fTotal_Pt) return;
      //fITS_TPC_TOF_Pt->Multiply(sec_f);
      fITS_TPC_TOF_Pt->Multiply(fSec);
      fITS_TPC_Pt->Multiply(fSec_TPC);
      fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(kNPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), kPtBins);
      fITS_TPC_TOF_Pt = (TH1D *)fITS_TPC_TOF_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_TOF_Pt", kAntimatterMatter[iMatt]), kPtBins);
      fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(kNPtBins, TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), kPtBins);
      TH1D fEffPt(TString::Format("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], cent_bin_lim_min, cent_bin_lim_max), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);
      TH1D fEffPt_TPC(TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], cent_bin_lim_min, cent_bin_lim_max), TString::Format("%s Efficiency #times Acceptance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);

      fEffPt.Divide(fITS_TPC_TOF_Pt, fTotal_Pt, 1, 1, "B");
      fEffPt_TPC.Divide(fITS_TPC_Pt, fTotal_Pt, 1, 1, "B");
      fEffPt.SetMarkerStyle(20);
      fEffPt.SetMarkerSize(0.8);
      fEffPt.GetYaxis()->SetRangeUser(0., 1.);
      fEffPt.GetXaxis()->SetRangeUser(1., 2.);
      fEffPt.GetXaxis()->SetTitle("#it{p}_{T}");
      fEffPt.GetYaxis()->SetTitle("#epsilon #times A");
      fEffPt.SetOption("PE");
      fEffPt_TPC.SetMarkerStyle(20);
      fEffPt_TPC.SetMarkerSize(0.8);
      fEffPt_TPC.GetYaxis()->SetRangeUser(0., 1.);
      fEffPt_TPC.GetXaxis()->SetRangeUser(1., 2.);
      fEffPt_TPC.GetXaxis()->SetTitle("#it{p}_{T}");
      fEffPt_TPC.GetYaxis()->SetTitle("#epsilon #times A");
      fEffPt_TPC.SetOption("PE");
      dirOutFile->cd();
      fEffPt.Write();
      fEffPt_TPC.Write();
      fTotal_Pt->Write();
      fITS_TPC_TOF_Pt->Write();
      fITS_TPC_Pt->Write();
      // save plot image
      TCanvas canv;
      fEffPt.Draw("");
      canv.Print(Form("%s/efficiency/%s_%s_/cent_%.0f_%.0f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, cent_bin_lim_min, cent_bin_lim_max));
    }
  }
  outFile.Close();
}