// EfficiencySec.cpp
// This macro computes the efficiency x acceptance correction (for secondaries)

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

const double kHyperTritonHe3BR = 0.25;

using utils::Eff;
using utils::EffErr;
using namespace he3;

const double hyperTritonToHe3Ratio[2] = {0.42,0.47};

void EfficiencySec(const float cutDCAz = 1.f, const int cutTPCcls = 89, const float cutDCAxy = 0.1f, const float cutChi2TPC = 2.5f, const char *inFileNameMC = "TreeOutMC", const char *outFileNameEff = "EfficiencyHe3SecWd")
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile inFile(Form("%s/%s.root", kResDir, inFileNameMC));
  TFile inFileEff(Form("%s/%s.root", kOutDir, "EfficiencyHe3"));
  TFile outFile(Form("%s/%s.root", kOutDir, outFileNameEff), "RECREATE");

  gStyle->SetOptStat(0);

  // define pt bins
  double pTbins[kNPtBins + 1] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 6.5f, 7.f, 8.f, 10.f};

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    // get histograms from file
    TH2F *fTotal = (TH2F *)inFile.Get(TString::Format("%1.1f_%d_%1.2f_%1.2f/f%sTotalWd", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, kAntimatterMatter[iMatt]));
    TH2F *fITS_TPC = (TH2F *)inFile.Get(TString::Format("%1.1f_%d_%1.2f_%1.2f/f%sITS_TPCwd", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, kAntimatterMatter[iMatt]));

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    { // loop over centrality
      TH1D *fEff = (TH1D*)inFileEff.Get(Form("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      if(!fEff){std::cout<<"No efficiency file!"<<std::endl;return;}
      TH1D *fTotal_Pt = fTotal->ProjectionY(TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      TH1D *fITS_TPC_Pt = fITS_TPC->ProjectionY(TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      TH1D fEffPt(TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), TString::Format("%s #it{f}_{#it{wd}}, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), kNPtBins, pTbins);

      fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(kNPtBins, fTotal_Pt->GetName(), pTbins);
      fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(kNPtBins, fITS_TPC_Pt->GetName(), pTbins);
      if (iCent == 2 || iCent==3)
      {
        int nPtBins = 13;
        double pTbinsNew[] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 7.f, 8.f, 10.f};
        fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(nPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), pTbinsNew);
        fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(nPtBins, TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), pTbinsNew);
        fEffPt = *(TH1D *)fEffPt.Rebin(nPtBins, TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), pTbinsNew);
      }
      if (iCent == 4)
      {
        int nPtBins = 7;
        double pTbinsNew[] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 5.f};
        fTotal_Pt = (TH1D *)fTotal_Pt->Rebin(nPtBins, TString::Format("f%sTotal_Pt", kAntimatterMatter[iMatt]), pTbinsNew);
        fITS_TPC_Pt = (TH1D *)fITS_TPC_Pt->Rebin(nPtBins, TString::Format("f%sITS_TPC_Pt", kAntimatterMatter[iMatt]), pTbinsNew);
        fEffPt = *(TH1D *)fEffPt.Rebin(nPtBins, TString::Format("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), pTbinsNew);
      }

      int nUsedPtBins = 0;
      (iCent <2) ? nUsedPtBins = 13 : nUsedPtBins = 13;
      if (iCent == 4) nUsedPtBins = 7;
      for (int iPtBin = 1; iPtBin < nUsedPtBins + 1; ++iPtBin)
      {
        fEffPt.SetBinContent(iPtBin, kHyperTritonHe3BR * hyperTritonToHe3Ratio[iMatt] * Eff(fITS_TPC_Pt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
        fEffPt.SetBinError(iPtBin, kHyperTritonHe3BR * EffErr(&fEffPt, fTotal_Pt, fEffPt.GetXaxis()->GetBinCenter(iPtBin)));
        fEffPt.SetBinContent(iPtBin, fEffPt.GetBinContent(iPtBin)/fEff->GetBinContent(iPtBin)/(fEffPt.GetBinContent(iPtBin)/fEff->GetBinContent(iPtBin)+1));
        fEffPt.SetBinError(iPtBin, TMath::Sqrt(fEffPt.GetBinError(iPtBin)*fEffPt.GetBinError(iPtBin)/fEffPt.GetBinContent(iPtBin)/fEffPt.GetBinContent(iPtBin)+fEff->GetBinError(iPtBin)*fEff->GetBinError(iPtBin)/fEff->GetBinContent(iPtBin)/fEff->GetBinContent(iPtBin)));
      }
      fEffPt.SetMarkerStyle(20);
      fEffPt.SetMarkerSize(0.8);
      fEffPt.GetXaxis()->SetTitle(kAxisTitlePt);
      fEffPt.GetXaxis()->SetRangeUser(2., 8.);
      fEffPt.GetYaxis()->SetRangeUser(0., 0.09);
      fEffPt.GetYaxis()->SetTitle("#it{f}_{#it{wd}}");
      fEffPt.SetTitle(Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fEffPt.SetOption("PE");
      outFile.cd();
      fEffPt.Write();

      system(Form("mkdir %s/efficiency_wd", kPlotDir));
      TCanvas cEffPt("cEffPt", "cEffPt");
      cEffPt.cd();
      fEffPt.Draw("");
      cEffPt.Print(Form("%s/efficiency_wd/%s.pdf", kPlotDir, fEffPt.GetName()));

      fTotal_Pt->Write();
      fITS_TPC_Pt->Write();
    }
  }
  outFile.Close();
  std::cout << " effWD " << std::endl;
}
