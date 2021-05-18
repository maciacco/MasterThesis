// SignalLoss.cpp
// This macro computes ratio between identified He3 and candidates in data and MC (primaries)

#include <string>

#include <Riostream.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TString.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "../utils/Config.h"
#include "../utils/Utils.h"

using namespace he3;

void SignalLoss(const float cutDCAz = 1.f, const int cutTPCcls = 89, const bool binCounting = true, const int bkg_shape = 1, const bool binomial = true, const char *inFileDatName = "SignalHe3", const char *inFileMCName = "TreeOutMC", const char *outFileName = "SignalLoss")
{
  TFile *inFileDatHe3PID = TFile::Open(Form("out/%s_He3PID.root", inFileDatName));
  TFile *inFileDatHe4PID = TFile::Open(Form("out/%s_AlphaPID.root", inFileDatName));
  TFile *inFileDatNoPID = TFile::Open(Form("out/%s_NoPID.root", inFileDatName));
  TFile *inFileMC = TFile::Open(Form("results/%s_He3PID.root", inFileMCName));
  TFile *outFile = TFile::Open(Form("out/%s.root", outFileName), "recreate");

  for (int iMatt = 0; iMatt < 1; ++iMatt)
  {
    // get mc histograms
    TH2F *fPIDMc = (TH2F *)inFileMC->Get(Form("%.0f_%d_/f%sPID", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));
    TH2F *fAlphaPIDMc = (TH2F *)inFileMC->Get(Form("%.0f_%d_/f%sAlphaPID", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));
    TH2F *fNoPIDMc = (TH2F *)inFileMC->Get(Form("%.0f_%d_/f%sNoPID", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      // get data histograms
      TH1D *fPIDdatPt = (TH1D *)inFileDatHe3PID->Get(Form("%.0f_%d_%d_%d/f%sTPCrawYield_%.0f_%.0f", cutDCAz, cutTPCcls, binCounting, bkg_shape, kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      TH1D *fAlphaPIDdatPt = (TH1D *)inFileDatHe4PID->Get(Form("%.0f_%d_%d_%d/f%sTPCrawYield_%.0f_%.0f", cutDCAz, cutTPCcls, binCounting, bkg_shape, kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      TH1D *fNoPIDdatPt = (TH1D *)inFileDatNoPID->Get(Form("%.0f_%d_%d_%d/f%sTPCrawYield_%.0f_%.0f", cutDCAz, cutTPCcls, binCounting, bkg_shape, kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      // projections
      TH1D *pidMcPt = fPIDMc->ProjectionY(Form("f%sPIDmc_Pt_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      TH1D *alphaPIDMcPt = fAlphaPIDMc->ProjectionY(Form("f%sAlphaPIDmc_Pt_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);
      TH1D *noPIDMcPt = fNoPIDMc->ProjectionY(Form("f%sNoPIDmc_Pt_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);

      int nUsedPtBins = 12;
      if (iCent == kNCentClasses - 1)
      {
        nUsedPtBins = 9;
        int nPtBins = 13;
        double pTbinsNew[] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 7.f, 8.f, 10.f};
        pidMcPt = (TH1D *)pidMcPt->Rebin(nPtBins, pidMcPt->GetName(), pTbinsNew);
        noPIDMcPt = (TH1D *)noPIDMcPt->Rebin(nPtBins, noPIDMcPt->GetName(), pTbinsNew);
      }

      // ratio
      TH1D fRatioDat(*fPIDdatPt);
      TH1D fRatioAlphaDat(*fPIDdatPt);
      TH1D fRatioMc(*fPIDdatPt);
      TH1D fRatioAlphaMc(*fPIDdatPt);
      fRatioDat.Reset();
      fRatioDat.SetTitle(Form("Centrality: %.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioDat.GetYaxis()->SetTitle("#it{N_{pid}} / #it{N_{tot}}");
      fRatioDat.SetName(Form("f%sRatioDat_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioDat.GetYaxis()->SetRangeUser(0., 1.15);
      fRatioAlphaDat.Reset();
      fRatioAlphaDat.SetTitle(Form("Centrality: %.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioAlphaDat.GetYaxis()->SetTitle("#it{N_{pid}} / #it{N_{tot}}");
      fRatioAlphaDat.SetName(Form("f%sRatioDat_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioAlphaDat.GetYaxis()->SetRangeUser(0., 1.15);
      fRatioMc.Reset();
      fRatioMc.SetTitle(Form("Centrality: %.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioMc.GetYaxis()->SetTitle("#it{N_{pid}} / #it{N_{tot}}");
      fRatioMc.SetName(Form("f%sRatioMc_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioMc.GetYaxis()->SetRangeUser(0., 1.15);
      fRatioAlphaMc.Reset();
      fRatioAlphaMc.SetTitle(Form("Centrality: %.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioAlphaMc.GetYaxis()->SetTitle("#it{N_{pid}} / #it{N_{tot}}");
      fRatioAlphaMc.SetName(Form("f%sRatioAlphaMc_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fRatioAlphaMc.GetYaxis()->SetRangeUser(0., 1.15);

      for (int iPtBin = 3; iPtBin < nUsedPtBins + 5; ++iPtBin)
      {
        if (fNoPIDdatPt->GetBinContent(iPtBin) > 0)
        {
          double dat_pid = (fPIDdatPt->GetBinContent(iPtBin) > 1.e-1) ? fPIDdatPt->GetBinContent(iPtBin) : 0.;
          double dat_alpha_pid = (fAlphaPIDdatPt->GetBinContent(iPtBin) > 1.e-1) ? fAlphaPIDdatPt->GetBinContent(iPtBin) : 0.;
          double dat_nopid = fNoPIDdatPt->GetBinContent(iPtBin);
          fRatioDat.SetBinContent(iPtBin, dat_pid / dat_nopid);
          fRatioAlphaDat.SetBinContent(iPtBin, dat_alpha_pid / dat_nopid);

          // errors
          (dat_pid > 1.e-1) ? fRatioDat.SetBinError(iPtBin, TMath::Sqrt(1 / dat_pid + 1 / dat_nopid)) : fRatioDat.SetBinError(iPtBin, 0);
          (dat_alpha_pid > 1.e-1) ? fRatioAlphaDat.SetBinError(iPtBin, TMath::Sqrt(1 / dat_alpha_pid + 1 / dat_nopid)) : fRatioAlphaDat.SetBinError(iPtBin, 0);
          if (binomial)
          {
            fRatioDat.SetBinError(iPtBin, TMath::Sqrt(TMath::Abs((dat_pid / dat_nopid) * (1 - dat_pid / dat_nopid)) / dat_nopid));
            fRatioAlphaDat.SetBinError(iPtBin, TMath::Sqrt(TMath::Abs((dat_alpha_pid / dat_nopid) * (1 - dat_alpha_pid / dat_nopid)) / dat_nopid));
          }
        }
        if (noPIDMcPt->GetBinContent(iPtBin) > 0)
        {
          double mc_pid = (pidMcPt->GetBinContent(iPtBin) > 1.e-1) ? pidMcPt->GetBinContent(iPtBin) : 0.;
          double mc_alpha_pid = (alphaPIDMcPt->GetBinContent(iPtBin) > 1.e-1) ? alphaPIDMcPt->GetBinContent(iPtBin) : 0.;
          double mc_nopid = noPIDMcPt->GetBinContent(iPtBin);
          fRatioMc.SetBinContent(iPtBin, mc_pid / mc_nopid);
          fRatioAlphaMc.SetBinContent(iPtBin, mc_alpha_pid / mc_nopid);

          // errors
          (mc_pid > 1.e-1) ? fRatioMc.SetBinError(iPtBin, TMath::Sqrt(1 / mc_pid + 1 / mc_nopid)) : fRatioMc.SetBinError(iPtBin, 0);
          (mc_alpha_pid > 1.e-1) ? fRatioAlphaMc.SetBinError(iPtBin, TMath::Sqrt(1 / mc_alpha_pid + 1 / mc_nopid)) : fRatioAlphaMc.SetBinError(iPtBin, 0);
          if (binomial)
          {
            fRatioMc.SetBinError(iPtBin, TMath::Sqrt(TMath::Abs((mc_pid / mc_nopid) * (1 - mc_pid / mc_nopid)) / mc_nopid));
            fRatioAlphaMc.SetBinError(iPtBin, TMath::Sqrt(TMath::Abs((mc_alpha_pid / mc_nopid) * (1 - mc_alpha_pid / mc_nopid)) / mc_nopid));
          }
        }
      }
      outFile->cd();
      fRatioDat.Write();
      fRatioAlphaDat.Write();
      fRatioMc.Write();
      fRatioAlphaMc.Write();

      // plot superposition
      TCanvas canv(Form("f%sRatio_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), Form("%sratio_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      TLegend leg(0.197708, 0.43949, 0.462751, 0.562633);
      leg.AddEntry(&fRatioMc, "MC, trackingPID == kHe3");
      leg.AddEntry(&fRatioDat, "Data, trackingPID == kHe3");
      leg.AddEntry(&fRatioAlphaMc, "MC, trackingPID == kHe4");
      leg.AddEntry(&fRatioAlphaDat, "Data, trackingPID == kHe4");

      fRatioMc.SetMarkerSize(0.8);
      fRatioMc.SetMarkerStyle(20);
      fRatioMc.SetMarkerColor(kRed);
      fRatioMc.SetLineColor(kRed);
      fRatioMc.Draw("pe");
      fRatioDat.SetMarkerSize(0.8);
      fRatioDat.SetMarkerStyle(20);
      fRatioDat.SetMarkerColor(kGreen + 2);
      fRatioDat.SetLineColor(kGreen + 2);
      fRatioDat.Draw("pesame");
      fRatioAlphaMc.SetMarkerSize(0.8);
      fRatioAlphaMc.SetMarkerStyle(20);
      fRatioAlphaMc.SetMarkerColor(kOrange);
      fRatioAlphaMc.SetLineColor(kOrange);
      fRatioAlphaMc.Draw("pesame");
      fRatioAlphaDat.SetMarkerSize(0.8);
      fRatioAlphaDat.SetMarkerStyle(20);
      fRatioAlphaDat.SetMarkerColor(kBlue + 2);
      fRatioAlphaDat.SetLineColor(kBlue + 2);
      fRatioAlphaDat.Draw("pesame");
      leg.Draw("same");
      canv.Write();
    }
  }
  outFile->Close();
}