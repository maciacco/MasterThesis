// Spectra.cpp
// This macro computes fully corrected spectra

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using utils::TTList;
using namespace triton;

double he3CorrectionPt(int iMatt, double pt){
  if (iMatt == 1)
    return 1.03776*TMath::Power(pt, -0.00576917);
  return 1.03887*TMath::Power(pt, -0.00675656);
};

void Spectra(const float cutDCAz = 1.f, const int cutTPCcls = 89, const float cutDCAxy = 0.10f, const float cutChi2TPC = 2.5f, const bool binCounting = true, const int bkg_shape = 1, const bool sigmoidCorrection = true, const char *histoNameDir = ".", const char *outFileName = "SpectraHe3", const char *outFileOption = "recreate", const char *dataFile = "AnalysisResults", const char *signalFile = "SignalHe3", const char *effFile = "EfficiencyHe3", const char *primFile = "PrimaryHe3")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTextFont(44);
  gStyle->SetTitleOffset(0);
  gStyle->SetTitleFont(4);
  gStyle->SetTitleFontSize(100);

  TH2F *fNevents;
  TFile *inFileDat = TFile::Open(Form("%s/%s.root", kDataDir, dataFile));
  TTList *fMultList = (TTList *)inFileDat->Get("mpuccio_t_");
  fNevents = (TH2F *)fMultList->Get("fNormalisationHist");

  TFile *inFileRaw = TFile::Open(Form("%s/%s.root", kOutDir, signalFile));
  TFile *inFileEff = TFile::Open(Form("%s/%s.root", kOutDir, effFile));
  TFile *inFileSec = TFile::Open(Form("%s/%s.root", kOutDir, primFile));
  if (!inFileRaw || !inFileEff || !inFileSec)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  outFile.mkdir(histoNameDir);

  TH1D *fRatio[kNCentClasses];
  for (int iCent = 0; iCent < kNCentClasses; ++iCent)
  {
    double pt[]={1.4f,1.6f,2.f,2.4f,3.f};
    fRatio[iCent] = new TH1D(Form("%1.1f_%d_%1.2f_%1.2f_%d_%d/fATOFrawYield_%.0f_%.0f", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, binCounting, bkg_shape, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]),";#it{p}_{T} (GeV/#it{c}); Ratio {}^{3}#bar{H} / ^{3}H",4,pt);
    fRatio[iCent]->Reset();
    fRatio[iCent]->SetName(Form("fRatio_%.0f_%.0f", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    fRatio[iCent]->SetTitle("");
  }

  TH1D *fSpectra[2];
  for (int iCent = 0; iCent < kNCentClasses; ++iCent)
  {
    TH1D *norm;
    norm = fNevents->ProjectionY("norm", kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);

    // compute corrected spectra
    for (int iMatt = 0; iMatt < 2; ++iMatt)
    {
      outFile.cd(histoNameDir);
      TH1D *eff = (TH1D *)inFileEff->Get(Form("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      TF1 *sec_f = nullptr;
      TH1D *sec = nullptr;
      if (iMatt==1)
      {
        sec_f=(TF1 *)inFileSec->Get(Form("f%sSigmoidFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
        sec=(TH1D *)inFileSec->Get(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      }
      TH1D *raw = (TH1D *)inFileRaw->Get(Form("%1.1f_%d_%1.2f_%1.2f_%d_%d/f%sTOFrawYield_%.0f_%.0f", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, binCounting, bkg_shape, kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      if (!eff || !raw){
        std::cout << "no data!" <<std::endl;
        return;
      }

      //sec->Fit(&fitFuncSec,"R");
      fSpectra[iMatt] = new TH1D(*eff);
      int pTbinMax = 4;
      fSpectra[iMatt]->Reset();
      for (int iPtBin = 2; iPtBin < pTbinMax + 1; ++iPtBin)
      {
        double correction = he3CorrectionPt(iMatt, fSpectra[0]->GetBinCenter(iPtBin));  
        double rawYield = raw->GetBinContent(iPtBin);
        std::cout << "raw = "<<rawYield<<std::endl;
        double rawYieldError = raw->GetBinError(iPtBin);
        double efficiency = eff->GetBinContent(iPtBin);
        double effError = eff->GetBinError(iPtBin);
        double primary = 1.;
        if (iMatt==1){
          (!sigmoidCorrection && raw->GetBinCenter(iPtBin) < 2) ? primary = sec->GetBinContent(iPtBin) : primary = sec_f->Eval(raw->GetXaxis()->GetBinCenter(iPtBin));
        }
        fSpectra[iMatt]->SetBinContent(iPtBin, rawYield * primary / efficiency / correction);
        fSpectra[iMatt]->SetBinError(iPtBin, rawYield * primary / efficiency / correction * TMath::Sqrt(rawYieldError * rawYieldError / rawYield / rawYield));
      }
      fSpectra[iMatt]->SetName(Form("f%sSpectra_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSpectra[iMatt]->SetTitle(Form("%s, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSpectra[iMatt]->GetYaxis()->SetTitle("1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}");
      fSpectra[iMatt]->GetXaxis()->SetTitle(kAxisTitlePt);
      fSpectra[iMatt]->GetXaxis()->SetTitleFont(50);

      // scale by number of events
      double events = norm->GetBinContent(4);
      fSpectra[iMatt]->Scale(1. / events, "width");

      // write to file
      fSpectra[iMatt]->Write();
    }

    // compute ratios
    int pTbinMax = 4;
    for (int iPtBin = 2; iPtBin < pTbinMax + 1; ++iPtBin)
    {
      double antiSpec = fSpectra[0]->GetBinContent(iPtBin);
      double spec = fSpectra[1]->GetBinContent(iPtBin);
      double antiSpecErr = fSpectra[0]->GetBinError(iPtBin);
      double specErr = fSpectra[1]->GetBinError(iPtBin);
      if (spec > 1.e-7 && antiSpec > 1.e-7)
      {
        fRatio[iCent]->SetBinContent(iPtBin, antiSpec / spec);
        fRatio[iCent]->SetBinError(iPtBin, antiSpec / spec * TMath::Sqrt(antiSpecErr * antiSpecErr / antiSpec / antiSpec + specErr * specErr / spec / spec));
      }
    }
    fRatio[iCent]->SetLineColor(centrality_colors[iCent]);
    fRatio[iCent]->SetMarkerColor(centrality_colors[iCent]);
    fRatio[iCent]->GetXaxis()->SetTitle(kAxisTitlePt);
    fRatio[iCent]->GetXaxis()->SetTitleOffset(1);
    fRatio[iCent]->GetXaxis()->SetTitleFont(44);
    fRatio[iCent]->GetXaxis()->SetTitleSize(28);
    fRatio[iCent]->GetYaxis()->SetTitle(Form("Ratio %s / %s", kAntimatterMatterLabel[0], kAntimatterMatterLabel[1]));
    fRatio[iCent]->SetTitle(Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    fRatio[iCent]->Fit("pol0");
    fRatio[iCent]->Write();

    TCanvas cRatio("cRatio", "cRatio");
    cRatio.SetTicks(1, 1);
    cRatio.cd();
    fRatio[iCent]->GetXaxis()->SetRangeUser(1.6, 3.);
    fRatio[iCent]->GetYaxis()->SetRangeUser(0., 1.8);
    fRatio[iCent]->Draw("");
    TLatex chi2(2.2, 1.45, Form("#chi^{2}/NDF = %.2f/%d", fRatio[iCent]->GetFunction("pol0")->GetChisquare(), fRatio[iCent]->GetFunction("pol0")->GetNDF()));
    chi2.SetTextSize(28);
    TLatex p0(2.2, 1.6, Form("R = %.3f #pm %.3f", fRatio[iCent]->GetFunction("pol0")->GetParameter(0), fRatio[iCent]->GetFunction("pol0")->GetParError(0)));
    p0.SetTextSize(28);
    chi2.Draw("same");
    p0.Draw("same");
    cRatio.Print(Form("%s/%s.pdf", kPlotDir, fRatio[iCent]->GetName()));
  }
  outFile.Close();
}
