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
using namespace he3;

double he3CorrectionPt(int iMatt, double pt){
  double fit_c_proton = 0.868419;
  double fit_c_proton_error = 0.0579547;
  /* if (iMatt == 1)
    return 0.99274*TMath::Power(pt,0.00143);
  return 1.04948*TMath::Power(pt,-0.01525); */
  //return 1;
  double f=1.-(1./0.029)* /* (0.738506/1.058)*(0.738506-1.) */ (1.058-1.)*0.00294*TMath::Power(pt,-0.19483);
  double f_note=0.99274*TMath::Power(pt,0.00143);
  double f_new = 0;
  if (iMatt == 1){
    f_new = 1.-(1./0.029)*(fit_c_proton/1.058)*(fit_c_proton-1.)*0.00294*TMath::Power(pt,-0.19483);
  }
  else f_new = 1.+(0.02088*TMath::Power(pt,-0.48766))/(0.084)*(1.-0.83);
  //return f_new*(1+(f-f_note)/f_note);
  return 1;
};

void Spectra(const float cutDCAz = 1.f, const int cutTPCcls = 89, const bool binCounting = true, const int bkg_shape = 1, const bool sigmoidCorrection = true, const char *histoNameDir = ".", const char *outFileName = "SpectraHe3", const char *outFileOption = "recreate", const char *dataFile = "AnalysisResults", const char *signalFile = "SignalHe3", const char *effFile = "EfficiencyHe3", const char *primFile = "PrimaryHe3")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTextFont(44);

  TH2F *fNevents[kNDataFiles];
  for (int iD = 0; iD < kNDataFiles; ++iD)
  {
    TFile *inFileDat = TFile::Open(Form("%s/%s_%s.root", kDataDir, dataFile, kDataFileLabel[iD]));
    TTList *fMultList = (TTList *)inFileDat->Get("mpuccio_he3_");
    fNevents[iD] = (TH2F *)fMultList->Get("fNormalisationHist");
  }
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
    fRatio[iCent] = new TH1D(*(TH1D *)inFileRaw->Get(Form("%1.1f_%d_%d_%d/fATPCrawYield_%.0f_%.0f", cutDCAz, cutTPCcls, binCounting, bkg_shape, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1])));
    fRatio[iCent]->Reset();
    fRatio[iCent]->SetName(Form("fRatio_%.0f_%.0f", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    fRatio[iCent]->SetTitle("");
  }

  TH1D *fSpectra[2];
  for (int iCent = 0; iCent < kNCentClasses; ++iCent)
  {
    TH1D *norm[2];
    for (int iD = 0; iD < kNDataFiles; ++iD)
      norm[iD] = fNevents[iD]->ProjectionY("norm", kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1]);

    // compute corrected spectra
    for (int iMatt = 0; iMatt < 2; ++iMatt)
    {
      outFile.cd(histoNameDir);
      TH1D *eff = (TH1D *)inFileEff->Get(Form("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      TF1 *sec_f = (TF1 *)inFileSec->Get(Form("f%sSigmoidFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      TH1D *sec = (TH1D *)inFileSec->Get(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      TH1D *raw = (TH1D *)inFileRaw->Get(Form("%1.1f_%d_%d_%d/f%sTPCrawYield_%.0f_%.0f", cutDCAz, cutTPCcls, binCounting, bkg_shape, kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));

      //sec->Fit(&fitFuncSec,"R");
      fSpectra[iMatt] = new TH1D(*eff);
      int pTbinMax = 11;
      if (iCent < kNCentClasses - 1)
        pTbinMax = 13;
      fSpectra[iMatt]->Reset();
      for (int iPtBin = 3; iPtBin < pTbinMax + 1; ++iPtBin)
      {
        double correction = he3CorrectionPt(iMatt, fSpectra[0]->GetBinCenter(iPtBin));  
        double rawYield = raw->GetBinContent(iPtBin);
        double rawYieldError = raw->GetBinError(iPtBin);
        double efficiency = eff->GetBinContent(iPtBin);
        double effError = eff->GetBinError(iPtBin);

        double primary = 0.;
        (!sigmoidCorrection && raw->GetBinCenter(iPtBin) < 6.5) ? primary = sec->GetBinContent(iPtBin) : primary = sec_f->Eval(raw->GetXaxis()->GetBinCenter(iPtBin));
        fSpectra[iMatt]->SetBinContent(iPtBin, rawYield * primary / efficiency / correction);
        fSpectra[iMatt]->SetBinError(iPtBin, rawYield * primary / efficiency / correction * TMath::Sqrt(effError * effError / efficiency / efficiency + rawYieldError * rawYieldError / rawYield / rawYield));
      }
      fSpectra[iMatt]->SetName(Form("f%sSpectra_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSpectra[iMatt]->SetTitle(Form("%s, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSpectra[iMatt]->GetYaxis()->SetTitle("1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}");
      fSpectra[iMatt]->GetXaxis()->SetTitle(kAxisTitlePt);

      // scale by number of events
      double events = 0.f;
      for (int iD = 0; iD < kNDataFiles; ++iD)
        events += norm[iD]->GetBinContent(4);
      fSpectra[iMatt]->Scale(1. / events, "width");

      // write to file
      fSpectra[iMatt]->Write();
    }

    // compute ratios
    int pTbinMax = 11;
    if (iCent < kNCentClasses - 1)
      pTbinMax = 13;
    for (int iPtBin = 3; iPtBin < pTbinMax + 1; ++iPtBin)
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
    fRatio[iCent]->GetXaxis()->SetTitle(kAxisTitlePt);
    fRatio[iCent]->GetYaxis()->SetTitle(Form("Ratio %s / %s", kAntimatterMatterLabel[0], kAntimatterMatterLabel[1]));
    fRatio[iCent]->SetTitle(Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    fRatio[iCent]->Fit("pol0");
    fRatio[iCent]->Write();

    TCanvas cRatio("cRatio", "cRatio");
    cRatio.SetTicks(1, 1);
    cRatio.cd();
    fRatio[iCent]->GetXaxis()->SetRangeUser(2., 8.);
    fRatio[iCent]->GetYaxis()->SetRangeUser(0., 1.8);
    fRatio[iCent]->Draw("");
    TLatex chi2(5.8, 1.45, Form("#chi^{2}/NDF = %.2f/%d", fRatio[iCent]->GetFunction("pol0")->GetChisquare(), fRatio[iCent]->GetFunction("pol0")->GetNDF()));
    chi2.SetTextSize(28);
    TLatex p0(5.8, 1.6, Form("R = %.3f #pm %.3f", fRatio[iCent]->GetFunction("pol0")->GetParameter(0), fRatio[iCent]->GetFunction("pol0")->GetParError(0)));
    p0.SetTextSize(28);
    chi2.Draw("same");
    p0.Draw("same");
    cRatio.Print(Form("%s/%s.pdf", kPlotDir, fRatio[iCent]->GetName()));
  }
  outFile.Close();
}
