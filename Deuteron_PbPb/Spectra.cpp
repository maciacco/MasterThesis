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
using namespace deuteron;

void Spectra(const char *cutSettings = "", const bool binCounting = false, const int bkg_shape = 1, const bool sigmoidCorrection = true, const char *histoNameDir = ".", const char *outFileName = "SpectraDeuteron", const char *outFileOption = "recreate", const char *dataFile = "AnalysisResults", const char *signalFile = "SignalDeuteron", const char *effFile = "EfficiencyDeuteron", const char *primFile = "PrimaryDeuteron")
{
  gStyle->SetOptFit(0);

  TH2F *fNevents;
  TFile *inFileDat = TFile::Open(Form("%s/%s.root", kDataDir, dataFile));
  TTList *fMultList = (TTList *)inFileDat->Get("mpuccio_deuterons_");
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
    fRatio[iCent] = new TH1D(*(TH1D *)inFileRaw->Get(Form("%s_%d_%d/fATOFrawYield_%.0f_%.0f", cutSettings, binCounting, bkg_shape, kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1])));
    fRatio[iCent]->Reset();
    fRatio[iCent]->SetName(Form("fRatio_%.0f_%.0f", kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
    fRatio[iCent]->SetTitle("");
  }

  TH1D *fSpectra[2];
  for (int iCent = 0; iCent < kNCentClasses; ++iCent)
  {
    TH1D *norm;
    norm = fNevents->ProjectionY("norm", kCentBinsDeuteron[iCent][0], kCentBinsDeuteron[iCent][1]);

    // compute corrected spectra
    for (int iMatt = 0; iMatt < 2; ++iMatt)
    {
      outFile.cd(histoNameDir);
      TH1D *eff = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      TF1 *sec_f = (TF1 *)inFileSec->Get(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      TH1D *sec = (TH1D *)inFileSec->Get(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      TH1D *raw = (TH1D *)inFileRaw->Get(Form("%s_%d_%d/f%sTOFrawYield_%.0f_%.0f", cutSettings, binCounting, bkg_shape, kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));

      //sec->Fit(&fitFuncSec,"R");
      fSpectra[iMatt] = new TH1D(*eff);
      int pTbinMax = 22;
      std::cout<<"entering pt loop..."<<std::endl;
      for (int iPtBin = 4; iPtBin < pTbinMax + 1; ++iPtBin)
      {
        double rawYield = raw->GetBinContent(iPtBin);
        double rawYieldError = raw->GetBinError(iPtBin);
        double efficiency = eff->GetBinContent(iPtBin);
        double effError = eff->GetBinError(iPtBin);

        double primary = 1.;
        double primaryError = 0.;
        if(iMatt == 1){
          primary = sec->GetBinContent(iPtBin);
          primaryError = sec->GetBinError(iPtBin);
          //primary = sec_f->Eval(raw->GetXaxis()->GetBinCenter(iPtBin));
        }
        fSpectra[iMatt]->SetBinContent(iPtBin, rawYield * primary / efficiency );
        fSpectra[iMatt]->SetBinError(iPtBin, (rawYield * primary / efficiency) * TMath::Sqrt(primaryError * primaryError / primary / primary + effError * effError / efficiency / efficiency + rawYieldError * rawYieldError / rawYield / rawYield));

        std::cout<<"eff="<<efficiency<<"; raw="<<rawYield<<"; rawError="<<rawYieldError<<"; primary="<<primary<<std::endl;
      }
      fSpectra[iMatt]->SetName(Form("f%sSpectra_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fSpectra[iMatt]->SetTitle(Form("%s, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fSpectra[iMatt]->GetYaxis()->SetTitle("1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}");
      fSpectra[iMatt]->GetXaxis()->SetTitle(kAxisTitlePt);

      // scale by number of events
      double events = norm->GetBinContent(4);
      fSpectra[iMatt]->Scale(1. / events, "width");

      // write to file
      fSpectra[iMatt]->Write();
    }

    // compute ratios
    int pTbinMax = 22;
    for (int iPtBin = 1; iPtBin < pTbinMax + 1; ++iPtBin)
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
    fRatio[iCent]->GetYaxis()->SetTitle(Form("Ratio %s/%s", kAntimatterMatterLabel[0], kAntimatterMatterLabel[1]));
    fRatio[iCent]->SetTitle(Form("%.0f-%.0f%%", kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
    fRatio[iCent]->Fit("pol0");
    fRatio[iCent]->Write();
    
    TCanvas cRatio("cRatio", "cRatio");
    cRatio.SetTicks(1, 1);
    cRatio.cd();
    fRatio[iCent]->GetXaxis()->SetRangeUser(1.,5.);
    fRatio[iCent]->GetYaxis()->SetRangeUser(0., 1.4);
    fRatio[iCent]->Draw("");
    TLatex chi2(3.5, 1.15, Form("#chi^{2}/NDF = %.2f/%d", fRatio[iCent]->GetFunction("pol0")->GetChisquare(), fRatio[iCent]->GetFunction("pol0")->GetNDF()));
    chi2.SetTextSize(0.035);
    TLatex p0(3.5, 1.3, Form("R = %.3f #pm %.3f", fRatio[iCent]->GetFunction("pol0")->GetParameter(0), fRatio[iCent]->GetFunction("pol0")->GetParError(0)));
    p0.SetTextSize(0.035);
    chi2.Draw("same");
    p0.Draw("same");
    cRatio.Print(Form("%s/%s.png", kPlotDir, fRatio[iCent]->GetName()));
  }
  outFile.Close();
}
