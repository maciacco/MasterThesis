// Systematics.cpp
// This macro computes systematic errors on pol0 fit parameter

// TODO: update the generation of random cuts

#include <stdlib.h>
#include <string.h>

#include <TFile.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <Riostream.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveStats.h>

#include "../utils/Config.h"

using namespace proton;

// #define USE_COUNTER
const bool sys_eff_error = true;

void Systematics(const int points = kNPoints, const bool cutVar = true, const bool binCountingVar = true, const bool expVar = true, const bool sigmoidVar = true, const char *outFileName = "SystematicsAllEff")
{
  gStyle->SetOptStat(110001110);
  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.85);
  TStopwatch swatch;
  swatch.Start(true);

  TFile *specFile = TFile::Open(Form("%s/SpectraProtonSys.root", kOutDir));
  TFile *effFile = TFile::Open(Form("%s/EfficiencyProtonSys.root", kOutDir));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");
  TDirectory *cdHist = outFile->mkdir("hist");

  for (int iC = 0; iC < kNCentClasses; ++iC) // TODO: extend the analysis to the third centrality class as well
  {
    TDirectory *cdFits = outFile->mkdir(Form("fits_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    TH1D fFitPar(Form("fFitPar_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), 3000, 0.9, 1.1);
    TH1D fProb(Form("fProb_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), 1000., 0., 1.0);
    TH1D fRatio("fRatio", "fRatio", kNPtBins, kPtBins);
    TH1D fChi2Scan("fChi2Scan","fChi2Scan",300,0.9750,1.0050);
    TH1D fChi2ScanEnvelope("fChi2ScanEnvelope","fChi2ScanEnvelope",300,0.9750,1.0050);
    double minR=0.975;
    for (int iScan=0;iScan<300;++iScan){
      double R=minR+iScan*0.0001;
      fChi2ScanEnvelope.SetBinContent(iScan+1,9999999999);
    }

    // fill and fit ratio histogram nPoint times
    int iP = 0;
    while (iP < 10000)
    {
      double nUsedPtBins = 24;
      for (int iPtBin = 5; iPtBin < nUsedPtBins+1; ++iPtBin)
      {
        // extract variable which the variation is applied to
        int cutVariable = gRandom->Rndm() * 3;

        // extract cut
        int cutIndex = 0;
        if (cutVariable == 0) cutIndex = gRandom->Rndm() * kNCutDCAz;
        else if (cutVariable == 1)
        {
          cutIndex = gRandom->Rndm() * kNTPCPidSigmas;
        }
        else
        {
          cutIndex = gRandom->Rndm() * kNCutTPCClusters;
        }
        // extract background flag
        int bkgFlag = 1;// gRandom->Rndm() * 2;

        int sigmoidFlagRnd = 1;
        if (sigmoidVar)
        {
          if (fRatio.GetBinCenter(iPtBin) < 1.6)
          {
            (gRandom->Rndm() > .5) ? sigmoidFlagRnd = 1 : sigmoidFlagRnd = 0;
          }
        }

        int iNsigma = gRandom->Rndm() * 3;

        auto tmpCutSettings = cutSettings[cutVariable];
        auto tmpCutIndex = Form("%d",cutIndex);
        if ( (cutVariable==0 && cutIndex==2) || (cutVariable==1 && cutIndex==0) || (cutVariable==2 && cutIndex==2))
        {
          tmpCutIndex = Form("");
          tmpCutSettings = Form("");
        }
        if ( (cutIndex > 2 && (cutVariable == 2 || cutVariable == 0)) || (cutIndex > 0 && cutVariable == 1))
        {
          cutIndex--;

          tmpCutIndex = Form("%d",cutIndex);
        }
        auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlagRnd, iNsigma);
        //std::cout<<"___fullCutSettings: "<<fullCutSettingsSigm<<std::endl;
        TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));

        // generate systematic variation on efficiency
        TH1D *h_a_eff = (TH1D *)effFile->Get(Form("%s%s_/fAEff_TOF_%.0f_%.0f", tmpCutSettings, tmpCutIndex, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        TH1D *h_m_eff = (TH1D *)effFile->Get(Form("%s%s_/fMEff_TOF_%.0f_%.0f", tmpCutSettings, tmpCutIndex, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        double sys_a_eff = 1., sys_m_eff = 1.;
        if (sys_eff_error){
          sys_a_eff = gRandom->Gaus(0,1)*h_a_eff->GetBinError(h_a_eff->FindBin(fRatio.GetBinCenter(iPtBin)));
          double a_eff = h_a_eff->GetBinContent(h_a_eff->FindBin(fRatio.GetBinCenter(iPtBin)));
          sys_a_eff = 1./(1+sys_a_eff/a_eff);
          //std::cout<<"sys_a_eff = "<<sys_a_eff<<std::endl;
          sys_m_eff = gRandom->Gaus(0,1)*h_a_eff->GetBinError(h_m_eff->FindBin(fRatio.GetBinCenter(iPtBin)));
          double m_eff = h_m_eff->GetBinContent(h_m_eff->FindBin(fRatio.GetBinCenter(iPtBin)));
          sys_m_eff = 1./(1+sys_m_eff/m_eff);
        }

        fRatio.SetBinContent(iPtBin, h->GetBinContent(h->FindBin(fRatio.GetBinCenter(iPtBin)))*sys_a_eff/sys_m_eff);
        fRatio.SetBinError(iPtBin, h->GetBinError(h->FindBin(fRatio.GetBinCenter(iPtBin)))*sys_a_eff/sys_m_eff);
        // if(iPtBin == nUsedPtBins)
        // std::cout<<"fullCutSettings = "<<fullCutSettingsSigm<<", bin center = "<<fRatio.GetBinCenter(iPtBin)<<", content = "<< h->GetBinContent(h->FindBin(fRatio.GetBinCenter(iPtBin)))<<std::endl;
      }
      if (iP%100 == 0)
        std::cout << "p=" << iP << std::endl;
      TF1 fitFunc("fitFunc", "pol0");
      auto fit = fRatio.Fit(&fitFunc, "QS");
      double minR=0.985;
      for (int iScan=0;iScan<300;++iScan){
        double R=minR+iScan*0.0001;
        double nUsedPtBins=24;
        double chi2=0;
        for (int iPtBin = 5; iPtBin < nUsedPtBins+1; ++iPtBin){
          chi2+=(fRatio.GetBinContent(iPtBin)-R)*(fRatio.GetBinContent(iPtBin)-R)/fRatio.GetBinError(iPtBin)/fRatio.GetBinError(iPtBin);
        }
        std::cout<<"R="<<R<<"; chi2="<<chi2<<std::endl;
        fChi2Scan.SetBinContent(iScan+1,chi2);
      }

      int ndf = 19;
      if (fit->Status() == 0 /* && (fit->Chi2()/fit->Ndf()) < 3. && */ && fit->Ndf() >= ndf)
      { // check chi2
        fFitPar.Fill(fitFunc.GetParameter(0));
        fProb.Fill(fitFunc.GetProb());
        cdFits->cd();
        fRatio.Write();
        fChi2Scan.Write();
        double minR=0.975;
        for (int iScan=0;iScan<300;++iScan){
          double R=minR+iScan*0.0001;
          if(fChi2Scan.GetBinContent(iScan+1)<fChi2ScanEnvelope.GetBinContent(iScan+1)){
            fChi2ScanEnvelope.SetBinContent(iScan+1,fChi2Scan.GetBinContent(iScan+1));
          }
        }
        ++iP;
      }
      fRatio.Reset();
      fChi2Scan.Reset();
    }

    cdHist->cd();
    fFitPar.GetYaxis()->SetTitle("Entries");
    fFitPar.GetXaxis()->SetTitle("R (#bar{p}/p)");
    fFitPar.SetDrawOption("histo");
    fFitPar.Rebin(8);
    fFitPar.SetFillStyle(3345);
    fFitPar.SetLineWidth(2);
    fFitPar.SetLineColor(kBlue);
    fFitPar.SetFillColor(kBlue);
    fFitPar.Write();

    fChi2ScanEnvelope.Write();

    TCanvas cFitPar("cFitPar", "cFitPar");
    cFitPar.cd();
    cFitPar.SetTicks(1, 1);
    fFitPar.GetXaxis()->SetRangeUser(fFitPar.GetMean()-5*fFitPar.GetRMS(),fFitPar.GetMean()+5*fFitPar.GetRMS());
    fFitPar.Draw("");
    cFitPar.Print(Form("%s/Systematics_%s.pdf", kPlotDir, fFitPar.GetName()));

    fProb.GetYaxis()->SetTitle("Entries");
    fProb.GetXaxis()->SetTitle("#it{P}( #chi^{2} > #chi^{2}_{#it{obs}} )");
    fProb.Write();
  }
  outFile->Close();
  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;
}
