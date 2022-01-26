// SystematicsPt.cpp
// This macro computes systematic errors as a function of pt

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

using namespace pion;

// #define USE_COUNTER
const bool sys_eff_error = true;

void SystematicsPt(const int points = kNPoints, const bool cutVar = true, const bool binCountingVar = true, const bool expVar = true, const bool sigmoidVar = true, const char *outFileName = "SystematicsAllEPt")
{
  gStyle->SetOptStat(110001110);
  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.85);
  TStopwatch swatch;
  swatch.Start(true);

  TFile *specFile = TFile::Open(Form("%s/SpectraPionSys.root", kOutDir));
  TFile *effFile = TFile::Open(Form("%s/EfficiencyPionSys.root", kOutDir));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iC = 0; iC < kNCentClasses; ++iC) // TODO: extend the analysis to the third centrality class as well
  {
    TH1D fSystematicUncertaintyTrackCuts(Form("fSystematicUncertaintyTrackCuts_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyROI(Form("fSystematicUncertaintyROI_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyEff(Form("fSystematicUncertaintyEff_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyPrim(Form("fSystematicUncertaintyPrim_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTot(Form("fSystematicUncertaintyTot_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTotal(Form("fSystematicUncertaintyTotal_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), kNPtBins, kPtBins);
    TH1D fRatioFromVariationsTot(Form("fRatioFromVariationsTot_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]), kNPtBins, kPtBins);
    double ratioBins[1000];
    for (int iRatioBins = 0; iRatioBins < 3000; iRatioBins++){
      ratioBins[iRatioBins]=0.85+iRatioBins*0.0001;
    }
    TH2D fRatiosVsPt(Form("fRatiosVsPt_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtROI(Form("fRatiosVsPt_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtPrim(Form("fRatiosVsPt_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtTot(Form("fRatiosVsPt_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    
    // total cuts
    int cutVariable=0;
    int bkgFlag = 1;
    int sigmoidFlag = 1;
    int iNsigma = 1;
    for (int iTrackCuts=0; iTrackCuts<kNTrackCuts-4; ++iTrackCuts){
      for (int iROI=0; iROI<3; ++iROI){
        for (int iSigmoid=1; iSigmoid<2; ++iSigmoid){
          auto tmpCutSettings = trackCutSettings[iTrackCuts];
          auto cutIndex = trackCutIndexes[iTrackCuts];
          auto tmpCutIndex = Form("%d",cutIndex);
          if ( iTrackCuts == 0 )
          {
            tmpCutIndex = Form("");
            tmpCutSettings = Form("");
          }
          auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, iSigmoid, iROI);
          std::cout << fullCutSettingsSigm <<std::endl;
          TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]));
          for(int iPtBins=4;iPtBins<24;++iPtBins){
            fRatiosVsPtTot.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
          }
        }
      }
    }

    // track cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 1;
    iNsigma = 1;
    for (int iTrackCuts=0; iTrackCuts<kNTrackCuts-4; ++iTrackCuts){
      auto tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      if ( iTrackCuts == 0 )
      {
        tmpCutIndex = Form("");
        tmpCutSettings = Form("");
      }
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPt.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    // efficiency error
    TH1D *h_a_eff = (TH1D *)effFile->Get(Form("_/fAEff_TOF_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]));
    TH1D *h_m_eff = (TH1D *)effFile->Get(Form("_/fMEff_TOF_%.0f_%.0f", kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]));
    for(int iPtBin=7;iPtBin<15;++iPtBin){
      double err_eff_a=h_a_eff->GetBinError(iPtBin);
      double eff_a=h_a_eff->GetBinContent(iPtBin);
      double err_eff_m=h_m_eff->GetBinError(iPtBin);
      double eff_m=h_m_eff->GetBinContent(iPtBin);
      fSystematicUncertaintyEff.SetBinContent(iPtBin,TMath::Sqrt(err_eff_a*err_eff_a/eff_a/eff_a+err_eff_m*err_eff_m/eff_m/eff_m));
      fSystematicUncertaintyEff.SetBinError(iPtBin,0);
    }

    // roi error
    cutVariable=0;
    iNsigma = 0;
    for (int iRoi=0; iRoi<3; ++iRoi){
      auto tmpCutSettings = Form("");
      auto tmpCutIndex = Form("");
      iNsigma=iRoi;
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtROI.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    // prim error
    cutVariable=0;
    iNsigma = 1;
    for (int iPrim=0; iPrim<2; ++iPrim){
      auto tmpCutSettings = Form("dcaz");
      auto tmpCutIndex = Form("2");
      sigmoidFlag=iPrim;
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtPrim.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    for(int iPtBins=7;iPtBins<15;++iPtBins){
      // track cuts
      TH1D *proj=fRatiosVsPt.ProjectionY("py",iPtBins,iPtBins);
      fSystematicUncertaintyTrackCuts.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyTrackCuts.SetBinError(iPtBins,0);

      // roi
      proj=fRatiosVsPtROI.ProjectionY("py",iPtBins,iPtBins);
      fSystematicUncertaintyROI.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyROI.SetBinError(iPtBins,0);

      // prim
      proj=fRatiosVsPtPrim.ProjectionY("py",iPtBins,iPtBins);
      fSystematicUncertaintyPrim.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyPrim.SetBinError(iPtBins,0);

      // tot
      proj=fRatiosVsPtTot.ProjectionY("py",iPtBins,iPtBins);
      double mean = proj->GetMean();
      double rms = proj->GetRMS();
      // reject outliers
      double rejection_criterion=3.; // 3 sigma rejection
      int count_outliers = 999;
      while (count_outliers>0){
        count_outliers = 0;
        mean=proj->GetMean();
        rms=proj->GetRMS();
        for(int iB=1;iB<=proj->GetNbinsX();++iB){
          double tmp=proj->GetBinCenter(iB);
          bool isFilled=proj->GetBinContent(iB)>0;
          if (std::abs(tmp-mean)>rejection_criterion*rms && isFilled){
            proj->SetBinContent(iB,0.);
            count_outliers++;
          }
        }
        std::cout << "outliers found = " << count_outliers << std::endl;
      }
      double totSys=TMath::Sqrt(fSystematicUncertaintyEff.GetBinContent(iPtBins)*fSystematicUncertaintyEff.GetBinContent(iPtBins)+proj->GetRMS()*proj->GetRMS()/proj->GetMean()/proj->GetMean());
      fSystematicUncertaintyTot.SetBinContent(iPtBins,totSys);
      fSystematicUncertaintyTot.SetBinError(iPtBins,0);
      
      fSystematicUncertaintyTotal.SetBinContent(iPtBins,proj->GetRMS());
      fSystematicUncertaintyTotal.SetBinError(iPtBins,0);
      
      fRatioFromVariationsTot.SetBinContent(iPtBins,proj->GetMean());
      fRatioFromVariationsTot.SetBinError(iPtBins,0);
    }

    fRatiosVsPt.Write();
    fRatiosVsPtROI.Write();
    fRatiosVsPtPrim.Write();
    fRatiosVsPtTot.Write();
    fSystematicUncertaintyTrackCuts.Write();
    fSystematicUncertaintyEff.Write();
    fSystematicUncertaintyROI.Write();
    fSystematicUncertaintyPrim.Write();
    fSystematicUncertaintyTot.Write();
    fSystematicUncertaintyTotal.Write();
    fRatioFromVariationsTot.Write();
    TCanvas cSys("cSys", "cSys");
    TLegend lSys(0.154381, 0.739079, 0.460362, 0.936246);
    fSystematicUncertaintyTot.SetMinimum(0);
    fSystematicUncertaintyTot.GetYaxis()->SetTitle("Systematic Uncertainty (%)");
    fSystematicUncertaintyTot.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fSystematicUncertaintyTot.GetXaxis()->SetRangeUser(1.,2.);
    fSystematicUncertaintyTot.SetLineColor(kBlue);
    fSystematicUncertaintyTot.Draw("histo");
    lSys.AddEntry(&fSystematicUncertaintyTot,"Total");
    fSystematicUncertaintyEff.SetLineColor(kGreen);
    fSystematicUncertaintyEff.Draw("histosame");
    lSys.AddEntry(&fSystematicUncertaintyEff,"Efficiency x Acceptance correction");
    fSystematicUncertaintyTrackCuts.SetLineColor(kRed);
    fSystematicUncertaintyTrackCuts.Draw("histosame");
    lSys.AddEntry(&fSystematicUncertaintyTrackCuts,"Track cuts");
    fSystematicUncertaintyROI.SetLineColor(kOrange);
    fSystematicUncertaintyROI.Draw("histosame");
    lSys.AddEntry(&fSystematicUncertaintyROI,"Signal extraction");
    fSystematicUncertaintyPrim.SetLineColor(kBlack);
    fSystematicUncertaintyPrim.Draw("histosame");
    lSys.AddEntry(&fSystematicUncertaintyPrim,"Primary correction");
    lSys.Draw("same");
    cSys.Write();
  }
  outFile->Close();
  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;
}
