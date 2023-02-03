// SystematicsPtNotCombined.cpp

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

bool use_tpc = true;
bool kSmooth = true;
bool barlow_criterion = true;

void remove_outliers(TH1D* proj, double rejection_criterion=3.){
  double mean = proj->GetMean();
  double rms = proj->GetRMS();
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
}

// #define USE_COUNTER
const bool sys_eff_error = true;
const int used_pt_bins = 21;
const int nTrials=10000;

void SystematicsPtNotCombined(const int points = kNPoints, const bool cutVar = true, const bool binCountingVar = true, const bool expVar = true, const bool sigmoidVar = true, const char *outFileName = "SystematicsAllEPtNotCombinedTOF_extend_4")
{

  const int nUsedPtBins = 43;
  gStyle->SetTextFont(44);
  gStyle->SetOptFit(0);
  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.85);
  TStopwatch swatch;
  swatch.Start(true);

  TFile *specFile = TFile::Open(Form("%s/SpectraProtonSysTOF_extend_4.root", kOutDir));
  TFile *hijingFile = TFile::Open("../HIJINGRatios.root");
  //TFile *fileEffFit=TFile::Open("out/SpectraProton_MC21l5_raw_fitEff.root");
  TFile *fG4 = TFile::Open("out/SpectraProton_MC21l5_raw_primaryInjected.root");      
  TFile *inFileSec = TFile::Open(Form("%s/PrimaryProton_large.root", kOutDir));
  TFile *inFileSec_ = TFile::Open(Form("%s/PrimaryProton_largeTOF.root", kOutDir));
  TFile *effFile = TFile::Open(Form("%s/EfficiencyProtonMC_21l5_LOWPT_TESTTOF.root", kOutDir));//EfficiencyProtonMC_21l5_false__
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iC = 0; iC < kNCentClasses; ++iC) // TODO: extend the analysis to the third centrality class as well
  {
    gStyle->SetOptStat(1);
    TH1D fSystematicUncertaintyDCAxy(Form("fSystematicUncertaintyDCAxy_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyDCAz(Form("fSystematicUncertaintyDCAz_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTPCCls(Form("fSystematicUncertaintyTPCCls_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyPID(Form("fSystematicUncertaintyPID_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyChi2TPC(Form("fSystematicUncertaintyChi2TPC_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyROI(Form("fSystematicUncertaintyROI_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyEff(Form("fSystematicUncertaintyEff_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTFF(Form("fSystematicUncertaintyTFF_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyPrim(Form("fSystematicUncertaintyPrim_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyPrimG3G4(Form("fSystematicUncertaintyPrimG3G4_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTotal(Form("fSystematicUncertaintyTotal_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTotalPtCorrelated(Form("fSystematicUncertaintyTotalPtCorrelated_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fRatioFromVariationsTot(Form("fRatioFromVariationsTot_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    double ratioBins[3000];
    for (int iRatioBins = 0; iRatioBins < 3000; iRatioBins++){
      ratioBins[iRatioBins]=0.85+iRatioBins*0.0001;
    }
    TH2D fRatiosVsPtDCAxy(Form("fRatiosVsPtDCAxy_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtDCAz(Form("fRatiosVsPtDCAz_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtTPCCls(Form("fRatiosVsPtTPCCls_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtPID(Form("fRatiosVsPtPID_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtChi2TPC(Form("fRatiosVsPtChi2TPC_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtROI(Form("fRatiosVsPtROI_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtPrim(Form("fRatiosVsPtPrim_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtPrimG3G4(Form("fRatiosVsPtPrimG3G4_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtEffFit(Form("fRatiosVsPtEffFit_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtTot(Form("fRatiosVsPt_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    
    // total cuts
    int cutVariable=0;
    int bkgFlag = 1;
    int sigmoidFlag = 0;
    int iNsigma = 1;
    for (int iTrackCuts=0; iTrackCuts<kNTrackCuts-2; ++iTrackCuts){
      for (int iROI=0; iROI<3; ++iROI){
        for (int iG3G4Prim=0; iG3G4Prim<1; ++iG3G4Prim){ // 0 = use G3 (DO NOT USE G4 in this case)
          for (int iSigmoid=0; iSigmoid<1; ++iSigmoid){ // 0 = directly from TFF -> do not use fit function
            auto tmpCutSettings = trackCutSettings[iTrackCuts];
            auto cutIndex = trackCutIndexes[iTrackCuts];
            auto tmpCutIndex = Form("%d",cutIndex);
            if ( iTrackCuts == 0 )
            {
              tmpCutIndex = Form("");
              tmpCutSettings = Form("");
            }
            auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, iSigmoid, iROI, iG3G4Prim);
            std::cout << fullCutSettingsSigm <<std::endl;
            TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
            for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
              fRatiosVsPtTot.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
            }
          }
        }
      }
    }

    TH1D *hDefault = (TH1D *)specFile->Get(Form("_1_0_1_0/fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));

    // DCAz cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 0;
    iNsigma = 1;
    for (int iTrackCuts=1; iTrackCuts<5; ++iTrackCuts){
      auto tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      if ( iTrackCuts == 0 )
      {
        tmpCutIndex = Form("");
        tmpCutSettings = Form("");
      }
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion)
          fRatiosVsPtDCAz.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    // PID cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 0;
    iNsigma = 1;
    for (int iTrackCuts=5; iTrackCuts<7; ++iTrackCuts){
      auto tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      if ( iTrackCuts == 0 )
      {
        tmpCutIndex = Form("");
        tmpCutSettings = Form("");
      }
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion)
          fRatiosVsPtPID.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    // TPCCls cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 0;
    iNsigma = 1;
    for (int iTrackCuts=7; iTrackCuts<11; ++iTrackCuts){
      auto tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      if ( iTrackCuts == 0 )
      {
        tmpCutIndex = Form("");
        tmpCutSettings = Form("");
      }
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion) fRatiosVsPtTPCCls.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    // DCAxy cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 0;
    iNsigma = 1;
    for (int iTrackCuts=11; iTrackCuts<kNTrackCuts-2-2; ++iTrackCuts){
      auto tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      if ( iTrackCuts == 0 )
      {
        tmpCutIndex = Form("");
        tmpCutSettings = Form("");
      }
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion) fRatiosVsPtDCAxy.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    // chi2tpc cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 0;
    iNsigma = 1;
    for (int iTrackCuts=15; iTrackCuts<kNTrackCuts-2; ++iTrackCuts){
      auto tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      if ( iTrackCuts == 0 )
      {
        tmpCutIndex = Form("");
        tmpCutSettings = Form("");
      }
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion) fRatiosVsPtChi2TPC.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    // efficiency error
    TH1D *h_a_eff = (TH1D *)effFile->Get(Form("_/fAEff_TOF_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    TH1D *h_m_eff = (TH1D *)effFile->Get(Form("_/fMEff_TOF_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    for(int iPtBin=5;iPtBin<nUsedPtBins;++iPtBin){
      double err_eff_a=h_a_eff->GetBinError(iPtBin);
      double eff_a=h_a_eff->GetBinContent(iPtBin);
      double err_eff_m=h_m_eff->GetBinError(iPtBin);
      double eff_m=h_m_eff->GetBinContent(iPtBin);
      fSystematicUncertaintyEff.SetBinContent(iPtBin,TMath::Sqrt(err_eff_a*err_eff_a/eff_a/eff_a+err_eff_m*err_eff_m/eff_m/eff_m));
      fSystematicUncertaintyEff.SetBinError(iPtBin,0);
    }

    // primary fraction (TFF) error
    for(int iPtBin=5;iPtBin<nUsedPtBins;++iPtBin){
      double primaryRelativeError[2];
      for (int iMatt = 0; iMatt < 2; ++iMatt){ // take uncertainties directly from TFF output
        TH1D* h_sec = (TH1D*)inFileSec_->Get(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        //if (h_sec->GetBinCenter(iPtBin)<1.) h_sec = (TH1D*)inFileSec_->Get(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        if (h_sec->GetBinCenter(iPtBin)>0.8){
          double primaryError=h_sec->GetBinError(iPtBin);
          double primary=h_sec->GetBinContent(iPtBin);
          primaryRelativeError[iMatt]=primaryError/primary;
        }
        else 
          primaryRelativeError[iMatt]=0;
      }
      fSystematicUncertaintyTFF.SetBinContent(iPtBin,TMath::Sqrt(primaryRelativeError[0]*primaryRelativeError[0]+primaryRelativeError[1]*primaryRelativeError[1]));
      fSystematicUncertaintyTFF.SetBinError(iPtBin,0);
    }

    // roi error
    cutVariable=0;
    iNsigma = 0;
    for (int iRoi=0; iRoi<3; ++iRoi){
      auto tmpCutSettings = Form("");
      auto tmpCutIndex = Form("");
      iNsigma=iRoi;
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion) fRatiosVsPtROI.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    // prim error
    cutVariable=0;
    iNsigma = 1;
    for (int iPrim=0; iPrim<1; ++iPrim){
      auto tmpCutSettings = Form("");
      auto tmpCutIndex = Form("");
      sigmoidFlag=iPrim;
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        fRatiosVsPtPrim.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    // prim error - G3 vs. G4
    cutVariable=0;
    iNsigma = 1;
    for (int iPrim=0; iPrim<1; ++iPrim){
      auto tmpCutSettings = Form("");
      auto tmpCutIndex = Form("");
      sigmoidFlag=iPrim;
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma, iPrim);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        fRatiosVsPtPrimG3G4.Fill(kPtBins[iPtBins-1],h->GetBinContent(iPtBins));
      }
    }

    for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
      // track cuts
      TH1D *proj=fRatiosVsPtDCAz.ProjectionY("py",iPtBins,iPtBins);
      if (proj->GetEntries()==0 || proj->GetMean()==0){
        fSystematicUncertaintyDCAz.SetBinContent(iPtBins,0);
      }
      else{
        fSystematicUncertaintyDCAz.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      }
      fSystematicUncertaintyDCAz.SetBinError(iPtBins,0);

      proj=fRatiosVsPtPID.ProjectionY("py",iPtBins,iPtBins);
      if (proj->GetEntries()==0 || proj->GetMean()==0){
        fSystematicUncertaintyPID.SetBinContent(iPtBins,0);
      }
      else{
        fSystematicUncertaintyPID.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      }fSystematicUncertaintyPID.SetBinError(iPtBins,0);

      proj=fRatiosVsPtTPCCls.ProjectionY("py",iPtBins,iPtBins);
      if (proj->GetEntries()==0 || proj->GetMean()==0){
        fSystematicUncertaintyTPCCls.SetBinContent(iPtBins,0);
      }
      else{
        fSystematicUncertaintyTPCCls.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      }
      fSystematicUncertaintyTPCCls.SetBinError(iPtBins,0);

      proj=fRatiosVsPtDCAxy.ProjectionY("py",iPtBins,iPtBins);
      if (proj->GetEntries()==0 || proj->GetMean()==0){
        fSystematicUncertaintyDCAxy.SetBinContent(iPtBins,0);
      }
      else{
        fSystematicUncertaintyDCAxy.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      }
      fSystematicUncertaintyDCAxy.SetBinError(iPtBins,0);

      proj=fRatiosVsPtChi2TPC.ProjectionY("py",iPtBins,iPtBins);
      if (proj->GetEntries()==0 || proj->GetMean()==0){
        fSystematicUncertaintyChi2TPC.SetBinContent(iPtBins,0);
      }
      else
        fSystematicUncertaintyChi2TPC.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyChi2TPC.SetBinError(iPtBins,0);

      // roi
      proj=fRatiosVsPtROI.ProjectionY("py",iPtBins,iPtBins);
      if (proj->GetEntries()==0 || proj->GetMean()==0){
        fSystematicUncertaintyROI.SetBinContent(iPtBins,0);
      }
      else
        fSystematicUncertaintyROI.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyROI.SetBinError(iPtBins,0);

      // prim
      proj=fRatiosVsPtPrim.ProjectionY("py",iPtBins,iPtBins);
      fSystematicUncertaintyPrim.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyPrim.SetBinError(iPtBins,0);
      
      // prim G3 vs. G4
      proj=fRatiosVsPtPrimG3G4.ProjectionY("py",iPtBins,iPtBins);
      fSystematicUncertaintyPrimG3G4.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyPrimG3G4.SetBinError(iPtBins,0);

      // tot
      proj=fRatiosVsPtTot.ProjectionY("py",iPtBins,iPtBins);
      remove_outliers(proj);
      
      fRatioFromVariationsTot.SetBinContent(iPtBins,proj->GetMean());
      fRatioFromVariationsTot.SetBinError(iPtBins,0);
    }

    outFile->cd();
    TCanvas cVarTot(fRatiosVsPtTot.GetName(),fRatiosVsPtTot.GetTitle());
    fRatiosVsPtTot.Draw("colz");
    fRatiosVsPtTot.GetXaxis()->SetRangeUser(1.,3.);
    fRatiosVsPtTot.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fRatiosVsPtTot.GetYaxis()->SetTitle("Ratio #bar{p}/p");
    cVarTot.Print(Form("%s.pdf",fRatiosVsPtTot.GetName()));
    fRatioFromVariationsTot.Write();
    fRatiosVsPtDCAxy.Write();
    fRatiosVsPtTot.Write();
    fRatiosVsPtDCAz.Write();
    fRatiosVsPtPID.Write();
    fRatiosVsPtTPCCls.Write();
    fRatiosVsPtChi2TPC.Write();
    fRatiosVsPtPrim.Write();
    fRatiosVsPtPrimG3G4.Write();
    if(kSmooth)
    {
      if (!use_tpc){
        fSystematicUncertaintyDCAxy.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyDCAz.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyPID.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyTPCCls.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyChi2TPC.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyPrim.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyROI.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyEff.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyTFF.GetXaxis()->SetRangeUser(1.,1.2);
        fSystematicUncertaintyPrimG3G4.GetXaxis()->SetRangeUser(1.,1.2);
      }
      else {
        fSystematicUncertaintyDCAxy.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyDCAz.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyPID.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyTPCCls.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyChi2TPC.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyPrim.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyROI.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyEff.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyTFF.GetXaxis()->SetRangeUser(1.2,2.);
        fSystematicUncertaintyTFF.Smooth(1,"R");
        fSystematicUncertaintyTFF.GetXaxis()->SetRangeUser(2.,2.6);
        fSystematicUncertaintyTFF.Smooth(1,"R");
        fSystematicUncertaintyTFF.GetXaxis()->SetRangeUser(1.2,3.);
        fSystematicUncertaintyPrimG3G4.GetXaxis()->SetRangeUser(1.2,3.);
      }
      fSystematicUncertaintyDCAxy.Smooth(1,"R");
      fSystematicUncertaintyDCAz.Smooth(1,"R");
      fSystematicUncertaintyPID.Smooth(1,"R");
      fSystematicUncertaintyTPCCls.Smooth(1,"R");
      fSystematicUncertaintyChi2TPC.Smooth(1,"R");
      fSystematicUncertaintyPrim.Smooth(1,"R");
      fSystematicUncertaintyROI.Smooth(1,"R");
      //fSystematicUncertaintyEff.Smooth(1,"R");
      //fSystematicUncertaintyTFF.Smooth(1,"R");
      fSystematicUncertaintyPrimG3G4.Smooth(1,"R");
    }
    fSystematicUncertaintyDCAxy.Write();
    fSystematicUncertaintyDCAz.Write();
    fSystematicUncertaintyPID.Write();
    fSystematicUncertaintyTPCCls.Write();
    fSystematicUncertaintyChi2TPC.Write();
    fSystematicUncertaintyPrim.Write();
    fSystematicUncertaintyROI.Write();
    fSystematicUncertaintyEff.Write();
    fSystematicUncertaintyTFF.Write();
    fSystematicUncertaintyPrimG3G4.Write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * //
    //   C O M P U T E   P T   C O R R E L A T I O N S   //
    // * * * * * * * * * * * * * * * * * * * * * * * * * //
    TH2D fSystematicsCorrelationsDCAxy(Form("fSystematicsCorrelationsDCAxy_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),5,0.095,0.145,1000,-.04,.04);
    TH2D fSystematicsCorrelationsDCAz(Form("fSystematicsCorrelationsDCAz_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),31,0.475,2.025,1000,-.04,.04);
    TH2D fSystematicsCorrelationsTPCCls(Form("fSystematicsCorrelationsTPCCls_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),24,58.5,82.5,1000,-.04,.04);
    TH2D fSystematicsCorrelationsChi2TPC(Form("fSystematicsCorrelationsChi2TPC_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),26,1.995,2.255,1000,-.04,.04);
    TH2D fSystematicsCorrelationsPID(Form("fSystematicsCorrelationsPID_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),25,3.245,3.505,1000,-.04,.04);
    TH2D fSystematicsCorrelationsROI(Form("fSystematicsCorrelationsROI_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),51,7.49,8.51,1000,-.04,.04);

    // roi
    cutVariable = 0;
    bkgFlag = 1;
    sigmoidFlag = 0;
    iNsigma = 1;
    hDefault = (TH1D *)specFile->Get(Form("_1_0_1_0/fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    for (int iROI=0;iROI<3;++iROI){
      if (iROI==1)continue;
      TString tmpCutSettings = Form("");
      auto cutIndex = 0;
      auto tmpCutIndex = Form("",cutIndex);
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings.Data(), tmpCutIndex, bkgFlag, sigmoidFlag, iROI);
      //std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion) fSystematicsCorrelationsROI.Fill(kROIVariation[iROI],h->GetBinContent(iPtBins)-hDefault->GetBinContent(iPtBins));
      }
    }
        
    // track cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 0;
    iNsigma = 1;
    hDefault = (TH1D *)specFile->Get(Form("_1_0_1_0/fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    for (int iTrackCuts=1; iTrackCuts<kNTrackCuts-2; ++iTrackCuts){
      TString tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_0", tmpCutSettings.Data(), tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      //std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
        double h_content = h->GetBinContent(iPtBins);
        double h_error = h->GetBinError(iPtBins);
        double h_0_content = hDefault->GetBinContent(iPtBins);
        double h_0_error = hDefault->GetBinError(iPtBins);
        double z = (h_content-h_error)/TMath::Sqrt(h_0_error*h_0_error+h_error*h_error);
        if ((std::abs(z)>2. && barlow_criterion) || !barlow_criterion){
          if (tmpCutSettings.EqualTo("dcaz"))
          {
            int cutVal = cutIndex;
            if (cutIndex > 1) cutVal = cutIndex+1;
            fSystematicsCorrelationsDCAz.Fill(kCutDCAz[cutVal],h->GetBinContent(iPtBins)-hDefault->GetBinContent(iPtBins));
          }
          if (tmpCutSettings.EqualTo("dcaxy"))
          {
            int cutVal = cutIndex;
            if (cutIndex > 1) cutVal = cutIndex+1;
            fSystematicsCorrelationsDCAxy.Fill(kCutDCAxy[cutVal],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
          }
          if (tmpCutSettings.EqualTo("tpc"))
          {
            int cutVal = cutIndex;
            if (cutIndex > 1) cutVal = cutIndex+1;
            fSystematicsCorrelationsTPCCls.Fill(kCutTPCClusters[cutVal],h->GetBinContent(iPtBins)-hDefault->GetBinContent(iPtBins));
          }
          if (tmpCutSettings.EqualTo("pid"))
          {
            int cutVal = cutIndex+1;
            fSystematicsCorrelationsPID.Fill(kTPCPidSigmas[cutVal],h->GetBinContent(iPtBins)-hDefault->GetBinContent(iPtBins));
          }
          if (tmpCutSettings.EqualTo("chisquare"))
          {
            int cutVal = cutIndex;
            fSystematicsCorrelationsChi2TPC.Fill(kCutChi2TPCVariations[cutVal],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
          }
        }
      }
    }

    std::cout << "* * * * * * * * * *"<<std::endl;
    std::cout << "Centrality = " << iC <<std::endl;
    double correlationDCAz = fSystematicsCorrelationsDCAz.GetCorrelationFactor();
    double covDCAz = fSystematicsCorrelationsDCAz.GetCovariance();
    double xStdDevDCAz = fSystematicsCorrelationsDCAz.GetRMS(1);
    double yStdDevDCAz = fSystematicsCorrelationsDCAz.GetRMS(2);
    std::cout << "CorrelationDCAz = " << correlationDCAz << std::endl;// "; mine = " << covDCAz/xStdDevDCAz/yStdDevDCAz <<std::endl;
    double correlationDCAxy = fSystematicsCorrelationsDCAxy.GetCorrelationFactor();
    std::cout << "CorrelationDCAxy = " << correlationDCAxy <<std::endl;
    double correlationTPCCls = fSystematicsCorrelationsTPCCls.GetCorrelationFactor();
    std::cout << "CorrelationTPCCls = " << correlationTPCCls <<std::endl;
    double correlationPID = fSystematicsCorrelationsPID.GetCorrelationFactor();
    std::cout << "CorrelationPID = " << correlationPID <<std::endl;
    double correlationChi2TPC = fSystematicsCorrelationsChi2TPC.GetCorrelationFactor();
    std::cout << "CorrelationChi2TPC = " << correlationChi2TPC <<std::endl;
    std::cout << "* * * * * * * * * *"<<std::endl;
    double correlationROI = fSystematicsCorrelationsROI.GetCorrelationFactor();
    std::cout << "CorrelationROI = " << correlationROI <<std::endl;
    std::cout << "* * * * * * * * * *"<<std::endl;

    TCanvas cCorrDCAz(fSystematicsCorrelationsDCAz.GetName(),fSystematicsCorrelationsDCAz.GetTitle());
    fSystematicsCorrelationsDCAz.GetXaxis()->SetTitle("DCA_{z} cut (cm)");
    fSystematicsCorrelationsDCAz.GetYaxis()->SetTitle("#it{p}_{T}-differential ratios - default");
    fSystematicsCorrelationsDCAz.GetYaxis()->SetRangeUser(-0.1,0.1);
    fSystematicsCorrelationsDCAz.Draw("colz");
    cCorrDCAz.Print(Form("%s.pdf",cCorrDCAz.GetName()));
    fSystematicsCorrelationsDCAz.Write();
    fSystematicsCorrelationsDCAxy.Write();

    TCanvas cCorrTPCCls(fSystematicsCorrelationsTPCCls.GetName(),fSystematicsCorrelationsTPCCls.GetTitle());
    fSystematicsCorrelationsTPCCls.GetXaxis()->SetTitle("n_{TPCcluster} cut");
    fSystematicsCorrelationsTPCCls.GetYaxis()->SetTitle("#it{p}_{T}-differential ratios - default");
    fSystematicsCorrelationsTPCCls.GetYaxis()->SetRangeUser(-0.05,0.05);
    fSystematicsCorrelationsTPCCls.Draw("colz");
    cCorrTPCCls.Print(Form("%s.pdf",cCorrTPCCls.GetName()));
    fSystematicsCorrelationsTPCCls.Write();
    
    TCanvas cCorrPID(fSystematicsCorrelationsPID.GetName(),fSystematicsCorrelationsPID.GetTitle());
    fSystematicsCorrelationsPID.GetXaxis()->SetTitle("TPC PID cut (#sigma_{TPC})");
    fSystematicsCorrelationsPID.GetYaxis()->SetTitle("#it{p}_{T}-differential ratios - default");
    fSystematicsCorrelationsPID.GetYaxis()->SetRangeUser(-0.03,0.03);
    fSystematicsCorrelationsPID.Draw("colz");
    cCorrPID.Print(Form("%s.pdf",cCorrPID.GetName()));
    fSystematicsCorrelationsPID.Write();

    TCanvas cCorrROI(fSystematicsCorrelationsROI.GetName(),fSystematicsCorrelationsROI.GetTitle());
    fSystematicsCorrelationsROI.GetXaxis()->SetTitle("ROI width (#sigma_{TOF})");
    fSystematicsCorrelationsROI.GetYaxis()->SetTitle("#it{p}_{T}-differential ratios - default");
    fSystematicsCorrelationsROI.GetYaxis()->SetRangeUser(-0.05,0.05);
    fSystematicsCorrelationsROI.Draw("colz");
    cCorrROI.Print(Form("%s.pdf",cCorrROI.GetName()));
    fSystematicsCorrelationsROI.Write();

    TCanvas cCorrDCAxy(fSystematicsCorrelationsDCAxy.GetName(),fSystematicsCorrelationsDCAxy.GetTitle());
    fSystematicsCorrelationsDCAxy.GetXaxis()->SetTitle("DCA_{xy}cut (cm)");
    fSystematicsCorrelationsDCAxy.GetYaxis()->SetTitle("#it{p}_{T}-differential ratios - default");
    fSystematicsCorrelationsDCAxy.GetYaxis()->SetRangeUser(-0.05,0.05);
    fSystematicsCorrelationsDCAxy.Draw("colz");
    cCorrDCAxy.Print(Form("%s.pdf",cCorrDCAxy.GetName()));
    fSystematicsCorrelationsDCAxy.Write();

    TCanvas cCorrChi2TPC(fSystematicsCorrelationsChi2TPC.GetName(),fSystematicsCorrelationsChi2TPC.GetTitle());
    fSystematicsCorrelationsChi2TPC.GetXaxis()->SetTitle("#chi^{2}_{TPC} cut");
    fSystematicsCorrelationsChi2TPC.GetYaxis()->SetTitle("#it{p}_{T}-differential ratios - default");
    fSystematicsCorrelationsChi2TPC.GetYaxis()->SetRangeUser(-0.05,0.05);
    fSystematicsCorrelationsChi2TPC.Draw("colz");
    cCorrChi2TPC.Print(Form("%s.pdf",cCorrChi2TPC.GetName()));
    fSystematicsCorrelationsChi2TPC.Write();

    fSystematicsCorrelationsDCAxy.Write();
    fSystematicsCorrelationsChi2TPC.Write();
    
    for(int iPtBins=5;iPtBins<nUsedPtBins;++iPtBins){
      // sigmas
      double sigma_DCAz=fSystematicUncertaintyDCAz.GetBinContent(iPtBins);
      double sigma_PID=fSystematicUncertaintyPID.GetBinContent(iPtBins);
      double sigma_TPCCls=fSystematicUncertaintyTPCCls.GetBinContent(iPtBins);
      double sigma_DCAxy=fSystematicUncertaintyDCAxy.GetBinContent(iPtBins);
      double sigma_Chi2TPC=fSystematicUncertaintyChi2TPC.GetBinContent(iPtBins);
      double sigma_ROI=fSystematicUncertaintyROI.GetBinContent(iPtBins);
      double sigma_Prim=fSystematicUncertaintyPrim.GetBinContent(iPtBins);
      double sigma_Prim_G3G4=fSystematicUncertaintyPrimG3G4.GetBinContent(iPtBins);
      double sigma_Eff=fSystematicUncertaintyEff.GetBinContent(iPtBins);
      double sigma_TFF=fSystematicUncertaintyTFF.GetBinContent(iPtBins);

      // fractions (uncorrelated)
      double fraction_uncorr_DCAz=1.-std::abs(correlationDCAz);
      double fraction_uncorr_PID=1.-std::abs(correlationPID);
      double fraction_uncorr_TPCCls=1.-std::abs(correlationTPCCls);
      double fraction_uncorr_ROI=1.-std::abs(correlationROI);
      double fraction_uncorr_DCAxy=1.-std::abs(correlationDCAxy);
      double fraction_uncorr_Chi2TPC=1.-std::abs(correlationChi2TPC);
      double totSys=0;
      double totSysPtCorrelated=0;
      if (iC>0){
        totSys=TMath::Sqrt(sigma_DCAz*sigma_DCAz*fraction_uncorr_DCAz*fraction_uncorr_DCAz
        +sigma_PID*sigma_PID*fraction_uncorr_PID*fraction_uncorr_PID
        +sigma_TPCCls*sigma_TPCCls*fraction_uncorr_TPCCls*fraction_uncorr_TPCCls
        +sigma_DCAxy*sigma_DCAxy*fraction_uncorr_DCAxy*fraction_uncorr_DCAxy
        +sigma_ROI*sigma_ROI*fraction_uncorr_ROI*fraction_uncorr_ROI
        +sigma_Chi2TPC*sigma_Chi2TPC*fraction_uncorr_Chi2TPC*fraction_uncorr_Chi2TPC
        //+sigma_Prim*sigma_Prim // -> do not include variation of fraction of primaries
        +sigma_Eff*sigma_Eff
        +sigma_TFF*sigma_TFF);
        fSystematicUncertaintyTotal.SetBinContent(iPtBins,totSys);
        fSystematicUncertaintyTotal.SetBinError(iPtBins,0);

        totSysPtCorrelated=TMath::Sqrt(sigma_DCAz*sigma_DCAz*correlationDCAz*correlationDCAz
        +sigma_PID*sigma_PID*correlationPID*correlationPID
        +sigma_TPCCls*sigma_TPCCls*correlationTPCCls*correlationTPCCls
        +sigma_ROI*sigma_ROI*correlationROI*correlationROI
        +sigma_Chi2TPC*sigma_Chi2TPC*correlationChi2TPC*correlationChi2TPC
        +sigma_DCAxy*sigma_DCAxy*correlationDCAxy*correlationDCAxy
        );
      }

      else {
        totSys=TMath::Sqrt(sigma_DCAz*sigma_DCAz*fraction_uncorr_DCAz*fraction_uncorr_DCAz
        //+sigma_PID*sigma_PID*fraction_uncorr_PID*fraction_uncorr_PID
        +sigma_TPCCls*sigma_TPCCls*fraction_uncorr_TPCCls*fraction_uncorr_TPCCls
        +sigma_DCAxy*sigma_DCAxy*fraction_uncorr_DCAxy*fraction_uncorr_DCAxy
        +sigma_ROI*sigma_ROI*fraction_uncorr_ROI*fraction_uncorr_ROI
        +sigma_Chi2TPC*sigma_Chi2TPC*fraction_uncorr_Chi2TPC*fraction_uncorr_Chi2TPC
        +sigma_Eff*sigma_Eff
        +sigma_TFF*sigma_TFF);
        fSystematicUncertaintyTotal.SetBinContent(iPtBins,totSys);
        fSystematicUncertaintyTotal.SetBinError(iPtBins,0);

        totSysPtCorrelated=TMath::Sqrt(sigma_DCAz*sigma_DCAz*correlationDCAz*correlationDCAz
        //+sigma_PID*sigma_PID*correlationPID*correlationPID
        +sigma_TPCCls*sigma_TPCCls*correlationTPCCls*correlationTPCCls
        +sigma_ROI*sigma_ROI*correlationROI*correlationROI
        +sigma_DCAxy*sigma_DCAxy*correlationDCAxy*correlationDCAxy
        +sigma_Chi2TPC*sigma_Chi2TPC*correlationChi2TPC*correlationChi2TPC
        );
      }
      fSystematicUncertaintyTotalPtCorrelated.SetBinContent(iPtBins,totSysPtCorrelated);
      fSystematicUncertaintyTotalPtCorrelated.SetBinError(iPtBins,0);
    }

    // Compute pt correlated systematic uncertainty
    TH1D hRatio(Form("fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins);
    for (int iPtBin=5;iPtBin<nUsedPtBins;++iPtBin){
      if (hRatio.GetBinCenter(iPtBin)<0.99 && use_tpc){
        TFile inFileTPC("out/SystematicsAllEPtNotCombinedTPC_extend_4.root");
        auto htmp_sys=(TH1D*)inFileTPC.Get(Form("fSystematicUncertaintyTotal_%.0f_%.0f;1", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        auto htmp=(TH1D*)inFileTPC.Get(Form("fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        fSystematicUncertaintyTotal.SetBinContent(iPtBin,htmp_sys->GetBinContent(iPtBin));
        fSystematicUncertaintyTotal.SetBinError(iPtBin,0.);
        hRatio.SetBinContent(iPtBin,htmp->GetBinContent(iPtBin));
        hRatio.SetBinError(iPtBin,htmp->GetBinError(iPtBin));
        inFileTPC.Close();
      }
      else if (hRatio.GetBinCenter(iPtBin)>0.99 && hRatio.GetBinCenter(iPtBin)<1.19 && use_tpc){
        TFile inFileTPC("out/SystematicsAllEPtNotCombinedTPC_extend_4TOF.root");
        auto htmp_sys=(TH1D*)inFileTPC.Get(Form("fSystematicUncertaintyTotal_%.0f_%.0f;1", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        auto htmp=(TH1D*)inFileTPC.Get(Form("fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        fSystematicUncertaintyTotal.SetBinContent(iPtBin,htmp_sys->GetBinContent(iPtBin));
        fSystematicUncertaintyTotal.SetBinError(iPtBin,0.);
        hRatio.SetBinContent(iPtBin,htmp->GetBinContent(iPtBin));
        hRatio.SetBinError(iPtBin,htmp->GetBinError(iPtBin));
        inFileTPC.Close();
      }
      else{
        //double scalingFactor = scaling_factor_antip(hRatio.GetBinCenter(iPtBin));
        hRatio.SetBinContent(iPtBin,fRatioFromVariationsTot.GetBinContent(iPtBin));
        hRatio.SetBinError(iPtBin,hRatio.GetBinContent(iPtBin)*fSystematicUncertaintyTotal.GetBinContent(iPtBin));
      }
    }
    outFile->cd();
    fSystematicUncertaintyTotal.Write();
    hRatio.Fit("pol0","R","",.5,3.);
    hRatio.Write();

    for (int iPtBin=3;iPtBin<nUsedPtBins;++iPtBin){
      if (hRatio.GetBinCenter(iPtBin)<0.99 && use_tpc){
        TFile inFileTPC("out/SystematicsAllEPtNotCombinedTPC_extend_4.root");
        auto htmp=(TH1D*)inFileTPC.Get(Form("fSystematicUncertaintyTotalPtCorrelated_%.0f_%.0f;1", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        fSystematicUncertaintyTotalPtCorrelated.SetBinContent(iPtBin,htmp->GetBinContent(iPtBin));
        fSystematicUncertaintyTotalPtCorrelated.SetBinError(iPtBin,0.);
        inFileTPC.Close();
      }
      else if (hRatio.GetBinCenter(iPtBin)>0.99 && hRatio.GetBinCenter(iPtBin)<1.19 && use_tpc){

        TFile inFileTPC("out/SystematicsAllEPtNotCombinedTPC_extend_4TOF.root");
        auto htmp=(TH1D*)inFileTPC.Get(Form("fSystematicUncertaintyTotalPtCorrelated_%.0f_%.0f;1", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        fSystematicUncertaintyTotalPtCorrelated.SetBinContent(iPtBin,htmp->GetBinContent(iPtBin));
        fSystematicUncertaintyTotalPtCorrelated.SetBinError(iPtBin,0.);
        inFileTPC.Close();
      }
    }

    std::cout<<"; err(0.525) = "<<fSystematicUncertaintyTotalPtCorrelated.GetBinContent(5)<<std::endl;
    outFile->cd();
    fSystematicUncertaintyTotalPtCorrelated.Write();

    TH1D h_trial("h_trial","h_trial",kNPtBins,kPtBins);
    TH1D fRatioDistributionTrials(Form("fRatioDistributionTrials_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),2000,0.95,1.05);
    for (int iTrial=0;iTrial<nTrials;++iTrial){
      double nsigma=gRandom->Gaus(0,1);
      for (int iPtBin=5;iPtBin<nUsedPtBins;++iPtBin){
        //std::cout<<"ratios from var = "<<hRatio.GetBinContent(iPtBin)<<std::endl;
        h_trial.SetBinContent(iPtBin,hRatio.GetBinContent(iPtBin)+nsigma*fSystematicUncertaintyTotalPtCorrelated.GetBinContent(iPtBin)*hRatio.GetBinContent(iPtBin));
        h_trial.SetBinError(iPtBin,hRatio.GetBinContent(iPtBin)*fSystematicUncertaintyTotal.GetBinContent(iPtBin));
      }
      //if (h_trial.GetEntries()) continue;
      h_trial.Fit("pol0","QR","",0.5,3.);
      //std::cout<<"trial: "<<iTrial<<std::endl;
      outFile->cd();
      //h_trial.Write();
      //if( /* h_trial.GetFunction("pol0")->GetProb()>0.025 &&  */h_trial.GetFunction("pol0")->GetProb()<0.95/* (h_trial.GetFunction("pol0")->GetChisquare()/h_trial.GetFunction("pol0")->GetNDF())<2. */)
      fRatioDistributionTrials.Fill(h_trial.GetFunction("pol0")->GetParameter(0));
    }
    
    for (int iPtBin=5;iPtBin<nUsedPtBins;++iPtBin){  
      fSystematicUncertaintyTotalPtCorrelated.SetBinContent(iPtBin,fSystematicUncertaintyTotalPtCorrelated.GetBinContent(iPtBin)*hRatio.GetBinContent(iPtBin));
    }

    TCanvas cRatio(Form("cRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), "cRatio");
    cRatio.SetTicks(1, 1);
    hRatio.SetMarkerStyle(20);
    hRatio.SetMarkerSize(0.8);
    hRatio.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio.GetYaxis()->SetTitle("Ratio #bar{p}/p");
    hRatio.GetXaxis()->SetRangeUser(.5,3.);
    hRatio.GetYaxis()->SetRangeUser(0.9, 1.1);
    hRatio.SetStats(0);
    hRatio.Draw("pe2");
    hRatio.SetLineColor(centrality_colors[iC]);
    hRatio.SetMarkerColor(centrality_colors[iC]);
    TGraphErrors gRatio(&hRatio);
    for (int iPtBin=3;iPtBin<nUsedPtBins;++iPtBin){  
      hRatio.SetBinError(iPtBin,fSystematicUncertaintyTotalPtCorrelated.GetBinContent(iPtBin)*fRatioFromVariationsTot.GetBinContent(iPtBin));
    }
    TGraphErrors gRatioCorr(&hRatio);
    //gRatio.Draw("P5 same");
    gRatio.GetXaxis()->SetRangeUser(.5,3.);
    gRatio.GetYaxis()->SetRangeUser(0.9, 1.1);
    gRatioCorr.SetFillStyle(3145);
    gRatioCorr.SetFillColor(centrality_colors[iC]);
    gRatio.Draw("e5 same");
    gRatioCorr.Draw("P5 same");
    hRatio.GetFunction("pol0")->Draw("same");
    TLatex chi2(2., 1.04, Form("#chi^{2}/NDF = %.2f/%d", hRatio.GetFunction("pol0")->GetChisquare(), hRatio.GetFunction("pol0")->GetNDF()));
    chi2.SetTextSize(28);
    TLatex p0(2., 1.07, Form("R = %.4f #pm %.4f", hRatio.GetFunction("pol0")->GetParameter(0), hRatio.GetFunction("pol0")->GetParError(0)));
    p0.SetTextSize(28);
    chi2.Draw("same");
    p0.Draw("same");
    cRatio.Modified();
    cRatio.Print(Form("%s/%s.pdf", kPlotDir, hRatio.GetName()));


    TCanvas cPtCorrelatedError(fRatioDistributionTrials.GetName(),fRatioDistributionTrials.GetTitle());
    fRatioDistributionTrials.GetXaxis()->SetRangeUser(fRatioDistributionTrials.GetMean()-5*fRatioDistributionTrials.GetRMS(),fRatioDistributionTrials.GetMean()+5*fRatioDistributionTrials.GetRMS());
    fRatioDistributionTrials.GetXaxis()->SetTitle("R(#bar{p}/p)");
    fRatioDistributionTrials.GetYaxis()->SetTitle("Entries");
    fRatioDistributionTrials.Draw("");
    gPad->Update();
    TPaveStats* paveStat = (TPaveStats *)fRatioDistributionTrials.FindObject("stats");
    paveStat->SetOptStat(110001110);
    cPtCorrelatedError.Modified();
    cPtCorrelatedError.Print(Form("plots/%s.pdf",fRatioDistributionTrials.GetName()));
    fRatioDistributionTrials.Write();
    //std::cout<<"Correlations fitEff = "<<fSystematicsCorrelationsEffFit.GetCorrelationFactor()<<std::endl;

    // save ratio plots
    TLegend lSys(0.2,0.6,0.4,0.8);
    fSystematicUncertaintyTotal.GetXaxis()->SetRangeUser(0.5,3.);
    fSystematicUncertaintyTotal.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fSystematicUncertaintyTotal.GetYaxis()->SetTitle("Systematic Uncertainty");
    fSystematicUncertaintyTotal.SetMinimum(0.);
    fSystematicUncertaintyTotal.GetYaxis()->SetRangeUser(0.,0.04);
    TCanvas cSysError(fSystematicUncertaintyTotal.GetName(),fSystematicUncertaintyTotal.GetTitle());
    fSystematicUncertaintyTotal.SetLineWidth(2);
    fSystematicUncertaintyTotal.Draw("histo");
    fSystematicUncertaintyTotalPtCorrelated.SetLineColor(kRed);
    fSystematicUncertaintyTotalPtCorrelated.SetLineWidth(2);
    fSystematicUncertaintyTotalPtCorrelated.Draw("histo same");
    lSys.AddEntry(&fSystematicUncertaintyTotal,"#it{p}_{T}-uncorrelated");
    lSys.AddEntry(&fSystematicUncertaintyTotalPtCorrelated,"#it{p}_{T}-correlated");
    lSys.Draw("same");
    gPad->Update();
    paveStat = (TPaveStats *)fSystematicUncertaintyTotal.FindObject("stats");
    paveStat->SetOptStat(0);
    cSysError.Modified();
    cSysError.Print(Form("plots/%s.pdf",fSystematicUncertaintyTotal.GetName()));
    cSysError.Write();
    
  }
  outFile->Close();
  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;
}
