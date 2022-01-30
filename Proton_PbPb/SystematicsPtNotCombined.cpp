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
const int used_pt_bins = 20;
const int nTrials=10000;

void SystematicsPtNotCombined(const int points = kNPoints, const bool cutVar = true, const bool binCountingVar = true, const bool expVar = true, const bool sigmoidVar = true, const char *outFileName = "SystematicsAllEPtNotCombined")
{
  gStyle->SetTextFont(44);
  gStyle->SetOptFit(0);
  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.85);
  TStopwatch swatch;
  swatch.Start(true);

  TFile *specFile = TFile::Open(Form("%s/SpectraProtonSys.root", kOutDir));
  //TFile *fileEffFit=TFile::Open("out/SpectraProton_MC21l5_raw_fitEff.root");
  TFile *fG4 = TFile::Open("out/SpectraProton_MC21l5_raw_primaryInjected.root");      
  TFile *inFileSec = TFile::Open(Form("%s/PrimaryProton.root", kOutDir));
  TFile *effFile = TFile::Open(Form("%s/EfficiencyProtonSys.root", kOutDir));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iC = 0; iC < kNCentClasses; ++iC) // TODO: extend the analysis to the third centrality class as well
  {
    gStyle->SetOptStat(1);
    //TH1D fSystematicUncertaintyDCAxy(Form("fSystematicUncertaintyDCAxy_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyDCAz(Form("fSystematicUncertaintyDCAz_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTPCCls(Form("fSystematicUncertaintyTPCCls_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyPID(Form("fSystematicUncertaintyPID_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyROI(Form("fSystematicUncertaintyROI_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyEff(Form("fSystematicUncertaintyEff_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTFF(Form("fSystematicUncertaintyTFF_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyPrim(Form("fSystematicUncertaintyPrim_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyPrimG3G4(Form("fSystematicUncertaintyPrimG3G4_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    //TH1D fSystematicUncertaintyEffFit(Form("fSystematicUncertaintyEffFit_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTotal(Form("fSystematicUncertaintyTotal_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fSystematicUncertaintyTotalPtCorrelated(Form("fSystematicUncertaintyTotalPtCorrelated_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    TH1D fRatioFromVariationsTot(Form("fRatioFromVariationsTot_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), kNPtBins, kPtBins);
    double ratioBins[1000];
    for (int iRatioBins = 0; iRatioBins < 3000; iRatioBins++){
      ratioBins[iRatioBins]=0.85+iRatioBins*0.0001;
    }
    //TH2D fRatiosVsPtDCAxy(Form("fRatiosVsPtDCAxy_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtDCAz(Form("fRatiosVsPtDCAz_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtTPCCls(Form("fRatiosVsPtTPCCls_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtPID(Form("fRatiosVsPtPID_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtROI(Form("fRatiosVsPtROI_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtPrim(Form("fRatiosVsPtPrim_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtPrimG3G4(Form("fRatiosVsPtPrimG3G4_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    //TH2D fRatiosVsPtEffFit(Form("fRatiosVsPtEffFit_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    TH2D fRatiosVsPtTot(Form("fRatiosVsPt_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins,2999,ratioBins);
    
    // total cuts
    int cutVariable=0;
    int bkgFlag = 1;
    int sigmoidFlag = 1;
    int iNsigma = 1;
    for (int iTrackCuts=0; iTrackCuts<kNTrackCuts; ++iTrackCuts){
      for (int iROI=0; iROI<3; ++iROI){
        for (int iG3G4Prim=1; iG3G4Prim<2; ++iG3G4Prim){
          for (int iSigmoid=0; iSigmoid<2; ++iSigmoid){
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
            for(int iPtBins=4;iPtBins<24;++iPtBins){
              fRatiosVsPtTot.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
            }
          }
        }
      }
    }

    // DCAz cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 1;
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
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_1", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtDCAz.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    // PID cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 1;
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
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_1", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtPID.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    // TPCCls cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 1;
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
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_1", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtTPCCls.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    // DCAxy cuts
    /* cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 1;
    iNsigma = 1;
    for (int iTrackCuts=11; iTrackCuts<kNTrackCuts; ++iTrackCuts){
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
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtDCAxy.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    } */

    TH1D *hDefault = (TH1D *)specFile->Get(Form("_1_1_1_1/fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));

    // efficiency error
    TH1D *h_a_eff = (TH1D *)effFile->Get(Form("_/fAEff_TOF_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    TH1D *h_m_eff = (TH1D *)effFile->Get(Form("_/fMEff_TOF_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    for(int iPtBin=5;iPtBin<25;++iPtBin){
      double err_eff_a=h_a_eff->GetBinError(iPtBin);
      double eff_a=h_a_eff->GetBinContent(iPtBin);
      double err_eff_m=h_m_eff->GetBinError(iPtBin);
      double eff_m=h_m_eff->GetBinContent(iPtBin);
      fSystematicUncertaintyEff.SetBinContent(iPtBin,TMath::Sqrt(err_eff_a*err_eff_a/eff_a/eff_a+err_eff_m*err_eff_m/eff_m/eff_m));
      fSystematicUncertaintyEff.SetBinError(iPtBin,0);
    }

    // primary fraction (TFF) error
    for(int iPtBin=5;iPtBin<25;++iPtBin){
      double primaryRelativeError[2];
      for (int iMatt = 0; iMatt < 2; ++iMatt){
        TH1D *h_sec = (TH1D*)inFileSec->Get(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        double primaryError=h_sec->GetBinError(iPtBin);
        double primary=h_sec->GetBinContent(iPtBin);
        /* TF1 *sec_f = (TF1 *)inFileSec->Get(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        TH2D *sec_f_cov = (TH2D *)inFileSec->Get(Form("f%sCovMat_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
        TH1D h_tmp("h_tmp","h_tmp",kNPtBins,kPtBins);
        double primary = sec_f->Eval(h_tmp.GetXaxis()->GetBinCenter(iPtBin));
        double pt_center = h_tmp.GetXaxis()->GetBinCenter(iPtBin);
        double par_0 = sec_f->GetParameter(0);
        double par_1 = sec_f->GetParameter(1);
        double var_par_0 = sec_f_cov->GetBinContent(1,1);
        double var_par_1 = sec_f_cov->GetBinContent(2,2);
        double cov = sec_f_cov->GetBinContent(1,2);
        double exponential = TMath::Exp(par_1*pt_center);
        double denominator = 1+par_0*exponential;
        double first_derivative_par_0 = -exponential/denominator/denominator;
        double first_derivative_par_1 = -par_0*exponential*pt_center/denominator/denominator;
        double primaryError = TMath::Sqrt(first_derivative_par_0*first_derivative_par_0*var_par_0+first_derivative_par_1*first_derivative_par_1*var_par_1+2*first_derivative_par_0*first_derivative_par_1*cov);
        */primaryRelativeError[iMatt]=primaryError/primary;
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
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_1", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtROI.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    // prim error
    cutVariable=0;
    iNsigma = 1;
    for (int iPrim=0; iPrim<2; ++iPrim){
      auto tmpCutSettings = Form("");
      auto tmpCutIndex = Form("");
      sigmoidFlag=iPrim;
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_1", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtPrim.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    // prim error - G3 vs. G4
    cutVariable=0;
    iNsigma = 1;
    for (int iPrim=0; iPrim<2; ++iPrim){
      auto tmpCutSettings = Form("");
      auto tmpCutIndex = Form("");
      sigmoidFlag=iPrim;
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_%d", tmpCutSettings, tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma, iPrim);
      std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fRatiosVsPtPrimG3G4.Fill(kPtBins[iPtBins],h->GetBinContent(iPtBins+1));
      }
    }

    for(int iPtBins=5;iPtBins<25;++iPtBins){
      // track cuts
      TH1D *proj=fRatiosVsPtDCAz.ProjectionY("py",iPtBins,iPtBins);
      remove_outliers(proj,1.55);
      fSystematicUncertaintyDCAz.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyDCAz.SetBinError(iPtBins,0);

      proj=fRatiosVsPtPID.ProjectionY("py",iPtBins,iPtBins);
      //remove_outliers(proj,2.);
      fSystematicUncertaintyPID.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyPID.SetBinError(iPtBins,0);

      proj=fRatiosVsPtTPCCls.ProjectionY("py",iPtBins,iPtBins);
      remove_outliers(proj,1.55);
      fSystematicUncertaintyTPCCls.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyTPCCls.SetBinError(iPtBins,0);

      /* proj=fRatiosVsPtDCAxy.ProjectionY("py",iPtBins,iPtBins);
      fSystematicUncertaintyDCAxy.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyDCAxy.SetBinError(iPtBins,0); */

      // roi
      proj=fRatiosVsPtROI.ProjectionY("py",iPtBins,iPtBins);
      remove_outliers(proj,1.4);
      fSystematicUncertaintyROI.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyROI.SetBinError(iPtBins,0);

      // prim
      proj=fRatiosVsPtPrim.ProjectionY("py",iPtBins,iPtBins);
      //remove_outliers(proj,2.);
      fSystematicUncertaintyPrim.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyPrim.SetBinError(iPtBins,0);
      
      // prim G3 vs. G4
      proj=fRatiosVsPtPrimG3G4.ProjectionY("py",iPtBins,iPtBins);
      //remove_outliers(proj,2.);
      fSystematicUncertaintyPrimG3G4.SetBinContent(iPtBins,proj->GetRMS()/proj->GetMean());
      fSystematicUncertaintyPrimG3G4.SetBinError(iPtBins,0);

      // tot
      proj=fRatiosVsPtTot.ProjectionY("py",iPtBins,iPtBins);
      remove_outliers(proj);
      
      fRatioFromVariationsTot.SetBinContent(iPtBins,proj->GetMean());
      fRatioFromVariationsTot.SetBinError(iPtBins,0);
    }
    fRatioFromVariationsTot.Write();
    //fRatiosVsPtDCAxy.Write();
    fRatiosVsPtDCAz.Write();
    fRatiosVsPtPID.Write();
    fRatiosVsPtTPCCls.Write();
    fRatiosVsPtPrim.Write();
    fRatiosVsPtPrimG3G4.Write();
    //fSystematicUncertaintyDCAxy.Write();
    fSystematicUncertaintyDCAz.Write();
    fSystematicUncertaintyPID.Write();
    fSystematicUncertaintyTPCCls.Write();
    fSystematicUncertaintyPrim.Write();
    fSystematicUncertaintyROI.Write();
    fSystematicUncertaintyEff.Write();
    fSystematicUncertaintyTFF.Write();
    fSystematicUncertaintyPrimG3G4.Write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * //
    //   C O M P U T E   P T   C O R R E L A T I O N S   //
    // * * * * * * * * * * * * * * * * * * * * * * * * * //
    //TH2D fSystematicsCorrelationsDCAxy(Form("fSystematicsCorrelationsDCAxy_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),5,0.045,0.095,1000,-.5,.5);
    TH2D fSystematicsCorrelationsDCAz(Form("fSystematicsCorrelationsDCAz_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),31,0.475,2.025,1000,-.5,.5);
    TH2D fSystematicsCorrelationsTPCCls(Form("fSystematicsCorrelationsTPCCls_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),24,58.5,82.5,1000,-.5,.5);
    TH2D fSystematicsCorrelationsPID(Form("fSystematicsCorrelationsPID_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),25,3.245,3.505,1000,-.5,.5);
    TH2D fSystematicsCorrelationsROI(Form("fSystematicsCorrelationsROI_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),51,7.49,8.51,1000,-.5,.5);

    // roi
    cutVariable = 0;
    bkgFlag = 1;
    sigmoidFlag = 1;
    iNsigma = 1;
    hDefault = (TH1D *)specFile->Get(Form("_1_1_1_1/fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    for (int iROI=0;iROI<3;++iROI){
      if (iROI==1)continue;
      TString tmpCutSettings = Form("");
      auto cutIndex = 0;
      auto tmpCutIndex = Form("",cutIndex);
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_1", tmpCutSettings.Data(), tmpCutIndex, bkgFlag, sigmoidFlag, iROI);
      //std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        fSystematicsCorrelationsROI.Fill(kROIVariation[iROI],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
      }
    }
        
    // track cuts
    cutVariable=0;
    bkgFlag = 1;
    sigmoidFlag = 1;
    iNsigma = 1;
    hDefault = (TH1D *)specFile->Get(Form("_1_1_1_1/fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    for (int iTrackCuts=1; iTrackCuts<kNTrackCuts; ++iTrackCuts){
      TString tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d_1", tmpCutSettings.Data(), tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      //std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        if (tmpCutSettings.EqualTo("dcaz"))
        {
          int cutVal = cutIndex;
          if (cutIndex > 1) cutVal = cutIndex+1;
          fSystematicsCorrelationsDCAz.Fill(kCutDCAz[cutVal],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
        }
        /* if (tmpCutSettings.EqualTo("dcaxy"))
        {
          int cutVal = cutIndex;
          if (cutIndex > 1) cutVal = cutIndex+1;
          fSystematicsCorrelationsDCAxy.Fill(kCutDCAxy[cutVal],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
        } */
        if (tmpCutSettings.EqualTo("tpc"))
        {
          int cutVal = cutIndex;
          if (cutIndex > 1) cutVal = cutIndex+1;
          fSystematicsCorrelationsTPCCls.Fill(kCutTPCClusters[cutVal],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
        }
        if (tmpCutSettings.EqualTo("pid"))
        {
          int cutVal = cutIndex+1;
          fSystematicsCorrelationsPID.Fill(kTPCPidSigmas[cutVal],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
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
    /* double correlationDCAxy = fSystematicsCorrelationsDCAxy.GetCorrelationFactor();
    std::cout << "CorrelationDCAxy = " << correlationDCAxy <<std::endl;
     */double correlationTPCCls = fSystematicsCorrelationsTPCCls.GetCorrelationFactor();
    std::cout << "CorrelationTPCCls = " << correlationTPCCls <<std::endl;
    double correlationPID = fSystematicsCorrelationsPID.GetCorrelationFactor();
    std::cout << "CorrelationPID = " << correlationPID <<std::endl;
    std::cout << "* * * * * * * * * *"<<std::endl;
    double correlationROI = fSystematicsCorrelationsROI.GetCorrelationFactor();
    std::cout << "CorrelationROI = " << correlationROI <<std::endl;
    std::cout << "* * * * * * * * * *"<<std::endl;

    fSystematicsCorrelationsDCAz.Write();
    //fSystematicsCorrelationsDCAxy.Write();
    fSystematicsCorrelationsTPCCls.Write();
    fSystematicsCorrelationsPID.Write();
    fSystematicsCorrelationsROI.Write();
    
    for(int iPtBins=5;iPtBins<25;++iPtBins){
      // sigmas
      double sigma_DCAz=fSystematicUncertaintyDCAz.GetBinContent(iPtBins);
      double sigma_PID=fSystematicUncertaintyPID.GetBinContent(iPtBins);
      double sigma_TPCCls=fSystematicUncertaintyTPCCls.GetBinContent(iPtBins);
      // double sigma_DCAxy=fSystematicUncertaintyDCAxy.GetBinContent(iPtBins);
      double sigma_ROI=fSystematicUncertaintyROI.GetBinContent(iPtBins);
      double sigma_Prim=fSystematicUncertaintyPrim.GetBinContent(iPtBins);
      double sigma_Prim_G3G4=fSystematicUncertaintyPrimG3G4.GetBinContent(iPtBins);
      double sigma_Eff=fSystematicUncertaintyEff.GetBinContent(iPtBins);
      double sigma_TFF=fSystematicUncertaintyTFF.GetBinContent(iPtBins);
      //double sigma_EffFit=fSystematicUncertaintyEffFit.GetBinContent(iPtBins);

      // fractions (uncorrelated)
      double fraction_uncorr_DCAz=1.-std::abs(correlationDCAz);
      double fraction_uncorr_PID=1.-std::abs(correlationPID);
      double fraction_uncorr_TPCCls=1.-std::abs(correlationTPCCls);
      double fraction_uncorr_ROI=1.-std::abs(correlationROI);
      // double fraction_uncorr_DCAxy=1.-std::abs(correlationDCAxy);
      double totSys=0;
      double totSysPtCorrelated=0;
      if (iC<2 || (iC==2 && iPtBins!=18)){
        totSys=TMath::Sqrt(sigma_DCAz*sigma_DCAz*fraction_uncorr_DCAz*fraction_uncorr_DCAz
        +sigma_PID*sigma_PID*fraction_uncorr_PID*fraction_uncorr_PID
        +sigma_TPCCls*sigma_TPCCls*fraction_uncorr_TPCCls*fraction_uncorr_TPCCls
        //+sigma_DCAxy*sigma_DCAxy*fraction_uncorr_DCAxy*fraction_uncorr_DCAxy
        +sigma_ROI*sigma_ROI*fraction_uncorr_ROI*fraction_uncorr_ROI
        +sigma_Prim*sigma_Prim
        +sigma_Eff*sigma_Eff
        +sigma_TFF*sigma_TFF
        /* +sigma_EffFit*sigma_EffFit */);
        fSystematicUncertaintyTotal.SetBinContent(iPtBins,totSys);
        fSystematicUncertaintyTotal.SetBinError(iPtBins,0);

        totSysPtCorrelated=TMath::Sqrt(sigma_DCAz*sigma_DCAz*correlationDCAz*correlationDCAz
        +sigma_PID*sigma_PID*correlationPID*correlationPID
        +sigma_TPCCls*sigma_TPCCls*correlationTPCCls*correlationTPCCls
        +sigma_ROI*sigma_ROI*correlationROI*correlationROI
        //+sigma_Prim_G3G4*sigma_Prim_G3G4
        //+sigma_DCAxy*sigma_DCAxy*correlationDCAxy*correlationDCAxy
        );
      }

      else {
        totSys=TMath::Sqrt(sigma_DCAz*sigma_DCAz*fraction_uncorr_DCAz*fraction_uncorr_DCAz
        //+sigma_PID*sigma_PID*fraction_uncorr_PID*fraction_uncorr_PID
        +sigma_TPCCls*sigma_TPCCls*fraction_uncorr_TPCCls*fraction_uncorr_TPCCls
        //+sigma_DCAxy*sigma_DCAxy*fraction_uncorr_DCAxy*fraction_uncorr_DCAxy
        +sigma_ROI*sigma_ROI*fraction_uncorr_ROI*fraction_uncorr_ROI
        +sigma_Prim*sigma_Prim
        +sigma_Prim_G3G4*sigma_Prim_G3G4
        +sigma_Eff*sigma_Eff
        +sigma_TFF*sigma_TFF
        /* +sigma_EffFit*sigma_EffFit */);
        fSystematicUncertaintyTotal.SetBinContent(iPtBins,totSys);
        fSystematicUncertaintyTotal.SetBinError(iPtBins,0);

        totSysPtCorrelated=TMath::Sqrt(sigma_DCAz*sigma_DCAz*correlationDCAz*correlationDCAz
        //+sigma_PID*sigma_PID*correlationPID*correlationPID
        +sigma_TPCCls*sigma_TPCCls*correlationTPCCls*correlationTPCCls
        +sigma_ROI*sigma_ROI*correlationROI*correlationROI
        //+sigma_Prim_G3G4*sigma_Prim_G3G4
        //+sigma_DCAxy*sigma_DCAxy*correlationDCAxy*correlationDCAxy
        );
      }
      fSystematicUncertaintyTotalPtCorrelated.SetBinContent(iPtBins,totSysPtCorrelated);
      fSystematicUncertaintyTotalPtCorrelated.SetBinError(iPtBins,0);
    }
    fSystematicUncertaintyTotal.Write();
    fSystematicUncertaintyTotalPtCorrelated.Write();

    // Compute pt correlated systematic uncertainty
    TH1D hRatio(Form("fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),kNPtBins,kPtBins);
    for (int iPtBin=5;iPtBin<kNPtBins;++iPtBin){
      hRatio.SetBinContent(iPtBin,fRatioFromVariationsTot.GetBinContent(iPtBin));
      hRatio.SetBinError(iPtBin,hRatio.GetBinContent(iPtBin)*fSystematicUncertaintyTotal.GetBinContent(iPtBin));
    }
    hRatio.Fit("pol0");
    hRatio.Write();

    TH1D h_trial("h_trial","h_trial",kNPtBins,kPtBins);
    TH1D fRatioDistributionTrials(Form("fRatioDistributionTrials_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),500,0.95,1.05);
    for (int iTrial=0;iTrial<nTrials;++iTrial){
      double nsigma=gRandom->Gaus(0,1);
      for (int iPtBin=5;iPtBin<kNPtBins;++iPtBin){
        h_trial.SetBinContent(iPtBin,fRatioFromVariationsTot.GetBinContent(iPtBin)+nsigma*fSystematicUncertaintyTotalPtCorrelated.GetBinContent(iPtBin)*fRatioFromVariationsTot.GetBinContent(iPtBin));
        h_trial.SetBinError(iPtBin,hRatio.GetBinContent(iPtBin)*fSystematicUncertaintyTotal.GetBinContent(iPtBin));
      }
      h_trial.Fit("pol0","Q");
      if( (h_trial.GetFunction("pol0")->GetChisquare()/h_trial.GetFunction("pol0")->GetNDF())<2.)
      fRatioDistributionTrials.Fill(h_trial.GetFunction("pol0")->GetParameter(0));
    }

    TCanvas cRatio(Form("cRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]), "cRatio");
    cRatio.SetTicks(1, 1);
    hRatio.SetMarkerStyle(20);
    hRatio.SetMarkerSize(0.8);
    hRatio.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatio.GetYaxis()->SetTitle("Systematic Uncertainty");
    hRatio.GetXaxis()->SetRangeUser(1.0,2.0);
    hRatio.GetYaxis()->SetRangeUser(0.88, 1.12);
    hRatio.SetStats(0);
    hRatio.Draw("");
    TLatex chi2(1.5, 1.08, Form("#chi^{2}/NDF = %.2f/%d", hRatio.GetFunction("pol0")->GetChisquare(), hRatio.GetFunction("pol0")->GetNDF()));
    chi2.SetTextSize(28);
    TLatex p0(1.5, 1.10, Form("R = %.4f #pm %.4f", hRatio.GetFunction("pol0")->GetParameter(0), hRatio.GetFunction("pol0")->GetParError(0)));
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
    fSystematicUncertaintyTotal.GetXaxis()->SetRangeUser(1.,2.);
    fSystematicUncertaintyTotal.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fSystematicUncertaintyTotal.GetYaxis()->SetTitle("Systematic Uncertainty");
    fSystematicUncertaintyTotal.SetMinimum(0.);
    TCanvas cSysError(fSystematicUncertaintyTotal.GetName(),fSystematicUncertaintyTotal.GetTitle());
    fSystematicUncertaintyTotal.Draw("histo");
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
