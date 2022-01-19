// SystematicsCorrelations.cpp

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

void SystematicsCorrelations(const int points = kNPoints, const bool cutVar = true, const bool binCountingVar = true, const bool expVar = true, const bool sigmoidVar = true, const char *outFileName = "SystematicsAllEPt")
{
  gStyle->SetOptStat(110001110);
  gStyle->SetStatX(0.87);
  gStyle->SetStatY(0.85);
  TStopwatch swatch;
  swatch.Start(true);

  TFile *specFile = TFile::Open(Form("%s/SpectraProtonSys.root", kOutDir));
  TFile *effFile = TFile::Open(Form("%s/EfficiencyProtonSys.root", kOutDir));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iC = 0; iC < kNCentClasses; ++iC)
  {
    TH2D fSystematicsCorrelationsDCAxy(Form("fSystematicsCorrelationsDCAxy_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),5,0.045,0.095,1000,-.5,.5);
    TH2D fSystematicsCorrelationsDCAz(Form("fSystematicsCorrelationsDCAz_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),31,0.475,2.025,1000,-.5,.5);
    TH2D fSystematicsCorrelationsTPCCls(Form("fSystematicsCorrelationsTPCCls_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),12,59,83,1000,-.5,.5);
    TH2D fSystematicsCorrelationsPID(Form("fSystematicsCorrelationsPID_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),Form("%.0f-%.0f%%", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]),16,3.24,3.56,1000,-.5,.5);
    //TH2D fSystematicsCorrelationsROI();
    fSystematicsCorrelationsDCAxy.Write();
    fSystematicsCorrelationsDCAz.Write();
    fSystematicsCorrelationsTPCCls.Write();
    fSystematicsCorrelationsPID.Write();

    // track cuts
    int cutVariable=0;
    int bkgFlag = 1;
    int sigmoidFlag = 1;
    int iNsigma = 1;
    TH1D *hDefault = (TH1D *)specFile->Get(Form("_1_1_1/fRatio_%.0f_%.0f", kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
    for (int iTrackCuts=1; iTrackCuts<kNTrackCuts; ++iTrackCuts){
      TString tmpCutSettings = trackCutSettings[iTrackCuts];
      auto cutIndex = trackCutIndexes[iTrackCuts];
      auto tmpCutIndex = Form("%d",cutIndex);
      auto fullCutSettingsSigm = Form("%s%s_%d_%d_%d", tmpCutSettings.Data(), tmpCutIndex, bkgFlag, sigmoidFlag, iNsigma);
      //std::cout << fullCutSettingsSigm <<std::endl;
      TH1D *h = (TH1D *)specFile->Get(Form("%s/fRatio_%.0f_%.0f", fullCutSettingsSigm, kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]));
      for(int iPtBins=4;iPtBins<24;++iPtBins){
        if (tmpCutSettings.EqualTo("dcaz"))
        {
          int cutVal = cutIndex;
          if (cutIndex > 1) cutVal = cutIndex+1;
          fSystematicsCorrelationsDCAz.Fill(kCutDCAz[cutVal],h->GetBinContent(iPtBins+1)-hDefault->GetBinContent(iPtBins+1));
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
    double correlationDCAxy = fSystematicsCorrelationsDCAxy.GetCorrelationFactor();
    std::cout << "CorrelationDCAxy = " << correlationDCAxy <<std::endl;
    double correlationTPCCls = fSystematicsCorrelationsTPCCls.GetCorrelationFactor();
    std::cout << "CorrelationTPCCls = " << correlationTPCCls <<std::endl;
    double correlationPID = fSystematicsCorrelationsPID.GetCorrelationFactor();
    std::cout << "CorrelationPID = " << correlationPID <<std::endl;
    std::cout << "* * * * * * * * * *"<<std::endl;

    fSystematicsCorrelationsDCAz.Write();
    fSystematicsCorrelationsDCAxy.Write();
    fSystematicsCorrelationsTPCCls.Write();
    fSystematicsCorrelationsPID.Write();

  }
  outFile->Close();
  swatch.Stop();
  std::cout << "Elapsed (real) time = " << swatch.RealTime() << " s" << std::endl;
}
