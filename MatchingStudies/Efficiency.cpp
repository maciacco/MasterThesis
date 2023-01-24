#include "Config.h"

void Efficiency (const char* inFileName = "out_q.root", const char* outFileName = "MatchingTPCTOF_q.root")
{
  TFile outFile(Form("%s/%s", kOutDir, outFileName), "recreate");
  TFile inFile(Form("%s/%s", kOutDir, inFileName));
  TH1D hEffRatioVsCent("hEffRatioVsCent", ";Centrality (%);#epsilon_{Matching}^{-}/#epsilon_{Matching}^{+}", kNCentBins, kCentBinLimits);
  for (int iCent{1}; iCent < kNCentBins+1; ++iCent)
  {
    TH1D* hEff[2];
    for (int iC{0}; iC < 2; ++iC)
    {
      TH2D* hTPC = (TH2D*)inFile.Get(Form("f%sTPC", kChargeFlag[iC][0]));
      TH2D* hTPCTOF = (TH2D*)inFile.Get(Form("f%sTPCTOF", kChargeFlag[iC][0]));
      TH1D* hTPCproj = (TH1D*)hTPC->ProjectionY(Form("f%sTPCproj", kChargeFlag[iC][0]), iCent, iCent);
      TH1D* hTPCTOFproj = (TH1D*)hTPCTOF->ProjectionY(Form("f%sTPCTOFproj", kChargeFlag[iC][0]), iCent, iCent);
      hEff[iC] = new TH1D(*hTPCTOFproj);
      hEff[iC]->SetName(Form("f%sEff", kChargeFlag[iC][0]));
      hEff[iC]->GetYaxis()->SetTitle("#epsilon_{Matching}");
      hEff[iC]->Divide(hTPCTOFproj, hTPCproj, 1, 1, "B");
      outFile.mkdir(Form("%.0f_%.0f", kCentBinLimits[iCent-1], kCentBinLimits[iCent]));
      outFile.cd(Form("%.0f_%.0f", kCentBinLimits[iCent-1], kCentBinLimits[iCent]));
      hTPCproj->Write();
      hTPCTOFproj->Write();
      hEff[iC]->Write();
    }
    TH1D* hEffRatio = new TH1D(*hEff[0]);
    hEffRatio->SetName("fEffRatio");
    hEffRatio->GetYaxis()->SetTitle("#epsilon_{Matching}^{-}/#epsilon_{Matching}^{+}");
    hEffRatio->Divide(hEff[0], hEff[1], 1, 1);
    hEffRatio->GetXaxis()->SetRangeUser(0.7, 1.6);
    hEffRatio->Fit("pol0");
    hEffRatioVsCent.SetBinContent(iCent, hEffRatio->GetFunction("pol0")->GetParameter(0));
    hEffRatioVsCent.SetBinError(iCent, hEffRatio->GetFunction("pol0")->GetParError(0));
    hEffRatio->Write();
  }
  outFile.cd();
  hEffRatioVsCent.Write();
  outFile.Close();
}