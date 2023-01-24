#include "Config.h"

constexpr const char* periods[] = {"q", "r"};
constexpr const char* dataMc[] = {"data", "mc"};
constexpr const char* drawOpt[] = {"", "same"};
constexpr Color_t colors[] = {kBlue, kRed};
constexpr double legPos[] = {0.1761038961038961, 0.2368003465285372, 0.4763636363636364, 0.4373387045366125};
constexpr int markerStyle[2] = {20, 21};

void CombineResults ()
{
  TFile outFile(Form("%s/CombineResults.root", kOutDirCombine), "recreate");

  TCanvas cEffRatioVsPeriod[2];
  TLegend *lEffRatioVsPeriod[2];
  for (int i{0}; i < 2; ++i)
  {
    cEffRatioVsPeriod[i].SetName(Form("QvsRmatching_%s", dataMc[i]));
    lEffRatioVsPeriod[i] = new TLegend(legPos[0], legPos[1], legPos[2], legPos[3]);
  }

  // combine plots of average efficiency ratios wrt period or production
  TH1D* hEffRatio[2][2];
  TFile* inFile[2][2];
  for (int iP{0}; iP < 2; ++iP)
  {
    TCanvas cEffRatioDataVsMC(Form("DataVsMcMatching_%s", periods[iP]), Form("DataVsMcMatching_%s", periods[iP]));
    TLegend lEffRatioDataVsMC(legPos[0], legPos[1], legPos[2], legPos[3]);
    for (int iD{0}; iD < 2; ++iD)
    {
      inFile[iP][iD] = new TFile(Form("results_%s_int/MatchingTPCTOF_%s.root", dataMc[iD], periods[iP]));
      hEffRatio[iP][iD] = (TH1D*)inFile[iP][iD]->Get("hEffRatioVsCent");
      hEffRatio[iP][iD]->GetYaxis()->SetRangeUser(0.94, 1.02);
      hEffRatio[iP][iD]->SetMarkerStyle(markerStyle[iP]);
      hEffRatio[iP][iD]->SetMarkerSize(1.5);
      hEffRatio[iP][iD]->SetMarkerColor(colors[iD]);
      hEffRatio[iP][iD]->SetLineColor(colors[iD]);
      cEffRatioDataVsMC.cd();
      hEffRatio[iP][iD]->Draw(drawOpt[iD]);
      lEffRatioDataVsMC.AddEntry(hEffRatio[iP][iD], Form("%s_%s", dataMc[iD], periods[iP]));
      cEffRatioVsPeriod[iD].cd();
      lEffRatioVsPeriod[iD]->AddEntry(hEffRatio[iP][iD], Form("%s_%s", dataMc[iD], periods[iP]));
      hEffRatio[iP][iD]->Draw(drawOpt[iP]);
    }
    cEffRatioDataVsMC.cd();
    lEffRatioDataVsMC.Draw("same");
    outFile.cd();
    cEffRatioDataVsMC.Write();

    // data/mc ratio
    TH1D* hDataMcRatio = new TH1D(*hEffRatio[iP][0]);
    hDataMcRatio->SetName(Form("hDataMcRatio_%s", periods[iP]));
    hDataMcRatio->GetYaxis()->SetTitle("#frac{(#epsilon^{-}_{match}/#epsilon^{+}_{match})_{Data}}{(#epsilon^{-}_{match}/#epsilon^{+}_{match})_{MC}}");
    hDataMcRatio->Divide(hEffRatio[iP][1]);
    hDataMcRatio->Write();
  }
  // write to file
  for (int i{0}; i < 2; ++i)
  {
    outFile.cd();
    cEffRatioVsPeriod[i].cd();
    lEffRatioVsPeriod[i]->Draw("same");
    cEffRatioVsPeriod[i].Write();
  }

  // combine efficiency plots
  TH1D* hEff[kNCentBins][2][2][2];
  for (int iCent{0}; iCent < kNCentBins; ++iCent)
  {
    // get all plots
    outFile.mkdir(Form("%.0f_%.0f", kCentBinLimits[iCent], kCentBinLimits[iCent+1]));
    for (int iP{0}; iP < 2; ++iP)
    {
      for (int iD{0}; iD < 2; ++iD)
      {
        for (int iC{0}; iC < 2; ++iC)
        {
          hEff[iCent][iP][iD][iC] = (TH1D*)inFile[iP][iD]->Get(Form("%.0f_%.0f/f%sEff", kCentBinLimits[iCent], kCentBinLimits[iCent+1], kChargeFlag[iC][0]));
        }
      }
    }
    // plot pos and neg + ratio
    for (int iP{0}; iP < 2; ++iP)
    {
      for (int iD{0}; iD < 2; ++iD)
      {
        TCanvas cPosNeg(Form("cPosNeg_%s_%s", dataMc[iD], periods[iP]), Form("cPosNeg_%s_%s", dataMc[iD], periods[iP]));
        TLegend l(0.5729870129870129, 0.2227611427535554, 0.8732467532467533, 0.42276114275355536);
        for (int iC{0}; iC < 2; ++iC)
        {
          hEff[iCent][iP][iD][iC]->SetMarkerStyle(0);
          hEff[iCent][iP][iD][iC]->SetMarkerSize(0);
          hEff[iCent][iP][iD][iC]->SetLineWidth(2);
          hEff[iCent][iP][iD][iC]->SetLineColor(colors[iC]);
          hEff[iCent][iP][iD][iC]->GetYaxis()->SetRangeUser(0., 1.);
          hEff[iCent][iP][iD][iC]->Draw(drawOpt[iC]);
          l.AddEntry(hEff[iCent][iP][iD][iC], Form("LHC18%s, %s, %s", periods[iP], dataMc[iD], kChargeFlag[iC][1]));
        }
        outFile.cd(Form("%.0f_%.0f", kCentBinLimits[iCent], kCentBinLimits[iCent+1]));
        l.Draw("same");
        cPosNeg.Write();

        // ratio
        TH1D* r = new TH1D(*hEff[iCent][iP][iD][0]);
        r->SetName(Form("fEffRatio_%s_%s", dataMc[iD], periods[iP]));
        r->Divide(hEff[iCent][iP][iD][0], hEff[iCent][iP][iD][1]);
        r->GetYaxis()->SetTitle("#epsilon^{-}_{Matching} / #epsilon^{+}_{Matching}");
        r->GetYaxis()->SetRangeUser(0.9, 1.1);
        r->Write();
      }
    }
  }

  outFile.Close();
}

void CombineResults (const char *period = "qr")
{
  TFile outFile(Form("%s/CombineResults.root", kOutDirCombine), "recreate");
  // combine plots of average efficiency ratios wrt period or production
  TH1D* hEffRatio[2];
  TFile* inFile[2];
  TCanvas cEffRatioDataVsMC("DataVsMcMatching_qr", "DataVsMcMatching_qr");
  TLegend lEffRatioDataVsMC(legPos[0], legPos[1], legPos[2], legPos[3]);
  for (int iD{0}; iD < 2; ++iD)
  {
    inFile[iD] = new TFile(Form("results_%s_int_merge/MatchingTPCTOF_qr.root", dataMc[iD]));
    hEffRatio[iD] = (TH1D*)inFile[iD]->Get("hEffRatioVsCent");
    hEffRatio[iD]->GetYaxis()->SetRangeUser(0.94, 1.02);
    hEffRatio[iD]->SetMarkerStyle(markerStyle[0]);
    hEffRatio[iD]->SetMarkerSize(1.5);
    hEffRatio[iD]->SetMarkerColor(colors[iD]);
    hEffRatio[iD]->SetLineColor(colors[iD]);
    cEffRatioDataVsMC.cd();
    hEffRatio[iD]->Draw(drawOpt[iD]);
    lEffRatioDataVsMC.AddEntry(hEffRatio[iD], Form("%s_qr", dataMc[iD]));
  }
  cEffRatioDataVsMC.cd();
  lEffRatioDataVsMC.Draw("same");
  outFile.cd();
  cEffRatioDataVsMC.Write();

  // data/mc ratio
  TH1D* hDataMcRatio = new TH1D(*hEffRatio[0]);
  hDataMcRatio->SetName("hDataMcRatio_qr");
  hDataMcRatio->GetYaxis()->SetTitle("#frac{(#epsilon^{-}_{match}/#epsilon^{+}_{match})_{Data}}{(#epsilon^{-}_{match}/#epsilon^{+}_{match})_{MC}}");
  hDataMcRatio->Divide(hEffRatio[1]);
  hDataMcRatio->Write();

  // systematics
  TH1D* hSystematics = new TH1D(*hDataMcRatio);
  hSystematics->SetName("hSystematics");
  hSystematics->GetYaxis()->SetTitle("Systematic uncertainty");
  hSystematics->GetYaxis()->SetRangeUser(0., 0.02);
  for (int iB{1}; iB < hSystematics->GetNbinsX() + 1; ++iB)
  {
    hSystematics->SetBinContent(iB, 1 - hSystematics->GetBinContent(iB));
    hSystematics->SetBinError(iB, hSystematics->GetBinError(iB));
  }
  hSystematics->Write();
}