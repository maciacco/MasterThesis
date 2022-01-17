#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include "../utils/Config.h"

using namespace proton;
const int kMarkerSize[] = {20,24};
const char* kFittingMethod[] = {"TFF","RooFit"};
const Color_t kHistColor[] = {kRed,kBlue};

void PlotSecondaryCorrection(const char *inFileRooFitName = "PrimaryProton_ROO", const char *inFileTFFName = "PrimaryProton", const char *outFileCompareName = "PrimaryCompareRooFitTFF", const char *inFileSpectraTFFName = "SpectraProton_MC21l5_raw", const char* inFileSpectraRooName = "SpectraProton_MC21l5_raw_ROO"){
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  // compare RooFit and ROOT::TFractionFitter
  TFile *InFile[2];
  InFile[0] = new TFile(Form("out/%s.root",inFileTFFName));
  InFile[1] = new TFile(Form("out/%s.root",inFileRooFitName));
  TFile inFileSpectraTFF(Form("out/%s.root",inFileSpectraTFFName));
  TFile inFileSpectraRoo(Form("out/%s.root",inFileSpectraRooName));
  TFile OutFileCompare(Form("%s.root",outFileCompareName),"recreate");
  for (int iC = 0; iC < kNCentClasses; ++iC){
    TH1D *hPrimary[2][2];
    for (int iM = 0; iM < 2; ++iM) {
      TCanvas c(Form("f%sCompare_%.0f_%.0f",kAntimatterMatter[iM],kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1]),Form("f%sCompare_%.0f_%.0f",kAntimatterMatter[iM],kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1]));
      TLegend l(0.137845,0.712281,0.887218,0.877193);
      c.cd();
      for (int iFit = 0; iFit < 2; ++iFit){
        hPrimary[iM][iFit] = (TH1D*)InFile[iFit]->Get(Form("f%sPrimFrac_%.0f_%.0f",kAntimatterMatter[iM],kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1])); 
        hPrimary[iM][iFit]->GetXaxis()->SetRangeUser(1.,2.);
        hPrimary[iM][iFit]->GetYaxis()->SetRangeUser(0.75,.95);
        hPrimary[iM][iFit]->Draw("");
        hPrimary[iM][iFit]->SetLineColor(kHistColor[iFit]);
        hPrimary[iM][iFit]->SetMarkerColor(kHistColor[iFit]);
        if (iFit == 1)
          hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->SetLineStyle(7);
        hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->SetLineColor(kHistColor[iFit]);
        hPrimary[iM][iFit]->SetMarkerStyle(kMarkerSize[iFit]);
        hPrimary[iM][iFit]->SetMarkerSize(0.9);
        double chi2 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->GetChisquare();
        int ndf = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->GetNDF();
        double p0 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->GetParameter(0);
        double errp0 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->GetParError(0);
        double p1 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->GetParameter(1);
        double errp1 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsProton[iC][0], kCentBinsLimitsProton[iC][1]))->GetParError(1);
        l.AddEntry(hPrimary[iM][iFit],Form("%s, #chi^{2}/NDF = %.2f/%.d, p_{0} = %.4f #pm %.4f, p_{1} = %.4f #pm %.4f",kFittingMethod[iFit],chi2,ndf,p0,errp0,p1,errp1));
      }
      TH1D hPrimaryRatio(*hPrimary[iM][0]);
      hPrimaryRatio.SetName(Form("f%sPrimaryRatio_%.0f_%.0f",kAntimatterMatter[iM],kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1]));
      hPrimaryRatio.Divide(hPrimary[iM][1],hPrimary[iM][0]); 
      hPrimaryRatio.GetYaxis()->SetTitle("RooFit/TFF");
      hPrimaryRatio.Fit("pol0");
      hPrimaryRatio.Write();
      hPrimary[iM][0]->Draw("");
      hPrimary[iM][1]->Draw("same");
      l.Draw("same");
      c.Write();
      c.Print(Form("plots/%s.pdf",hPrimary[iM][0]->GetName()));
    }
    TH1D hDivideTFF(*hPrimary[0][0]);
    hDivideTFF.SetName(Form("fDivideTFF_%.0f_%.0f",kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1]));
    hDivideTFF.Divide(hPrimary[1][0],hPrimary[0][0]);
    hDivideTFF.GetYaxis()->SetTitle("#it{f}_{prim}^{antip}/#it{f}_{prim}^{p}");
    hDivideTFF.Fit("pol0");
    hDivideTFF.Write();
    TH1D hDivideRooFit(*hPrimary[0][0]);
    hDivideRooFit.SetName(Form("fDivideRooFit_%.0f_%.0f",kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1]));
    hDivideRooFit.Divide(hPrimary[1][1],hPrimary[0][1]);
    hDivideRooFit.GetYaxis()->SetTitle("#it{f}_{prim}^{antip}/#it{f}_{prim}^{p}");
    hDivideRooFit.Fit("pol0");
    hDivideRooFit.Write();
    TCanvas cCompareDivide("cCompareDivide","cCompareDivide");
    TLegend lCompareDivide(0.555138,0.710526,0.779449,0.864912);
    cCompareDivide.cd();
    hDivideTFF.GetYaxis()->SetRangeUser(0.98,1.03);
    hDivideTFF.Draw("");
    hDivideRooFit.SetMarkerStyle(24);
    hDivideRooFit.SetMarkerColor(kBlue);
    hDivideRooFit.SetLineColor(kBlue);
    hDivideRooFit.GetFunction("pol0")->SetLineColor(kBlue);
    hDivideRooFit.GetFunction("pol0")->SetLineStyle(7);
    hDivideRooFit.Draw("same");
    lCompareDivide.AddEntry(&hDivideTFF,"TFF");
    lCompareDivide.AddEntry(&hDivideRooFit,"RooFit");
    lCompareDivide.Draw("same");
    cCompareDivide.Write();
    // plot ratios
    TCanvas cRatio("cRatio","cRatio");
    TLegend lRatio(0.666667,0.701754,0.890977,0.85614);
    TH1D *hRatioTFF = (TH1D*)inFileSpectraTFF.Get(Form("fRatio_%.0f_%.0f",kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1]));
    TH1D *hRatioRoo = (TH1D*)inFileSpectraRoo.Get(Form("fRatio_%.0f_%.0f",kCentBinsLimitsProton[iC][0],kCentBinsLimitsProton[iC][1]));
    hRatioTFF->SetMarkerSize(0.9);
    hRatioTFF->SetMarkerColor(kRed);
    hRatioTFF->SetLineColor(kRed);
    hRatioTFF->GetYaxis()->SetRangeUser(0.90,1.10);
    hRatioTFF->GetXaxis()->SetRangeUser(1.,2.);
    hRatioTFF->Draw("");
    hRatioRoo->SetMarkerStyle(24);
    hRatioRoo->SetMarkerSize(0.9);
    hRatioRoo->SetLineColor(kBlue);
    hRatioRoo->SetMarkerColor(kBlue);
    hRatioRoo->GetFunction("pol0")->SetLineColor(kBlue);
    hRatioRoo->GetFunction("pol0")->SetLineStyle(7);
    hRatioRoo->Draw("same");
    lRatio.AddEntry(hRatioTFF,"TFF");
    lRatio.AddEntry(hRatioRoo,"RooFit");
    lRatio.Draw("same");
    TLatex ratioTFF(1.1,1.08,Form("R_{TFF} = %.4f #pm %.4f",hRatioTFF->GetFunction("pol0")->GetParameter(0),hRatioTFF->GetFunction("pol0")->GetParError(0)));
    TLatex ratioRoo(1.1,1.06,Form("R_{RooFit} = %.4f #pm %.4f",hRatioRoo->GetFunction("pol0")->GetParameter(0),hRatioRoo->GetFunction("pol0")->GetParError(0)));
    ratioTFF.Draw("same");
    ratioRoo.Draw("same");
    cRatio.Write();
  }
  OutFileCompare.Close();
}