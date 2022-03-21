#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include "../../utils/Config.h"

using namespace pion;
const int kMarkerSize[] = {20,24};
const char* kFittingMethod[] = {"TFF","RooFit"};
const Color_t kHistColor[] = {kRed,kBlue};

void PlotSecondaryCorrection(const char *inFileRooFitName = "PrimaryPionRoo", const char *inFileTFFName = "PrimaryPion", const char *outFileCompareName = "PrimaryCompareRooFitTFF", const char *inFileSpectraTFFName = "SpectraPion_MC21l5_raw", const char* inFileSpectraRooName = "SpectraPion_MC21l5_raw_ROO"){
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  // compare RooFit and ROOT::TFractionFitter
  TFile *InFile[2];
  InFile[0] = new TFile(Form("../out/%s.root",inFileTFFName));
  InFile[1] = new TFile(Form("../out/%s.root",inFileRooFitName));
  TFile OutFileCompare(Form("%s.root",outFileCompareName),"recreate");
  for (int iC = 0; iC < kNCentClasses; ++iC){
    TH1D *hPrimary[2][2];
    for (int iM = 0; iM < 2; ++iM) {
      TCanvas c(Form("f%sCompare_%.0f_%.0f",kAntimatterMatter[iM],kCentBinsLimitsPion[iC][0],kCentBinsLimitsPion[iC][1]),Form("f%sCompare_%.0f_%.0f",kAntimatterMatter[iM],kCentBinsLimitsPion[iC][0],kCentBinsLimitsPion[iC][1]));
      TLegend l(0.137845,0.712281,0.887218,0.877193);
      c.cd();
      for (int iFit = 0; iFit < 2; ++iFit){
        hPrimary[iM][iFit] = (TH1D*)InFile[iFit]->Get(Form("f%sPrimFrac_%.0f_%.0f",kAntimatterMatter[iM],kCentBinsLimitsPion[iC][0],kCentBinsLimitsPion[iC][1])); 
        hPrimary[iM][iFit]->GetXaxis()->SetRangeUser(0.5,1.25);
        hPrimary[iM][iFit]->GetYaxis()->SetRangeUser(0.96,1.0);
        hPrimary[iM][iFit]->Draw("");
        hPrimary[iM][iFit]->SetLineColor(kHistColor[iFit]);
        hPrimary[iM][iFit]->SetMarkerColor(kHistColor[iFit]);
        if (iFit == 1)
          hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->SetLineStyle(7);
        hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->SetLineColor(kHistColor[iFit]);
        hPrimary[iM][iFit]->SetMarkerStyle(kMarkerSize[iFit]);
        hPrimary[iM][iFit]->SetMarkerSize(0.9);
        double chi2 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->GetChisquare();
        int ndf = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->GetNDF();
        double p0 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->GetParameter(0);
        double errp0 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->GetParError(0);
        double p1 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->GetParameter(1);
        double errp1 = hPrimary[iM][iFit]->GetFunction(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iM], kCentBinsLimitsPion[iC][0], kCentBinsLimitsPion[iC][1]))->GetParError(1);
        l.AddEntry(hPrimary[iM][iFit],Form("%s, #chi^{2}/NDF = %.2f/%.d, p_{0} = %.4f #pm %.4f, p_{1} = %.4f #pm %.4f",kFittingMethod[iFit],chi2,ndf,p0,errp0,p1,errp1));
      }

      hPrimary[iM][0]->Draw("");
      hPrimary[iM][1]->Draw("same");
      l.Draw("same");
      c.Write();
      c.Print(Form("../plots/%s.pdf",hPrimary[iM][0]->GetName()));
    }
  }
  OutFileCompare.Close();
}