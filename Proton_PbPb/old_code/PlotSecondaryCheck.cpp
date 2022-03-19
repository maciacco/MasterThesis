// PlotSecondaryCheck.cpp

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TSystem.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using utils::TTList;
using namespace proton;

void PlotSecondaryCheck(const char *cutSettings="", const char *outFileName = "RatioCheckPlot", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFilePrimAntiName = "PrimaryProton_2018_NoActiveLengthCut_AntiProtonsAsPrimaries", const char *inFilePrimaryName = "PrimaryProton_2018_NoActiveLengthCut")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gSystem->Exec(Form("mkdir %s/prim_check_plots", kPlotDir));

  TFile *inFilePrimAnti = TFile::Open(Form("%s/%s.root", kOutDir, inFilePrimAntiName));
  if (!inFilePrimAnti)
  {
    std::cout << "Input file does not exist!" << std::endl;
    return;
  }
    
  TFile *inFilePrimary = TFile::Open(Form("%s/%s.root", kOutDir, inFilePrimaryName));
  if (!inFilePrimary)
  {
    std::cout << "Input file does not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  outFile.mkdir(histoNameDir);

  for (int iCent = 0; iCent < kNCentClasses; ++iCent)
  {
    TCanvas cPrimCheck(Form("cPrimCheck_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), Form("cPrimCheck_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    cPrimCheck.SetTicks(1, 1);
    TLegend lPrimCheck(0.484241, 0.208511+0.1, 0.795129, 0.389362+0.1);
    lPrimCheck.SetTextSize(0.035);
    lPrimCheck.SetBorderSize(0);
    
    auto fPrimaryFraction = (TH1D*)inFilePrimary->Get(Form("fMPrimFrac_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    fPrimaryFraction->SetTitle(Form("%.0f-%.0f%%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    auto f2 = fPrimaryFraction->GetFunction(Form("fMFunctionFit_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    f2->SetLineStyle(7);
    f2->SetLineWidth(3);
    f2->SetLineColor(kGray+1);
    fPrimaryFraction->GetYaxis()->SetTitle("#it{f_{prim}}");
    fPrimaryFraction->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fPrimaryFraction->GetYaxis()->SetRangeUser(0., 1.2);
    fPrimaryFraction->GetXaxis()->SetRangeUser(0.8, 1.6);
    fPrimaryFraction->SetMarkerColor(kGray+1);
    fPrimaryFraction->SetLineColor(kGray+1);
    fPrimaryFraction->SetMarkerStyle(24);
    fPrimaryFraction->SetMarkerSize(1.0);
    fPrimaryFraction->Draw("pe");
    lPrimCheck.AddEntry(fPrimaryFraction, "primaries from MC");

    auto fPrimAnti = (TH1D *)inFilePrimAnti->Get(Form("fMPrimFrac_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    auto f = fPrimAnti->GetFunction(Form("fMFunctionFit_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    f->SetLineStyle(7);
    f->SetLineWidth(3);
    f->SetLineColor(kRed);
    fPrimAnti->SetMarkerStyle(24);
    fPrimAnti->SetMarkerSize(1.0);
    fPrimAnti->SetLineColor(kRed);
    fPrimAnti->SetMarkerColor(kRed);
    fPrimAnti->Draw("same");
    lPrimCheck.AddEntry(fPrimAnti, "antiProton data as primaries");

    lPrimCheck.Draw("same");
    cPrimCheck.Write();
    cPrimCheck.Print(Form("%s/prim_check_plots/%s.png", kPlotDir, cPrimCheck.GetName()));
  }

  outFile.Close();
}
