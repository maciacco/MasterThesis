// PlotRatioCheck.cpp

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

void PlotRatioCheck(const char *cutSettings="", const char *outFileName = "RatioCheckPlot", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFileRatioName = "SpectraProtonNoPrimary", const char *inFilePrimaryName = "PrimaryProton")
{
  gStyle->SetOptStat(0);

  gSystem->Exec(Form("mkdir %s/ratio_check_plots", kPlotDir));

  TFile *inFileRatio = TFile::Open(Form("%s/%s.root", kOutDir, inFileRatioName));
  if (!inFileRatio)
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
    TCanvas cRatioCheck(Form("cRatioCheck_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), Form("cRatioCheck_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    cRatioCheck.SetTicks(1, 1);
    TLegend lRatioCheck(0.484241, 0.208511+0.2, 0.795129, 0.389362+0.2);
    lRatioCheck.SetTextSize(0.035);
    lRatioCheck.SetBorderSize(0);
    auto fRatioCheck = (TH1D *)inFileRatio->Get(Form("./fRatio_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    fRatioCheck->SetTitle(Form("%.0f-%.0f%%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    auto f = fRatioCheck->GetFunction("pol0");
    f->Delete("");
    fRatioCheck->GetYaxis()->SetRangeUser(0., 1.2);
    fRatioCheck->GetXaxis()->SetRangeUser(0.0, 5.0);
    fRatioCheck->GetYaxis()->SetTitle("#it{f_{prim}}");
    fRatioCheck->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fRatioCheck->SetMarkerStyle(24);
    fRatioCheck->SetMarkerSize(0.8);
    fRatioCheck->SetLineColor(kBlack);
    fRatioCheck->SetMarkerColor(kBlack);
    fRatioCheck->Draw("pe");
    lRatioCheck.AddEntry(fRatioCheck, "#frac{[#it{N_{raw}} / (#epsilon #times Acc.)]_{ #bar{d}}}{[#it{N_{raw}} / (#epsilon #times Acc.)]_{ d}}");
  
    auto fPrimaryFraction = (TH1D*)inFilePrimary->Get(Form("fMPrimFrac_%.0f_%.0f", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    fPrimaryFraction->SetMarkerColor(kRed);
    fPrimaryFraction->SetLineColor(kRed);
    fPrimaryFraction->SetMarkerStyle(20);
    fPrimaryFraction->Draw("same");
    lRatioCheck.AddEntry(fPrimaryFraction, "#it{f_{prim}} from TFF");

    lRatioCheck.Draw("same");
    cRatioCheck.Write();
    cRatioCheck.Print(Form("%s/ratio_check_plots/%s.png", kPlotDir, cRatioCheck.GetName()));
  }

  outFile.Close();
}
