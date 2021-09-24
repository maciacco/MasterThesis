// PlotEfficiencies.cpp

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
using namespace he3;

void PlotEfficiencies(const float cutDCAz = 1.f, const int cutTPCcls = 89, const char *outFileName = "EfficiencyPlots", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFile = "EfficiencyHe3")
{
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  gSystem->Exec(Form("mkdir %s/efficiency_plots", kPlotDir));

  TFile *inFileEff = TFile::Open(Form("%s/%s.root", kOutDir, inFile));
  if (!inFileEff)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  outFile.mkdir(histoNameDir);

  Color_t colors[3] = {kRed, kBlue, kGreen + 2};

  for (int iMatt = 0; iMatt < 2; iMatt++)
  {
    TCanvas cEff(Form("cEff_%s", kAntimatterMatter[iMatt]), Form("cEff_%s", kAntimatterMatterLabel[iMatt]));
    cEff.SetTicks(1, 1);
    TLegend lEff(0.484241, 0.208511, 0.795129, 0.389362);
    lEff.SetHeader(Form("%s, ITS + TPC",kAntimatterMatterLabel[iMatt]));
    lEff.SetTextSize(0.035);
    lEff.SetBorderSize(0);
    TH1D *fEff[kNCentClasses];
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      fEff[iCent] = (TH1D *)inFileEff->Get(Form("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fEff[iCent]->SetTitle("");
      fEff[iCent]->GetYaxis()->SetRangeUser(0., 1.1);
      fEff[iCent]->GetXaxis()->SetRangeUser(1.0, 8.0);
      fEff[iCent]->SetMarkerStyle(24);
      fEff[iCent]->SetMarkerSize(0.8);
      fEff[iCent]->SetLineColor(colors[iCent]);
      fEff[iCent]->SetMarkerColor(colors[iCent]);
      if (iCent == 0)
        fEff[iCent]->Draw("pe");
      else
        fEff[iCent]->Draw("pesame");
      lEff.AddEntry(fEff[iCent], Form("%0.f - %0.f %%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
    }
    lEff.Draw("same");
    cEff.Write();
    cEff.Print(Form("%s/efficiency_plots/%s.pdf", kPlotDir, cEff.GetName()));
  }

  outFile.Close();
}
