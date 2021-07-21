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
using namespace proton;

void PlotEfficiencies(const float cutDCAz = 1.f, const int cutTPCcls = 89, const char *outFileName = "EfficiencyPlots", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFile = "EfficiencyProton")
{
  gStyle->SetOptFit(0);
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

  Color_t colors[3] = {kRed, kBlack};

  TCanvas cEff("cEff", "cEff");
  cEff.SetTicks(1, 1);
  TLegend lEff(0.484241, 0.208511+0.4, 0.795129, 0.389362+0.4);
  lEff.SetHeader("ITS + TPC + TOF");
  lEff.SetTextSize(0.035);
  lEff.SetBorderSize(0);
  TH1D *fEff[kNCentClasses];
  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    fEff[iMatt] = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_0_90", kAntimatterMatter[iMatt]));
    fEff[iMatt]->SetTitle("");
    fEff[iMatt]->GetYaxis()->SetRangeUser(0., 1.1);
    fEff[iMatt]->GetXaxis()->SetRangeUser(0.0, 6.0);
    fEff[iMatt]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fEff[iMatt]->SetMarkerStyle(20);
    fEff[iMatt]->SetMarkerSize(0.8);
    fEff[iMatt]->SetLineColor(colors[iMatt]);
    fEff[iMatt]->SetMarkerColor(colors[iMatt]);
    if (iMatt == 0)
      fEff[iMatt]->Draw("pe");
    else
      fEff[iMatt]->Draw("pesame");
    lEff.AddEntry(fEff[iMatt], Form("%s", kAntimatterMatterLabelExtended[iMatt]));
  }
  lEff.Draw("same");
  cEff.Write();
  cEff.Print(Form("%s/efficiency_plots/%s.png", kPlotDir, cEff.GetName()));

  outFile.Close();
}
