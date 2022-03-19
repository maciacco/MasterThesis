// PlotEfficienciesDifferential.cpp

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

void PlotEfficienciesDifferential1(const char *cutSettings = "", const char *outFileName = "EfficiencyPlotsDifferential1", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFile = "EfficiencyProton")
{
  gStyle->SetOptStat(0000000000000);

  gSystem->Exec(Form("mkdir %s/efficiency_plots", kPlotDir));

  TFile *inFileEff = TFile::Open(Form("%s/%s.root", kOutDir, inFile));
  if (!inFileEff)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  outFile.mkdir(histoNameDir);
  TH1D *hTemp[2][3]{nullptr};
  Color_t colors[4] = {kRed, kBlue, kGreen + 2, kBlack};

  for (int iMatt = 0; iMatt < 2; iMatt++)
  {
    // efficiency plots in centrality classes
    TCanvas cEff(Form("cEff_%s", kAntimatterMatter[iMatt]), Form("cEff_%s", kAntimatterMatterLabel[iMatt]));
    cEff.SetTicks(1, 1);
    TLegend lEff(0.484241, 0.208511 + 0.4, 0.795129, 0.389362 + 0.4);
    lEff.SetHeader(Form("%s, ITS + TPC + TOF", kAntimatterMatterLabel[iMatt]));
    lEff.SetTextSize(0.035);
    lEff.SetBorderSize(0);

    // efficiency comparison to MB
    TCanvas cEffCompare(Form("cEffCompare_%s", kAntimatterMatter[iMatt]), Form("cEffCompare_%s", kAntimatterMatterLabel[iMatt]));
    cEffCompare.SetTicks(1, 1);
    TLegend lEffCompare(0.484241, 0.208511 + 0.47, 0.795129, 0.389362 + 0.47);
    lEffCompare.SetHeader(Form("%s, ITS + TPC + TOF", kAntimatterMatterLabel[iMatt]));
    lEffCompare.SetTextSize(0.035);
    lEffCompare.SetBorderSize(0);

    TH1D *fEff[kNCentClasses];
    TH1D *fEffMB = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_0_90", kAntimatterMatter[iMatt]));
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      fEff[iCent] = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fEff[iCent]->SetTitle("");
      fEff[iCent]->GetYaxis()->SetRangeUser(0., 1.1);
      fEff[iCent]->GetXaxis()->SetRangeUser(0.0, 6.0);
      fEff[iCent]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fEff[iCent]->SetMarkerStyle(24);
      fEff[iCent]->SetMarkerSize(0.8);
      fEff[iCent]->SetLineColor(colors[iCent]);
      fEff[iCent]->SetMarkerColor(colors[iCent]);
      cEff.cd();
      if (iCent == 0)
        fEff[iCent]->Draw("pe");
      else
        fEff[iCent]->Draw("pesame");
      lEff.AddEntry(fEff[iCent], Form("%0.f - %0.f %%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));

      // compare with MB
      hTemp[iMatt][iCent] = (TH1D *)fEff[iCent]->Clone();
      hTemp[iMatt][iCent]->Divide(fEffMB);
      hTemp[iMatt][iCent]->GetYaxis()->SetRangeUser(0.7, 1.3);
      hTemp[iMatt][iCent]->GetYaxis()->SetTitle("#frac{#bar{p} ratio to MB}{p ratio to MB}");
      hTemp[iMatt][iCent]->SetTitle("LHC18qr_pass3_anchored");
      lEffCompare.AddEntry(hTemp[iMatt][iCent], Form("%0.f - %0.f %%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));

      if (iMatt == 1)
      {
        cEffCompare.cd();
        hTemp[0][iCent]->Divide(hTemp[1][iCent]);
        if (iCent == 1)
          hTemp[0][iCent]->Draw("pe");
        else
          hTemp[0][iCent]->Draw("pesame");
        lEffCompare.Draw("same");
      }
    }

    if (iMatt == 1)
    {
      cEffCompare.Write();
      cEffCompare.Print(Form("%s/efficiency_plots/%s.pdf", kPlotDir, cEffCompare.GetName()));
    }
    cEff.cd();
    lEff.Draw("same");
    cEff.Write();
    cEff.Print(Form("%s/efficiency_plots/%s.pdf", kPlotDir, cEff.GetName()));
  }

  outFile.Close();
}
