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

void PlotEfficienciesDifferential(const char *cutSettings="", const char *outFileName = "EfficiencyPlotsDifferential", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFile = "EfficiencyProton")
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

  Color_t colors[3] = {kRed, kBlue, kGreen + 2};

  for (int iMatt = 0; iMatt < 2; iMatt++)
  {
    // efficiency plots in centrality classes
    TCanvas cEff(Form("cEff_%s", kAntimatterMatter[iMatt]), Form("cEff_%s", kAntimatterMatterLabel[iMatt]));
    cEff.SetTicks(1, 1);
    TLegend lEff(0.484241, 0.208511+0.4, 0.795129, 0.389362+0.4);
    lEff.SetHeader(Form("%s, ITS + TPC + TOF",kAntimatterMatterLabel[iMatt]));
    lEff.SetTextSize(0.035);
    lEff.SetBorderSize(0);
    
    // efficiency comparison to MB
    TCanvas cEffCompare(Form("cEffCompare_%s", kAntimatterMatter[iMatt]), Form("cEffCompare_%s", kAntimatterMatterLabel[iMatt]));
    cEffCompare.SetTicks(1, 1);
    TLegend lEffCompare(0.484241, 0.208511+0.47, 0.795129, 0.389362+0.47);
    lEffCompare.SetHeader(Form("%s, ITS + TPC + TOF",kAntimatterMatterLabel[iMatt]));
    lEffCompare.SetTextSize(0.035);
    lEffCompare.SetBorderSize(0);

    TH1D *fEff[kNCentClasses];
    TH1D *fEffMB = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_0_90", kAntimatterMatter[iMatt]));
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      fEff[iCent] = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fEff[iCent]->SetTitle("");
      fEff[iCent]->GetYaxis()->SetRangeUser(0., 1.1);
      fEff[iCent]->GetXaxis()->SetRangeUser(1., 2.0);
      fEff[iCent]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fEff[iCent]->GetYaxis()->SetTitle("#epsilon #times A");
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
      cEffCompare.cd();
      TH1D *hTemp = (TH1D*)fEff[iCent]->Clone();
      hTemp->Divide(fEffMB);
      hTemp->GetYaxis()->SetRangeUser(0.7, 1.3);
      if (iCent == 0)
        hTemp->Draw("pe");
      else
        hTemp->Draw("pesame");
      hTemp->GetYaxis()->SetTitle("Efficiency ratio to MB");
      //hTemp->SetTitle("LHC16h7_nucleiInjected_LHC15o");
      lEffCompare.AddEntry(hTemp, Form("%0.f - %0.f %%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
    }
    
    cEff.cd();
    lEff.Draw("same");
    cEff.Write();
    cEff.Print(Form("%s/efficiency_plots/%s.pdf", kPlotDir, cEff.GetName()));

    cEffCompare.cd();
    lEffCompare.Draw("same");
    cEffCompare.Write();
    cEffCompare.Print(Form("%s/efficiency_plots/%s.pdf", kPlotDir, cEffCompare.GetName()));
  }

  outFile.Close();
}
