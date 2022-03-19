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
#include <TLine.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using utils::TTList;
using namespace proton;

void PlotEfficienciesDifferential2(const char *cutSettings = "", const char *outFileName = "EfficiencyPlotsDifferential2_20g7", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFile = "EfficiencyProton_LongMCTracks")
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
  TH1D *fEff[2][4]{nullptr};
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

    TH1D *fEffMB = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_0_90", kAntimatterMatter[iMatt]));
    for (int iCent = 3; iCent > -1; --iCent)
    {
      fEff[iMatt][iCent] = (TH1D *)inFileEff->Get(Form("f%sEff_TOF_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fEff[iMatt][iCent]->SetTitle("");
      fEff[iMatt][iCent]->GetYaxis()->SetRangeUser(0., 1.1);
      fEff[iMatt][iCent]->GetXaxis()->SetRangeUser(0.8, 5.0);
      fEff[iMatt][iCent]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fEff[iMatt][iCent]->SetMarkerStyle(24);
      fEff[iMatt][iCent]->SetMarkerSize(0.8);
      fEff[iMatt][iCent]->SetLineColor(colors[iCent]);
      fEff[iMatt][iCent]->SetMarkerColor(colors[iCent]);

      // compare with MB
      hTemp[iMatt][iCent] = (TH1D *)fEff[iMatt][iCent]->Clone();
      hTemp[iMatt][iCent]->Divide(fEffMB);
      hTemp[iMatt][iCent]->GetYaxis()->SetRangeUser(0.96, 1.04);
      hTemp[iMatt][iCent]->GetXaxis()->SetRangeUser(0.8, 6.);
      hTemp[iMatt][iCent]->GetYaxis()->SetTitle("#frac{#bar{p} ratio to MB}{p ratio to MB}");
      hTemp[iMatt][iCent]->SetTitle("LHC18qr_pass3_anchored");
      
      if (iMatt == 1)
      {
        fEff[0][iCent]->Divide(fEff[1][iCent]);
        fEff[0][iCent]->GetYaxis()->SetTitle("(#epsilon#times A)_{antiproton}/(#epsilon#times A)_{proton}");
        fEff[0][iCent]->GetYaxis()->SetTitleSize(0.035);
        fEff[0][iCent]->GetYaxis()->SetTitleOffset(1.5);
        cEff.cd();
        if (iCent == 3)
          fEff[0][iCent]->Draw("pe");
        else
          fEff[0][iCent]->Draw("pesame");
        lEff.AddEntry(fEff[0][iCent], Form("%0.f - %0.f %%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      }

      if (iCent < 3) {
        auto fEff1=(TH1D*)fEff[0][iCent]->Clone();
        fEff1->Divide(fEff[0][3]);
        cEffCompare.cd();
        if (iCent == 2)
          fEff1->Draw("pe");
        else
          fEff1->Draw("pesame");
        fEff1->GetYaxis()->SetTitle("#frac{(#epsilon#times A)_{antiproton}}{(#epsilon#times A)_{proton}} ratio to MB");
        lEffCompare.AddEntry(fEff1, Form("%0.f - %0.f %%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
        
        TLine l(1.0,1.,6.0,1.);
        l.SetLineStyle(kDashed);
        l.Draw("same");
        lEffCompare.Draw("same");

      }
    }

    if (iMatt == 1)
    {
      cEff.cd();
      lEff.Draw("same");
      cEff.Write();
      cEff.Print(Form("%s/efficiency_plots/%s.png", kPlotDir, cEff.GetName()));

      cEffCompare.Write();
      cEffCompare.Print(Form("%s/efficiency_plots/%s.png", kPlotDir, cEffCompare.GetName()));
    }
  }

  outFile.Close();
}
