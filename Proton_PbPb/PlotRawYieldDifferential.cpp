// PlotRawYieldicienciesDifferential.cpp

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

void PlotRawYieldDifferential(const char *cutSettings = "", const char *outFileName = "RawYieldPlotsDifferential1", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFile = "SignalProtonGausDExpSignal1_finePtBinning")
{
  gStyle->SetOptStat(0000000000000);

  gSystem->Exec(Form("mkdir %s/rawYield_plots", kPlotDir));

  TFile *inFileRaw = TFile::Open(Form("%s/%s.root", kOutDir, inFile));
  if (!inFileRaw)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  outFile.mkdir(histoNameDir);
  TH1D *fRaw[2][4]{nullptr};
  Color_t colors[4] = {kRed, kBlue, kGreen + 2, kBlack};

  for (int iMatt = 0; iMatt < 2; iMatt++)
  {
    // raw yield plots in centrality classes
    TCanvas cRaw(Form("cRaw_%s", kAntimatterMatter[iMatt]), Form("cRaw_%s", kAntimatterMatterLabel[iMatt]));
    cRaw.SetTicks(1, 1);
    TLegend lRaw(0.484241, 0.208511 + 0.4, 0.795129, 0.389362 + 0.4);
    lRaw.SetHeader(Form("%s, Raw yields (TOF)", kAntimatterMatterLabel[iMatt]));
    lRaw.SetTextSize(0.035);
    lRaw.SetBorderSize(0);

    // raw yield comparison to MB
    TCanvas cRawCompare(Form("cRawCompare_%s", kAntimatterMatter[iMatt]), Form("cRawCompare_%s", kAntimatterMatterLabel[iMatt]));
    cRawCompare.SetTicks(1, 1);
    TLegend lRawCompare(0.484241, 0.208511 + 0.47, 0.795129, 0.389362 + 0.47);
    lRawCompare.SetHeader(Form("%s, Raw yields (TOF)", kAntimatterMatterLabel[iMatt]));
    lRawCompare.SetTextSize(0.035);
    lRawCompare.SetBorderSize(0);

    for (int iCent = 3; iCent > -1; --iCent)
    {
      fRaw[iMatt][iCent] = (TH1D *)inFileRaw->Get(Form("_1_1/f%sTOFrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fRaw[iMatt][iCent]->SetTitle("");
      fRaw[iMatt][iCent]->GetYaxis()->SetRangeUser(0., 1.1);
      fRaw[iMatt][iCent]->GetXaxis()->SetRangeUser(0.8, 5.0);
      fRaw[iMatt][iCent]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fRaw[iMatt][iCent]->SetMarkerStyle(24);
      fRaw[iMatt][iCent]->SetMarkerSize(0.8);
      fRaw[iMatt][iCent]->SetLineColor(colors[iCent]);
      fRaw[iMatt][iCent]->SetMarkerColor(colors[iCent]);
      
      if (iMatt == 1)
      {
        fRaw[0][iCent]->Divide(fRaw[1][iCent]);
        fRaw[0][iCent]->GetYaxis()->SetTitle("(rawYield)_{antiproton}/(rawYield)_{proton}");
        cRaw.cd();
        if (iCent == 3)
          fRaw[0][iCent]->Draw("pe");
        else
          fRaw[0][iCent]->Draw("pesame");
        fRaw[0][iCent]->GetYaxis()->SetTitleSize(0.035);
        fRaw[0][iCent]->GetYaxis()->SetTitleOffset(1.5);
        lRaw.AddEntry(fRaw[0][iCent], Form("%0.f - %0.f %%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      }

      if (iCent < 3) {
        auto fRaw1=(TH1D*)fRaw[0][iCent]->Clone();
        fRaw1->Divide(fRaw[0][3]);
        cRawCompare.cd();
        if (iCent == 2)
          fRaw1->Draw("pe");
        else
          fRaw1->Draw("pesame");
        fRaw1->GetYaxis()->SetTitle("#frac{(RawYield)_{antiproton}}{(RawYield)_{proton}} ratio to MB");
        lRawCompare.AddEntry(fRaw1, Form("%0.f - %0.f %%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
        lRawCompare.Draw("same");
      }
    }

    if (iMatt == 1)
    {
      cRaw.cd();
      lRaw.Draw("same");
      cRaw.Write();
      cRaw.Print(Form("%s/rawYield_plots/%s.png", kPlotDir, cRaw.GetName()));

      cRawCompare.Write();
      cRawCompare.Print(Form("%s/rawYield_plots/%s.png", kPlotDir, cRawCompare.GetName()));
    }
  }

  outFile.Close();
}
