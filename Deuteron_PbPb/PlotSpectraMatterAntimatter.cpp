// PlotSpectraMatterAntimatter.cpp

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
using namespace deuteron;

void PlotSpectraMatterAntimatter(const char *cutSettings = "", const char *outFileName = "SpectraMatterAntimatter", const char *histoNameDir = "", const char *outFileOption = "recreate", const char *inFile = "SpectraDeuteron")
{
  gStyle->SetOptStat(0);

  gSystem->Exec(Form("mkdir %s/spectra_AM_plot", kPlotDir));

  TFile *inFileSpectra = TFile::Open(Form("%s/%s.root", kOutDir, inFile));
  if (!inFileSpectra)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  TH1D *spectra[2];
  int specColor[]={kRed, kGreen+2};
  for (int iC = 0; iC < 3; ++iC)
  {
    TCanvas c(Form("fSpectraCompareAntimatterMatter_%.0f_%.0f", kCentBinsLimitsDeuteron[iC][0], kCentBinsLimitsDeuteron[iC][1]),Form("fSpectraCompareAntimatterMatter_%.0f_%.0f", kCentBinsLimitsDeuteron[iC][0], kCentBinsLimitsDeuteron[iC][1]));
    TLegend l(0.575188, 0.696491, 0.913534, 0.861404);
    l.SetHeader(Form("Centrality %.0f-%.0f%%", kCentBinsLimitsDeuteron[iC][0], kCentBinsLimitsDeuteron[iC][1]));
    c.cd();
    for (int iMatt = 0; iMatt < 2; ++iMatt)
    {
      spectra[iMatt]=(TH1D*)inFileSpectra->Get(Form("_1_1_1/f%sSpectra_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iC][0], kCentBinsLimitsDeuteron[iC][1]));
      spectra[iMatt]->SetMarkerStyle(20);
      spectra[iMatt]->SetMarkerSize(0.9);
      spectra[iMatt]->SetMarkerColor(specColor[iMatt]);
      spectra[iMatt]->SetLineColor(specColor[iMatt]);
      spectra[iMatt]->SetTitle("");
      spectra[iMatt]->GetXaxis()->SetRangeUser(1., 5.);
      std::string drawOption = "pe";
      if (iMatt == 1){
        drawOption = drawOption + "same";
      }
      spectra[iMatt]->Draw(drawOption.data());
      l.AddEntry(spectra[iMatt],kAntimatterMatterLabelExtended[iMatt]);
    }
    l.Draw("same");
    c.Write();
    c.Print(Form("%s/spectra_AM_plot/%s.png", kPlotDir, c.GetName()));
  }
  outFile.Close();
}
