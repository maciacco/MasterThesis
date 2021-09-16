// PlotCentralitySelection.cpp

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TSystem.h>

#include "utils/Utils.h"
#include "utils/Config.h"

using utils::TTList;

void PlotCentralitySelection(const char *outFileName = "CentralitySelectionPlot", const char *inFileName = "AnalysisResults")
{
  gStyle->SetOptStat(0);

  TFile *inFile = TFile::Open(Form("data/Deuteron_PbPb/%s.root", inFileName));
  if (!inFile)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }
  auto inList = (TTList *)inFile->Get("mpuccio_deuterons_");
  auto fNormalisationHist = (TH2F *)inList->Get("fNormalisationHist");
  TH1D *fCent = fNormalisationHist->ProjectionX("fCent", 4, 4);

  TFile outFile(Form("%s.root", outFileName), "recreate");
  TCanvas cCent("cCent_LHC18qr", "cCent_LHC18qr");
  cCent.SetTicks(1, 1);
  fCent->SetDrawOption("histo");
  fCent->GetXaxis()->SetTitle("Centrality (%)");
  fCent->GetXaxis()->SetTitleSize(0.05);
  fCent->GetXaxis()->SetRangeUser(0,100);
  fCent->GetYaxis()->SetTitle("Events");
  fCent->SetLineWidth(2);
  fCent->SetLineColor(kRed);
  fCent->SetFillStyle(3345);
  fCent->SetFillColor(kRed);
  fCent->Draw("");
  cCent.Write();
  cCent.Print("cCent_LHC18qr.pdf");

  outFile.Close();
}
