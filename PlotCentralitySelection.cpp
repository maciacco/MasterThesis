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

void PlotCentralitySelection(const char *outFileName = "CentralitySelectionPlotProton", const char *inFileName = "AnalysisResults_LHC18qr")
{
  gStyle->SetOptStat(0);

  TFile *inFile_q = TFile::Open(Form("data/He3_PbPb/%s_q.root", inFileName));
  TFile *inFile_r = TFile::Open(Form("data/He3_PbPb/%s_r.root", inFileName));
  if (!inFile_q || !inFile_r)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }
  auto inList_q = (TTList *)inFile_q->Get("mpuccio_he3_");
  auto inList_r = (TTList *)inFile_r->Get("mpuccio_he3_");
  auto fNormalisationHist_q = (TH2F *)inList_q->Get("fNormalisationHist");
  auto fNormalisationHist_r = (TH2F *)inList_r->Get("fNormalisationHist");
  TH1D *fCent_q = fNormalisationHist_q->ProjectionX("fCent_1", 4, 4);
  TH1D *fCent_r = fNormalisationHist_r->ProjectionX("fCent_2", 4, 4);
  TH1D *fCent = new TH1D(*fCent_q);
  fCent->Clear();
  fCent->Add(fCent_q,fCent_r,1,1);

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
