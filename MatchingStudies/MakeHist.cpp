#include "Config.h"

void MakeHist (const char* inFileName = "AnalysisResults_q_K0s.root", const char* outFileName = "out_q.root", const char* treeName = "K0sTree", const char* inFileName2 = "")
{
  TFile outFile(Form("%s/%s", kOutDir, outFileName), "recreate");
  std::vector<std::string> inFiles = {inFileName};
  if (!split_polarities)
  {
    inFiles.push_back(inFileName2);
  }
  ROOT::RDataFrame df(treeName, inFiles);
  auto df_filter = df.Filter(Form("cosPA > %f && dcaV0PV < %f && dcaV0tracks < %f && radius > %f && %s", kCutCosPA, kCutDcaV0PV, kCutDcaV0tracks, kCutRadius, kReconstructed));
  if (kDebug)
  {
    auto h = df_filter.Histo1D({"h", "h", 1000, kCutMass[0], kCutMass[1]}, "mass");
    h->Write();
  }
  for (int iC{0}; iC < 2; ++iC)
  {
    auto df_filter_daugh = df_filter.Filter(Form("dca%sPV < %f && std::abs(tpcNsigma%s) < %f", kChargeFlag[iC][0], kCutDCADaughPV, kChargeFlag[iC][0], kCutTPCNsigmaDaugh));
    auto hTPC = df_filter_daugh.Histo2D({Form("f%sTPC", kChargeFlag[iC][0]), "; Centrality (%); #it{p}_{T} (GeV/#it{c})", kNCentBins, kCentBinLimits, kNPtBins, kPtBinLimits}, "centrality", Form("pt%s", kChargeFlag[iC][0]));
    auto hTPCTOF = df_filter_daugh.Filter(Form("hasTOF%s", kChargeFlag[iC][1])).Histo2D({Form("f%sTPCTOF", kChargeFlag[iC][0]), "; Centrality (%); #it{p}_{T} (GeV/#it{c})", kNCentBins, kCentBinLimits, kNPtBins, kPtBinLimits}, "centrality", Form("pt%s", kChargeFlag[iC][0]));
    hTPC->Write();
    hTPCTOF->Write();
  }
  outFile.Close();
}