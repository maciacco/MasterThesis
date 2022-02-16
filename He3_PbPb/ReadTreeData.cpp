// ReadTreeData.cpp
// This macro reads data tree and creates output data file

#include <string>

#include <Riostream.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TH1D.h>
#include <TH3F.h>
#include <TString.h>
#include <TLatex.h>
#include <TDirectory.h>
#include <TTree.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using namespace utils;
using namespace he3;

void ReadTreeData(const float cutDCAz = 1.f, const int cutTPCcls = 89, const float cutDCAxy = 0.10f, const char *outFileName = "TreeOutData", const char *outFileOption = "recreate", const char *flagSelections = "( ( (std::abs(pt)<2.5f) && (trackingPID==7) ) || !(std::abs(pt)<2.5f) )", const bool data = true)
{
  TFile *outFile = TFile::Open(Form("%s/%s.root", kResDir, outFileName), outFileOption);
  TDirectory *dirOutFile = outFile->mkdir(Form("%1.1f_%d_%1.2f", cutDCAz, cutTPCcls, cutDCAxy));
  dirOutFile->cd();

  // define dcaxy track selections
  auto trackSelectionsDCAxy = Form("std::abs(dcaxy)<%f",cutDCAxy);

  // define pt bins
  double pTbins[kNPtBins + 1] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 6.5f, 7.f, 8.f, 10.f};

  // define nSigma bins
  double sigmaBins[kNSigmaBins + 1];
  for (int iSigmaBins = 0; iSigmaBins < kNSigmaBins + 1; ++iSigmaBins)
  {
    sigmaBins[iSigmaBins] = kLowestNSigma + iSigmaBins * kDeltaNSigma;
  }

  // read tree
  ROOT::EnableImplicitMT(); // use all cores
  RDataFrame *dataFrame;
  if (data)
    dataFrame = new RDataFrame("RTree", {Form("%s/AnalysisResults_LHC18q.root", kDataDir), Form("%s/AnalysisResults_LHC18r.root", kDataDir)}); // get tree from file
  else
    dataFrame = new RDataFrame("RTree", Form("%s/mc.root", kDataDir));                                                                                                                 // get tree from mc                                                                                                            // get tree from file
  auto trackSelect = dataFrame->Filter(Form("(%s) && (std::abs(dcaz)<%f) && (tpcPIDcls>%d)", kTrackSelectionsEta, cutDCAz, cutTPCcls)).Define("pT", "std::abs(pt*2.f)"); // apply track selections

  dirOutFile->cd();
  // TPC counts
  auto trackSelectCutDCA = trackSelect.Filter(Form("(%s) && (%s)", trackSelectionsDCAxy, flagSelections)); // apply track selections
  auto antiHe3CutDCA = trackSelectCutDCA.Filter("pt<0");                                                    // select antimatter and change pt sign
  auto he3CutDCA = trackSelectCutDCA.Filter("pt>0");                                                        // select matter
  auto fATPCcounts = antiHe3CutDCA.Histo3D({"fATPCcounts", "fATPCcounts", kNCentBins, kCentBins, kNPtBins, pTbins, kNSigmaBins, sigmaBins}, "centrality", "pT", "tpcNsigma");
  auto fMTPCcounts = he3CutDCA.Histo3D({"fMTPCcounts", "fMTPCcounts", kNCentBins, kCentBins, kNPtBins, pTbins, kNSigmaBins, sigmaBins}, "centrality", "pT", "tpcNsigma");

  // DCA
  auto antiHe3 = trackSelect.Filter(Form("pt<0 && (%s) && (std::abs(tpcNsigma)<3.f)", flagSelections)); // select antimatter and change pt sign
  auto he3 = trackSelect.Filter(Form("pt>0 && (%s) && (std::abs(tpcNsigma)<3.f)", flagSelections));     // select matter
  auto fADCA = antiHe3.Histo3D({"fADCAxyTPC", "fADCAxyTPC", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");
  auto fMDCA = he3.Histo3D({"fMDCAxyTPC", "fMDCAxyTPC", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");

  // use snapshots to write trees
  ROOT::RDF::RSnapshotOptions opts;
  opts.fLazy = true;
  using SnapRet_t = ROOT::RDF::RResultPtr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>;
  std::vector<SnapRet_t> rets;
  rets.emplace_back(antiHe3CutDCA.Snapshot("fATreeTrackCuts", Form("%s/%s_anti.root", kResDir, outFileName), {"tpcNsigma", "pT", "centrality"}, opts));
  rets.emplace_back(he3CutDCA.Snapshot("fMTreeTrackCuts", Form("%s/%s_matt.root", kResDir, outFileName), {"tpcNsigma", "pT", "centrality"}, opts));
  *he3CutDCA.Count();

  dirOutFile->cd();

  // set axis labels
  fATPCcounts->GetXaxis()->SetTitle(kAxisTitleCent);
  fATPCcounts->GetYaxis()->SetTitle(kAxisTitlePt);
  fATPCcounts->GetZaxis()->SetTitle(kAxisTitleNSigma);
  fMTPCcounts->GetXaxis()->SetTitle(kAxisTitleCent);
  fMTPCcounts->GetYaxis()->SetTitle(kAxisTitlePt);
  fMTPCcounts->GetZaxis()->SetTitle(kAxisTitleNSigma);

  fADCA->GetXaxis()->SetTitle(kAxisTitleCent);
  fADCA->GetYaxis()->SetTitle(kAxisTitlePt);
  fADCA->GetZaxis()->SetTitle(kAxisTitleDCA);
  fMDCA->GetXaxis()->SetTitle(kAxisTitleCent);
  fMDCA->GetYaxis()->SetTitle(kAxisTitlePt);
  fMDCA->GetZaxis()->SetTitle(kAxisTitleDCA);

  // write
  fATPCcounts->Write();
  fMTPCcounts->Write();

  fADCA->Write();
  fMDCA->Write();

  outFile->Close();

  // merge files
  TFile *outFile2 = TFile::Open(Form("%s/%s.root", kResDir, outFileName), "update");
  TDirectory *dirOut = outFile2->mkdir(Form("%1.1f_%d_%1.2f", cutDCAz, cutTPCcls, cutDCAxy), "", true);
  TFile *outFileAnti = TFile::Open(Form("%s/%s_anti.root", kResDir, outFileName));
  TFile *outFileMatt = TFile::Open(Form("%s/%s_matt.root", kResDir, outFileName));
  TTree *fTreeAnti = (TTree *)outFileAnti->Get("fATreeTrackCuts");
  TTree *fTreeMatt = (TTree *)outFileMatt->Get("fMTreeTrackCuts");
  dirOut->cd();
  fTreeAnti->CloneTree()->Write();
  fTreeMatt->CloneTree()->Write();
  outFile2->Close();
  outFileAnti->Close();
  outFileMatt->Close();
}