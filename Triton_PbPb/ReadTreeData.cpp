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
#include "../Config.h"

using namespace utils;
using namespace triton;

const bool split_qr = false;

void ReadTreeData(const float cutDCAz = 1.f, const int cutTPCcls = 69, const float cutDCAxy = 0.10f, const float cutChi2TPC = 2.5f, const char *outFileName = "TreeOutData", const char *outFileOption = "recreate", const char *flagSelections = "true", const bool data = true)
{
  TFile *outFile = TFile::Open(Form("%s/%s.root", kResDir, outFileName), outFileOption);
  TDirectory *dirOutFile = outFile->mkdir(Form("%1.1f_%d_%1.2f_%1.2f", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC));
  dirOutFile->cd();

  // define dcaxy track selections
  auto trackSelectionsDCAxy = Form("std::abs(dcaxy)<%f",cutDCAxy);
  auto trackSelectionsChi2TPC = Form("(chi2tpc<%f)",cutChi2TPC);

  // define pt bins
  double pTbins[kNPtBins + 1] = {1.4f,1.6f,2.f,2.4f,3.f};

  // define nSigma bins
  double sigmaBins[kNSigmaBins + 1];
  for (int iSigmaBins = 0; iSigmaBins < kNSigmaBins + 1; ++iSigmaBins)
  {
    sigmaBins[iSigmaBins] = kLowestNSigma + iSigmaBins * kDeltaNSigma;
  }

  // read tree
  ROOT::EnableImplicitMT(50); // use all cores
  RDataFrame *dataFrame;
  if (data && split_qr)
    dataFrame = new RDataFrame("RTree", {Form("%s/AnalysisResults_LHC18qr_q.root", kDataDir), Form("%s/AnalysisResults_LHC18qr_r.root", kDataDir)}); // get tree from file
  else if (data && !split_qr)
    dataFrame = new RDataFrame("RTree", Form("%s/AnalysisResults_LHC18qr_reduced_new.root", kDataDir)); // get tree from file
  else
    dataFrame = new RDataFrame("RTree", Form("%s/mc.root", kDataDir));                                                                                                                 // get tree from mc                                                                                                            // get tree from file
  auto trackSelect = dataFrame->Filter(Form("(%s) && (std::abs(dcaz)<%f) && (tpcPIDcls>%d) && %s", kTrackSelectionsEta, cutDCAz, cutTPCcls, trackSelectionsChi2TPC)).Define("pT", "std::abs(pt)").Define("matter","pt>0"); // apply track selections

  dirOutFile->cd();
  // TPC counts
  auto trackSelectCutDCA = trackSelect.Filter(Form("(%s) && (%s) && (std::abs(tpcNsigma)<2.f)", trackSelectionsDCAxy, flagSelections)); // apply track selections && (std::abs(tpcNsigma)<3.f)
  auto antiTritonCutDCA = trackSelectCutDCA.Filter("pt<0");                                                    // select antimatter and change pt sign
  auto TritonCutDCA = trackSelectCutDCA.Filter("pt>0");                                                        // select matter
  auto fATOFcounts = antiTritonCutDCA.Histo3D({"fATOFcounts", "fATOFcounts", kNCentBins, kCentBins, kNPtBins, pTbins, kNSigmaBins, sigmaBins}, "centrality", "pT", "tofNsigma");
  auto fMTOFcounts = TritonCutDCA.Histo3D({"fMTOFcounts", "fMTOFcounts", kNCentBins, kCentBins, kNPtBins, pTbins, kNSigmaBins, sigmaBins}, "centrality", "pT", "tofNsigma");

  // DCA
  auto antiTriton = trackSelect.Filter(Form("pt<0 && (%s) && (std::abs(tpcNsigma)<2.f) && (std::abs(tofNsigma)<5.f)", flagSelections)); // select antimatter and change pt sign
  auto Triton = trackSelect.Filter(Form("pt>0 && (%s) && (std::abs(tpcNsigma)<2.f) && (std::abs(tofNsigma)<5.f)", flagSelections));     // select matter
  auto fADCA = antiTriton.Histo3D({"fADCAxyTOF", "fADCAxyTOF", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");
  auto fMDCA = Triton.Histo3D({"fMDCAxyTOF", "fMDCAxyTOF", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");

  // use snapshots to write trees
  // ROOT::RDF::RSnapshotOptions opts;
  // opts.fLazy = true;
  // using SnapRet_t = ROOT::RDF::RResultPtr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>;
  // std::vector<SnapRet_t> rets;
  // rets.emplace_back(antiTritonCutDCA.Snapshot("fATreeTrackCuts", Form("%s/%s_anti.root", kResDir, outFileName), {"tofNsigma", "pT", "centrality"}, opts));
  // rets.emplace_back(TritonCutDCA.Snapshot("fMTreeTrackCuts", Form("%s/%s_matt.root", kResDir, outFileName), {"tofNsigma", "pT", "centrality"}, opts));
  // *TritonCutDCA.Count();
  trackSelectCutDCA.Filter("pT<3. && pT>1. && std::abs(tofNsigma)<12.f").Snapshot("fTreeTrackCuts", Form("%s/%s_tree.root", kResDir, outFileName), {"tofNsigma", "pT", "centrality", "matter"}); // _%1.1f_%d_%1.2f_%1.2f

  dirOutFile->cd();

  // set axis labels
  fATOFcounts->GetXaxis()->SetTitle(kAxisTitleCent);
  fATOFcounts->GetYaxis()->SetTitle(kAxisTitlePt);
  fATOFcounts->GetZaxis()->SetTitle(kAxisTitleNSigma);
  fMTOFcounts->GetXaxis()->SetTitle(kAxisTitleCent);
  fMTOFcounts->GetYaxis()->SetTitle(kAxisTitlePt);
  fMTOFcounts->GetZaxis()->SetTitle(kAxisTitleNSigma);

  fADCA->GetXaxis()->SetTitle(kAxisTitleCent);
  fADCA->GetYaxis()->SetTitle(kAxisTitlePt);
  fADCA->GetZaxis()->SetTitle(kAxisTitleDCA);
  fMDCA->GetXaxis()->SetTitle(kAxisTitleCent);
  fMDCA->GetYaxis()->SetTitle(kAxisTitlePt);
  fMDCA->GetZaxis()->SetTitle(kAxisTitleDCA);

  // write
  fATOFcounts->Write();
  fMTOFcounts->Write();

  fADCA->Write();
  fMDCA->Write();

  outFile->Close();

  // merge files
  TFile *outFile2 = TFile::Open(Form("%s/%s.root", kResDir, outFileName), "update");
  TDirectory *dirOut = outFile2->mkdir(Form("%1.1f_%d_%1.2f_%1.2f", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC), "", true);
  TFile *outFileTree = TFile::Open(Form("%s/%s_tree.root", kResDir, outFileName));
  TTree *fTree = (TTree *)outFileTree->Get("fTreeTrackCuts");
  dirOut->cd();
  fTree->CloneTree()->Write();
  outFile2->Close();
  outFileTree->Close();
}