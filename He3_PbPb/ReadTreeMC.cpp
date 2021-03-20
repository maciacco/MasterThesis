// ReadTreeMC.cpp
// This macro reads MC tree

#include <string>

#include <Riostream.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>
#include <TLatex.h>

#include "Utils.h"
#include "Config.h"

using namespace utils;

void ReadTreeMC(const float cutDCAz = 1.f, const int cutTPCcls = 89, const char *outFileName = "TreeOutMC", const char *outFileOption = "recreate", const char *flagSelections = "( ( (std::abs(pt)<2.5f) && (trackingPID==7) ) || !(std::abs(pt)<2.5f) )")
{
  TFile *outFile = TFile::Open(Form("%s/%s.root", kResDir, outFileName), outFileOption);
  TDirectory *dirOutFile = outFile->mkdir(Form("%1.1f_%d_", cutDCAz, cutTPCcls));
  dirOutFile->cd();

  // define pt bins
  double pTbins[kNPtBins + 1] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 6.5f, 7.f, 8.f, 10.f};

  // read tree
  ROOT::EnableImplicitMT();                                                  // use all cores
  ROOT::RDataFrame dataFrameR("RTree", Form("%s/mc.root", kDataDir));        // get reconstructed tree from file
  ROOT::RDataFrame dataFrameRSec("RTree", Form("%s/mc_sec.root", kDataDir)); // get reconstructed (deuteron) tree from file
  ROOT::RDataFrame dataFrameS("STree", Form("%s/mc.root", kDataDir));        // get simulated tree from file

  dirOutFile->cd();
  // track selection
  auto trackSelectR = dataFrameR.Filter(Form("(%s) && (std::abs(dcaz)<%f) && (tpcPIDcls>%d)", kTrackSelectionsEta, cutDCAz, cutTPCcls)).Define("pT", "std::abs(2.f*pt)");                          // apply track selections
  auto trackSelectRflag = trackSelectR.Filter(Form("(%s) && (std::abs(tpcNsigma)<3.f)", flagSelections));                                                                                          // apply flag selections and tpcNsigma selection
  auto trackSelectRSec = dataFrameRSec.Filter(Form("(%s) && (std::abs(dcaz)<%f) && (tpcPIDcls>%d) && (std::abs(tpcNsigma)<3.f)", kTrackSelectionsEta, cutDCAz, cutTPCcls)).Define("pT", "2.f*pt"); // apply track selections (deuteron -> scale momentum by 2)
  auto trackSelectRPrimary = trackSelectRflag.Filter("(flag & BIT(1))==2");                                                                                                                        // select primary particles
  auto trackSelectRSecondary = trackSelectRSec.Filter("(flag & BIT(2))==4");                                                                                                                       // select secondary particles
  auto trackSelectRSecondaryWeak = trackSelectRflag.Filter("(flag & BIT(3))==8");                                                                                                                  // select secondary particles from weak decays (hypertriton feed-down)

  // select matter/antimatter
  auto antiHe3R = trackSelectRPrimary.Filter("pt<0"); // select antimatter and change pt sign (primary)
  auto he3R = trackSelectRPrimary.Filter("pt>0");     // select matter (primary)
  auto antiHe3RSecWeak = trackSelectR.Filter(Form("pt<0 && (%s)", flagSelections));
  auto he3RSecWeak = trackSelectR.Filter(Form("pt>0 && (%s)", flagSelections));

  // DCAxy histograms
  auto fADCAPrimary = antiHe3R.Histo3D({"fADCAPrimary", "fADCAPrimary", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");
  auto fMDCAPrimary = he3R.Histo3D({"fMDCAPrimary", "fMDCAPrimary", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");
  auto fMDCASecondary = trackSelectRSecondary.Histo3D({"fMDCASecondary", "fMDCASecondary", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");
  auto fADCASecondaryWeak = antiHe3RSecWeak.Histo3D({"fADCASecondaryWeak", "fADCASecondaryWeak", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");
  auto fMDCASecondaryWeak = he3RSecWeak.Histo3D({"fMDCASecondaryWeak", "fMDCASecondaryWeak", kNCentBins, kCentBins, kNPtBins, pTbins, kNDCABins, kDCABins}, "centrality", "pT", "dcaxy");

  // ITS-TPC (all data track selections applied) - also for wd efficiency
  auto antiHe3RCutDCA = antiHe3R.Filter(kTrackSelectionsDCAxy); // cut on dcaxy
  auto he3RCutDCA = he3R.Filter(kTrackSelectionsDCAxy);
  auto antiHe3RSecWeakCutDCA = antiHe3RSecWeak.Filter(kTrackSelectionsDCAxy);
  auto he3RSecWeakCutDCA = he3RSecWeak.Filter(kTrackSelectionsDCAxy);
  auto fAITS_TPC = antiHe3RCutDCA.Histo2D({"fAITS_TPC", "fAITS_TPC", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fMITS_TPC = he3RCutDCA.Histo2D({"fMITS_TPC", "fMITS_TPC", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fAITS_TPCwd = antiHe3RCutDCA.Histo2D({"fAITS_TPCwd", "fAITS_TPCwd", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fMITS_TPCwd = he3RCutDCA.Histo2D({"fMITS_TPCwd", "fMITS_TPCwd", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");

  // total generated (also for wd efficiency)
  auto antiHe3Swd = dataFrameS.Filter("pdg<0 && ((flag & BIT(2))==4)").Define("pT", "pt");
  auto he3Swd = dataFrameS.Filter("pdg>0 && ((flag & BIT(2))==4)").Define("pT", "pt");
  auto antiHe3S = dataFrameS.Filter("pdg<0 && ((flag & BIT(0))>0)").Define("pT", "pt");
  auto he3S = dataFrameS.Filter("pdg>0 && ((flag & BIT(0))>0)").Define("pT", "pt");

  auto fATotal = antiHe3S.Histo2D({"fATotal", "fATotal", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fMTotal = he3S.Histo2D({"fMTotal", "fMTotal", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fATotalWd = antiHe3Swd.Histo2D({"fATotalWd", "fATotalWd", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fMTotalWd = he3Swd.Histo2D({"fMTotalWd", "fMTotalWd", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");

#ifdef PID
  auto trackSelectRNoPID = trackSelectR.Filter(Form("(%s)", kTrackSelectionsDCAxy)); // no tpcNsigma cut
  auto trackSelectRPID = trackSelectRNoPID.Filter("trackingPID==7");
  auto trackSelectRalphaPID = trackSelectRNoPID.Filter("trackingPID==8");
  auto antiHe3RPID = trackSelectRPID.Filter("pt<0");
  auto he3RPID = trackSelectRPID.Filter("pt>0");
  auto antiHe4RPID = trackSelectRalphaPID.Filter("pt<0");
  auto he4RPID = trackSelectRalphaPID.Filter("pt>0");
  auto antiHe3RnoPID = trackSelectRNoPID.Filter("pt<0");
  auto he3RnoPID = trackSelectRNoPID.Filter("pt>0");
  auto fAPID = antiHe3RPID.Histo2D({"fAPID", "fAPID", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fMPID = he3RPID.Histo2D({"fMPID", "fMPID", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fAAlphaPID = antiHe4RPID.Histo2D({"fAAlphaPID", "fAAlphaPID", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fMAlphaPID = he4RPID.Histo2D({"fMAlphaPID", "fMAlphaPID", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fANoPID = antiHe3RnoPID.Histo2D({"fANoPID", "fANoPID", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
  auto fMNoPID = he3RnoPID.Histo2D({"fMNoPID", "fMNoPID", kNCentBins, kCentBins, kNPtBins, pTbins}, "centrality", "pT");
#endif // PID

  // set axis labels
  fADCAPrimary->GetXaxis()->SetTitle(kAxisTitleCent);
  fADCAPrimary->GetYaxis()->SetTitle(kAxisTitlePt);
  fADCAPrimary->GetZaxis()->SetTitle(kAxisTitleDCA);
  fMDCAPrimary->GetXaxis()->SetTitle(kAxisTitleCent);
  fMDCAPrimary->GetYaxis()->SetTitle(kAxisTitlePt);
  fMDCAPrimary->GetZaxis()->SetTitle(kAxisTitleDCA);
  fMDCASecondary->GetXaxis()->SetTitle(kAxisTitleCent);
  fMDCASecondary->GetYaxis()->SetTitle(kAxisTitlePt);
  fMDCASecondary->GetZaxis()->SetTitle(kAxisTitleDCA);
  fADCASecondaryWeak->GetXaxis()->SetTitle(kAxisTitleCent);
  fADCASecondaryWeak->GetYaxis()->SetTitle(kAxisTitlePt);
  fADCASecondaryWeak->GetZaxis()->SetTitle(kAxisTitleDCA);
  fMDCASecondaryWeak->GetXaxis()->SetTitle(kAxisTitleCent);
  fMDCASecondaryWeak->GetYaxis()->SetTitle(kAxisTitlePt);
  fMDCASecondaryWeak->GetZaxis()->SetTitle(kAxisTitleDCA);

  fAITS_TPC->GetXaxis()->SetTitle(kAxisTitleCent);
  fAITS_TPC->GetYaxis()->SetTitle(kAxisTitlePt);
  fMITS_TPC->GetXaxis()->SetTitle(kAxisTitleCent);
  fMITS_TPC->GetYaxis()->SetTitle(kAxisTitlePt);
  fATotal->GetXaxis()->SetTitle(kAxisTitleCent);
  fATotal->GetYaxis()->SetTitle(kAxisTitlePt);
  fMTotal->GetXaxis()->SetTitle(kAxisTitleCent);
  fMTotal->GetYaxis()->SetTitle(kAxisTitlePt);

#ifdef PID
  fAPID->GetXaxis()->SetTitle(kAxisTitleCent);
  fAPID->GetYaxis()->SetTitle(kAxisTitlePt);
  fAAlphaPID->GetXaxis()->SetTitle(kAxisTitleCent);
  fAAlphaPID->GetYaxis()->SetTitle(kAxisTitlePt);
  fANoPID->GetXaxis()->SetTitle(kAxisTitleCent);
  fANoPID->GetYaxis()->SetTitle(kAxisTitlePt);
  fMPID->GetXaxis()->SetTitle(kAxisTitleCent);
  fMPID->GetYaxis()->SetTitle(kAxisTitlePt);
  fMAlphaPID->GetXaxis()->SetTitle(kAxisTitleCent);
  fMAlphaPID->GetYaxis()->SetTitle(kAxisTitlePt);
  fMNoPID->GetXaxis()->SetTitle(kAxisTitleCent);
  fMNoPID->GetYaxis()->SetTitle(kAxisTitlePt);
#endif // PID

  // wd efficiency
  fAITS_TPCwd->GetXaxis()->SetTitle(kAxisTitleCent);
  fAITS_TPCwd->GetYaxis()->SetTitle(kAxisTitlePt);
  fMITS_TPCwd->GetXaxis()->SetTitle(kAxisTitleCent);
  fMITS_TPCwd->GetYaxis()->SetTitle(kAxisTitlePt);
  fATotalWd->GetXaxis()->SetTitle(kAxisTitleCent);
  fATotalWd->GetYaxis()->SetTitle(kAxisTitlePt);

  // write to file
  fAITS_TPC->Write();
  fMITS_TPC->Write();
  fADCAPrimary->Write();
  fMDCAPrimary->Write();
  fMDCASecondary->Write();
  fADCASecondaryWeak->Write();
  fMDCASecondaryWeak->Write();
  fATotal->Write();
  fMTotal->Write();

#ifdef PID
  fAPID->Write();
  fAAlphaPID->Write();
  fANoPID->Write();
  fMPID->Write();
  fMAlphaPID->Write();
  fMNoPID->Write();
#endif // PID

  fAITS_TPCwd->Write();
  fMITS_TPCwd->Write();
  fATotalWd->Write();
  fMTotalWd->Write();

  outFile->Close();
}