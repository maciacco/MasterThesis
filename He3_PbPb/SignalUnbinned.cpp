// SignalUnbinned.cpp
// This macro extracts raw yields via an unbinned analysis + cut over trackingPID flag
// - apply cut over trackingPID (= 7 for He3)
// - define signal region after an unbinned fit to nsigma ditribution (gaussian+expo)
// - extract yields by bin counting in the signal region (background subtraction)

#include <string>

#include <Riostream.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooBinning.h>
#include <RooFitResult.h>

#include "Utils.h"
#include "Config.h"

using namespace utils;

const double kNSigma = 3; // define interval for bin counting

void SignalUnbinned(const float cutDCAz = 1.f, const int cutTPCcls = 89, const bool binCounting = true, const int bkg_shape = 1, const char *inFileDat = "TreeOutData", const char *outFileName = "SignalHe3", const char *outFileOption = "recreate", const bool extractSignal = true, const bool binCountingNoFit = false)
{
  // killing RooFit output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel = kError; // Suppressing warning outputs
  gStyle->SetOptStat(0);

  TTree *inTree;

  TFile *outFile = TFile::Open(TString::Format("%s/%s.root", kOutDir, outFileName), outFileOption); // output file
  if (outFile->GetDirectory(Form("%1.1f_%d_%d_%d", cutDCAz, cutTPCcls, binCounting, bkg_shape)))
    return;
  TDirectory *dirOutFile = outFile->mkdir(Form("%1.1f_%d_%d_%d", cutDCAz, cutTPCcls, binCounting, bkg_shape));
  TFile *dataFile = TFile::Open(TString::Format("%s/%s.root", kResDir, inFileDat)); // open data TFile
  if (!dataFile)
  {
    std::cout << "File not found!" << std::endl; // check data TFile opening
    return;
  }

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  { // loop on antimatter/matter
    // Get trees from file
    std::string treeName = Form("%1.1f_%d_/f", cutDCAz, cutTPCcls);
    treeName.append(kAntimatterMatter[iMatt]);
    treeName += "TreeTrackCuts";
    inTree = (TTree *)dataFile->Get(treeName.c_str());
    if (!inTree)
    {
      std::cout << "Tree not found!" << std::endl; // check data TFile opening
      return;
    }

    dirOutFile->cd();
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    { // loop on centrality classes
      RooRealVar tpcNsigma("tpcNsigma", "n#sigma", kNSigmaMin, kNSigmaMax, "a.u.");
      RooRealVar pT("pT", "p_{T}", 1.f, 10.f, "GeV/#it{c}");
      RooRealVar centrality("centrality", "Centrality", 0.f, 90.f, "%");
      RooDataSet data("data", "data", RooArgSet(tpcNsigma, pT, centrality), RooFit::Import(*inTree));

      // define pt bins
      double pTbins[kNPtBins + 1] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 6.5f, 7.f, 8.f, 10.f};

      TH1D fTPCrawYield("fRawYield", "fRawYield", kNPtBins, pTbins);          // declare TPC raw yields
      TH1D fSignificance("fSignificance", "fSignificance", kNPtBins, pTbins); // declare significance
      TH1D fSigma("fSigma", "fSigma", kNPtBins, pTbins);
      TH1D fMean("fMean", "fMean", kNPtBins, pTbins);
      int nUsedPtBins = 12;
      if (iCent == kNCentClasses - 1)
      {
        int nPtBins = 13;
        double pTbinsNew[] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 7.f, 8.f, 10.f};
        fTPCrawYield.SetBins(nPtBins, pTbinsNew);
        fSignificance.SetBins(nPtBins, pTbinsNew);
        fSigma.SetBins(nPtBins, pTbinsNew);
        nUsedPtBins = 9;
      }
      for (int iPtBin = 3; iPtBin < nUsedPtBins + 3; ++iPtBin)
      { // loop on pT bins
        double minPt = fTPCrawYield.GetXaxis()->GetBinLowEdge(iPtBin);
        double maxPt = fTPCrawYield.GetXaxis()->GetBinUpEdge(iPtBin);

        // apply cuts to dataset
        RooDataSet dataCut = *(RooDataSet *)data.reduce(RooFit::Cut(Form("pT>%f && pT<%f && centrality>%.0f && centrality<%.0f", minPt, maxPt, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1])));
        RooDataSet dataPt = *(RooDataSet *)dataCut.reduce(RooArgSet(tpcNsigma));

        // build composite model
        RooRealVar mean("#mu", "mean", 0., -1., 1.);
        RooRealVar *sigma = new RooRealVar("#sigma", "sigma", 0.8, 0.6, 0.9);
        RooGaussian signal("signal", "signal", tpcNsigma, mean, *sigma);
        RooAbsPdf *background;
        RooRealVar *slope;
        if (bkg_shape == 1)
        { // expo
          slope = new RooRealVar("slope", "slope", -0.4, 0.);
          background = (RooAbsPdf *)new RooExponential("background", "background", tpcNsigma, *slope);
        }
        else
        { // pol1
          slope = new RooRealVar("slope", "slope", -0.3, 0.);
          background = (RooAbsPdf *)new RooPolynomial("background", "background", tpcNsigma, RooArgList(*slope));
        }
        RooRealVar nSignal("N_{sig}", "nSignal", 1., 10000.);
        RooRealVar nBackground("N_{bkg}", "nBackground", 0., 1000.);
        RooAddPdf model("model", "model", RooArgList(signal, *background), RooArgList(nSignal, nBackground));

        if (extractSignal)
        {
          // fit model
          RooFitResult *r = model.fitTo(dataPt, RooFit::Save());

          if (r->status() > 0)
            continue;
          // signal range
          double signal_min = mean.getVal() - kNSigma * sigma->getVal(), signal_max = mean.getVal() + kNSigma * sigma->getVal();
          tpcNsigma.setRange("signalRange", signal_min, signal_max);

          // background integral
          auto components = model.getComponents();
          auto bkgPdf = (RooAbsPdf *)components->find("background");
          RooRealVar *bkgIntegral = (RooRealVar *)bkgPdf->createIntegral(RooArgSet(tpcNsigma), RooArgSet(tpcNsigma), "signalRange");
          double bkgIntegral_val = nBackground.getVal() * bkgIntegral->getVal();
          // std::cout<<bkgIntegral_val<<std::endl;

          double rawYield, rawYieldError, counts;
          if (binCounting)
          {
            // total counts
            counts = dataPt.sumEntries(Form("tpcNsigma>%f && tpcNsigma<%f", signal_min, signal_max));
            rawYield = counts - bkgIntegral_val; // signal=counts-error
            rawYieldError = TMath::Sqrt(counts); // counts=signal+bkg => correct error from poisson statistics
            if (binCountingNoFit)
            {
              rawYield = dataPt.sumEntries(Form("tpcNsigma>%f && tpcNsigma<%f", -0.7815 * kNSigma, 0.7815 * kNSigma)); // only for He4 in signal loss studies
            }
          }
          else
          {
            rawYield = nSignal.getVal();
            rawYieldError = nSignal.getError();
          }

          // fill raw yield histogram
          fTPCrawYield.SetBinContent(fTPCrawYield.FindBin((maxPt + minPt) / 2.), rawYield);
          fTPCrawYield.SetBinError(fTPCrawYield.FindBin((maxPt + minPt) / 2.), rawYieldError);

          // fill significance histogram
          if (binCounting)
            fSignificance.SetBinContent(fSignificance.FindBin((maxPt + minPt) / 2.), rawYield / rawYieldError);
          else
            fSignificance.SetBinContent(fSignificance.FindBin((maxPt + minPt) / 2.), rawYield / TMath::Sqrt(rawYield + bkgIntegral_val));
          fSignificance.SetBinError(fSignificance.FindBin((maxPt + minPt) / 2.), 0.);

          // fill sigma histogram
          fSigma.SetBinContent(fSigma.FindBin((maxPt + minPt) / 2.), sigma->getVal());
          fSigma.SetBinError(fSigma.FindBin((maxPt + minPt) / 2.), sigma->getError());

          // fill sigma histogram
          fMean.SetBinContent(fMean.FindBin((maxPt + minPt) / 2.), mean.getVal());
          fMean.SetBinError(fMean.FindBin((maxPt + minPt) / 2.), mean.getError());
        }

        // frame
        TString plotTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", minPt, maxPt, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]);
        RooPlot *xframe = tpcNsigma.frame(RooFit::Bins(50), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigma_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1], minPt, maxPt)));
        dataPt.plotOn(xframe, RooFit::Name("dataNsigma"));

        if (extractSignal)
        {
          model.plotOn(xframe, RooFit::Components("background"), RooFit::Name("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
          model.plotOn(xframe, RooFit::Components("signal"), RooFit::Name("signal"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
          model.plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue));
          model.paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.4f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(0.68, 0.96, 0.96));
        }

        xframe->Write();

        // save to pdf
        TCanvas canv;
        canv.SetName(plotTitle);
        xframe->Draw("");
        if (iMatt == 0 && iCent == 0 && iPtBin == 3)
        {
          canv.Print(Form("%s/%s.pdf[", kPlotDir, outFileName));
          canv.Print(Form("%s/%s.pdf", kPlotDir, outFileName));
        }
        else if (iMatt == 1 && iCent == (kNCentClasses-1) && iPtBin == (nUsedPtBins+2))
        {
          canv.Print(Form("%s/%s.pdf", kPlotDir, outFileName));
          canv.Print(Form("%s/%s.pdf]", kPlotDir, outFileName));
        }
        else
          canv.Print(Form("%s/%s.pdf", kPlotDir, outFileName));

      } // end of loop on pt bins

      // set raw yield histogram style and write to file
      fTPCrawYield.SetTitle(TString::Format("%s raw yield, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fTPCrawYield.SetName(TString::Format("f%sTPCrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fTPCrawYield.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fTPCrawYield.GetYaxis()->SetTitle("#it{N_{raw}}");
      fTPCrawYield.SetMarkerStyle(20);
      fTPCrawYield.SetMarkerSize(0.8);
      fTPCrawYield.SetOption("PE");
      fTPCrawYield.Write();

      // set significance histogram style and write to file
      fSignificance.SetTitle(TString::Format("%s significance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSignificance.SetName(TString::Format("f%sSignificance_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSignificance.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSignificance.GetYaxis()->SetTitle("S/#sqrt{N}");
      fSignificance.SetMarkerStyle(20);
      fSignificance.SetMarkerSize(0.8);
      fSignificance.SetOption("PE");
      fSignificance.Write();

      // set sigma histogram style and write to file
      fSigma.SetTitle(TString::Format("%s sigma, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSigma.SetName(TString::Format("f%sSigma_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fSigma.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSigma.GetYaxis()->SetTitle("#sigma (a.u.)");
      fSigma.SetMarkerStyle(20);
      fSigma.SetMarkerSize(0.8);
      fSigma.SetOption("PE");
      fSigma.Write();

      // set sigma histogram style and write to file
      fMean.SetTitle(TString::Format("%s mean, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fMean.SetName(TString::Format("f%sMean_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fMean.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fMean.GetYaxis()->SetTitle("#mu (a.u.)");
      fMean.SetMarkerStyle(20);
      fMean.SetMarkerSize(0.8);
      fMean.SetOption("PE");
      fMean.Write();
    } // end of loop on centrality bin
    inTree->Reset();
  } // end of loop on antimatter/matter
  outFile->Close();
} // end of macro Signal.cpp