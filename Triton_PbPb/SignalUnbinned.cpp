// SignalUnbinned.cpp
// This macro extracts raw yields via an unbinned analysis + cut over trackingPID flag
// - apply cut over trackingPID (= 7 for He3)
// - define signal region after an unbinned fit to nsigma ditribution (gaussian+expo)
// - extract yields by bin counting in the signal region (background subtraction)

#include <stdlib.h>
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

#include "../utils/Utils.h"
#include "../utils/Config.h"

using namespace utils;
using namespace triton;

const char *sign[] = {"<", ">"};

const double kNSigma = 3; // define interval for bin counting
const char *amString[] = {"anti", "matt"};

void SignalUnbinned(const float cutDCAz = 1.f, const int cutTPCcls = 69, const float cutDCAxy = 0.10f, const float cutChi2TPC = 2.5f, const bool binCounting = true, const int bkg_shape = 1, const char *inFileDat = "TreeOutData", const char *outFileName = "SignalHe3", const char *outFileOption = "recreate", const bool extractSignal = true, const bool binCountingNoFit = false)
{
  // make signal extraction plots directory
  system(Form("mkdir %s/signal_extraction", kPlotDir));

  // killing RooFit output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel = kError; // Suppressing warning outputs
  gStyle->SetOptStat(0);

  TTree *inTree;

  TFile *outFile = TFile::Open(TString::Format("%s/%s.root", kOutDir, outFileName), outFileOption); // output file
  if (outFile->GetDirectory(Form("%1.1f_%d_%1.2f_%1.2f_%d_%d", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, binCounting, bkg_shape)))
    return;
  TDirectory *dirOutFile = outFile->mkdir(Form("%1.1f_%d_%1.2f_%1.2f_%d_%d", cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, binCounting, bkg_shape));

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  { // loop on antimatter/matter

    TFile *dataFile = TFile::Open(TString::Format("%s/%s.root", kResDir, inFileDat)); // open data TFile
    if (!dataFile)
    {
      std::cout << "File not found!" << std::endl; // check data TFile opening
      return;
    }

    // Get trees from file
    inTree = nullptr;
    std::string treeName = Form("%1.1f_%d_%1.2f_%1.2f/fTreeTrackCuts",cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC);
    inTree = (TTree *)dataFile->Get(treeName.c_str());
    if (!inTree)
    {
      std::cout << "Tree not found!" << std::endl; // check data TFile opening
      return;
    }

    // make plot subdirectory
    system(Form("mkdir %s/signal_extraction/%s_%1.1f_%d_%1.2f_%1.2f_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, binCounting, bkg_shape));

    dirOutFile->cd();
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    { // loop on centrality classes
      RooRealVar tofNsigma("tofNsigma", "n#sigma", kNSigmaMin, kNSigmaMax, "a.u.");
      RooRealVar pT("pT", "p_{T}", 1.f, 10.f, "GeV/#it{c}");
      RooRealVar centrality("centrality", "Centrality", 0.f, 90.f, "%");
      RooRealVar matter("matter", "matter", 0.f, 2.f);
      RooDataSet data("data", "data", RooArgSet(tofNsigma, pT, centrality, matter), RooFit::Import(*inTree));

      // define pt bins
      double pTbins[kNPtBins + 1] = {1.4f,1.6f,2.f,2.4f,3.f};

      TH1D fTOFrawYield("fRawYield", "fRawYield", kNPtBins, pTbins);          // declare TPC raw yields
      TH1D fSignificance("fSignificance", "fSignificance", kNPtBins, pTbins); // declare significance
      TH1D fSigma("fSigma", "fSigma", kNPtBins, pTbins);
      TH1D fMean("fMean", "fMean", kNPtBins, pTbins);
      int nUsedPtBins = 4;

      for (int iPtBin = 1; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        double minPt = fTOFrawYield.GetXaxis()->GetBinLowEdge(iPtBin);
        double maxPt = fTOFrawYield.GetXaxis()->GetBinUpEdge(iPtBin);

        // apply cuts to dataset
        RooDataSet dataCut = *(RooDataSet *)data.reduce(RooFit::Cut(Form("pT>%f && pT<%f && centrality>%.0f && centrality<%.0f && matter %s 0.5", minPt, maxPt, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1], sign[iMatt])));
        RooDataSet dataPt = *(RooDataSet *)dataCut.reduce(RooArgSet(tofNsigma));

        // data anti + matt
        RooDataSet dataCut_all = *(RooDataSet *)data.reduce(RooFit::Cut(Form("pT>%f && pT<%f && centrality>%.0f && centrality<%.0f", minPt, maxPt, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1])));
        RooDataSet dataPt_all = *(RooDataSet *)dataCut_all.reduce(RooArgSet(tofNsigma));

        // build composite model
        RooRealVar mean("#mu", "mean", 0., -1., 4.);
        RooRealVar *sigma = new RooRealVar("#sigma", "sigma", 0.8, 0.6, 2.);
        RooGaussian signal("signal", "signal", tofNsigma, mean, *sigma);
        RooAbsPdf *background;
        RooRealVar *slope;
        RooRealVar *a;
        RooRealVar *b;
        if (bkg_shape == 1)
        { // expo
          slope = new RooRealVar("#tau", "tau", -1., 1.);
          background = (RooAbsPdf *)new RooExponential("background", "background", tofNsigma, *slope);
        }
        else if (bkg_shape == 0)
        { // pol1
          slope = new RooRealVar("slope", "slope", -6, 6.);
          background = (RooAbsPdf *)new RooPolynomial("background", "background", tofNsigma, RooArgList(*slope));
        }
        else
        { // pol2
          a = new RooRealVar("a", "a", 0.03, 0.02, 0.04);
          b = new RooRealVar("b", "b", -0.2, -0.28, 0.1);
          background = (RooAbsPdf *)new RooPolynomial("background", "background", tofNsigma, RooArgList(*b, *a));
        }
        RooRealVar nSignal("N_{sig}", "nSignal", 1., 100000.);
        RooRealVar nBackground("N_{bkg}", "nBackground", 0., 100000.);
        RooAddPdf model("model", "model", RooArgList(signal, *background), RooArgList(nSignal, nBackground));

        if (extractSignal)
        {
          // fit model
          for (int i = 0; i < 2; ++i) model.fitTo(dataPt_all);
          mean.setConstant();
          sigma->setConstant();
          RooFitResult *r = nullptr;
          for (int i = 0; i < 2; ++i)
          {
            r = model.fitTo(dataPt, RooFit::Save());
          }
          mean.setConstant(false);
          sigma->setConstant(false);

          if (r->status() > 0)
            continue;
          // signal range
          double signal_min = mean.getVal() - kNSigma * sigma->getVal(), signal_max = mean.getVal() + kNSigma * sigma->getVal();
          tofNsigma.setRange("signalRange", signal_min, signal_max);

          // background integral
          auto components = model.getComponents();
          auto bkgPdf = (RooAbsPdf *)components->find("background");
          double bkgIntegral = (((RooAbsPdf *)model.pdfList().at(1))->createIntegral(RooArgSet(tofNsigma), RooFit::NormSet(RooArgSet(tofNsigma)), RooFit::Range("signalRange")))->getVal();
          double bkgIntegral_val = nBackground.getVal() * bkgIntegral;
          // std::cout << bkgIntegral << std::endl;

          double rawYield, rawYieldError, counts;
          if (binCounting)
          {
            // total counts
            counts = dataPt.sumEntries(Form("tofNsigma>%f && tofNsigma<%f", signal_min, signal_max));
            rawYield = counts - bkgIntegral_val; // signal=counts-error
            rawYieldError = TMath::Sqrt(counts); // counts=signal+bkg => correct error from poisson statistics
            if (binCountingNoFit)
            {
              rawYield = dataPt.sumEntries(Form("tofNsigma>%f && tofNsigma<%f", -0.7815 * kNSigma, 0.7815 * kNSigma)); // only for He4 in signal loss studies
            }
          }
          else
          {
            rawYield = nSignal.getVal();
            rawYieldError = nSignal.getError();
          }

          // fill raw yield histogram
          fTOFrawYield.SetBinContent(fTOFrawYield.FindBin((maxPt + minPt) / 2.), rawYield);
          fTOFrawYield.SetBinError(fTOFrawYield.FindBin((maxPt + minPt) / 2.), rawYieldError);

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
        int nBins = 40;
        TString plotTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", minPt, maxPt, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]);
        RooPlot *xframe = tofNsigma.frame(RooFit::Bins(nBins), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigma_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1], minPt, maxPt)));
        dataPt.plotOn(xframe, RooFit::Name("dataNsigma"));

        if (extractSignal)
        {
          model.plotOn(xframe, RooFit::Components("background"), RooFit::Name("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
          model.plotOn(xframe, RooFit::Components("signal"), RooFit::Name("signal"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
          model.plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue));
          model.paramOn(xframe, RooFit::Parameters(RooArgSet(mean,*sigma,*slope,nSignal,nBackground)), RooFit::Label(TString::Format("#chi^{2}/NDF = %2.4f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(0.70, 0.92, 0.86));
          xframe->getAttLine()->SetLineWidth(0);
          xframe->getAttText()->SetTextFont(44);
          xframe->getAttText()->SetTextSize(20);
        }

        xframe->Write();

        // save to png
        TCanvas canv;
        canv.SetName(plotTitle);
        xframe->Draw("");
        canv.Print(Form("%s/signal_extraction/%s_%1.1f_%d_%1.2f_%1.2f_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls, cutDCAxy, cutChi2TPC, binCounting, bkg_shape, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1], minPt, maxPt));

      } // end of loop on pt bins

      // set raw yield histogram style and write to file
      fTOFrawYield.SetTitle(TString::Format("%s raw yield, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fTOFrawYield.SetName(TString::Format("f%sTOFrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      fTOFrawYield.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fTOFrawYield.GetYaxis()->SetTitle("#it{N_{raw}}");
      fTOFrawYield.SetMarkerStyle(20);
      fTOFrawYield.SetMarkerSize(0.8);
      fTOFrawYield.SetOption("PE");
      fTOFrawYield.Write();

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
    dataFile->Close();
  } // end of loop on antimatter/matter
  outFile->Close();
} // end of macro Signal.cpp
