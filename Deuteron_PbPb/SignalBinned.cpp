// SignalUnbinned.cpp

#include <stdlib.h>
#include <string>

#include <Riostream.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooBinning.h>
#include <RooFitResult.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"
#include "../utils/RooGausExp.h"
#include "../utils/RooDSCBShape.h"
#include "../utils/RooGausDExp.h"

using namespace utils;
using namespace deuteron;

const double kNSigma = 3; // define interval for bin counting

void SignalBinned(const char *cutSettings = "", const bool binCounting = false, const int bkg_shape = 1, const char *inFileDat = "AnalysisResults", const char *outFileName = "SignalDeuteron", const char *outFileOption = "recreate", const bool extractSignal = true, const bool useDSCB = false, const bool binCountingNoFit = false)
{

  // make signal extraction plots directory
  system(Form("mkdir %s/signal_extraction", kPlotDir));
  system(Form("mkdir %s/signal_extraction_preliminary", kPlotDir));

  // killing RooFit output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel = kError; // Suppressing warning outputs
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *outFile = TFile::Open(TString::Format("%s/%s.root", kOutDir, outFileName), outFileOption); // output file
  if (outFile->GetDirectory(Form("%s_%d_%d", cutSettings, binCounting, bkg_shape)))
    return;
  TDirectory *dirOutFile = outFile->mkdir(Form("%s_%d_%d", cutSettings, binCounting, bkg_shape));
  TFile *dataFile = TFile::Open(TString::Format("%s/%s.root", kDataDir, inFileDat)); // open data TFile
  if (!dataFile)
  {
    std::cout << "File not found!" << std::endl; // check data TFile opening
    return;
  }

  // get TTList
  std::string listName = Form("mpuccio_deuterons_%s", cutSettings);
  TTList *list = (TTList *)dataFile->Get(listName.data());

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  { // loop on antimatter/matter
    // Get histograms from file
    std::string histName = Form("f%sTOFsignal", kAntimatterMatter[iMatt]);
    TH3F *fTOFSignal = (TH3F *)list->Get(histName.data());
    if (!fTOFSignal)
    {
      std::cout << "Hstogram not found!" << std::endl; // check data TFile opening
      return;
    }

    // make plot subdirectory
    system(Form("mkdir %s/signal_extraction/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));
    system(Form("mkdir %s/signal_extraction_preliminary/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));

    dirOutFile->cd();
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    { // loop over centrality classes
      TH1D fTOFrawYield("fRawYield", "fRawYield", kNPtBins, kPtBins);          // declare raw yields histogram
      TH1D fSignificance("fSignificance", "fSignificance", kNPtBins, kPtBins); // declare significance
      TH1D fSigma("fSigma", "fSigma", kNPtBins, kPtBins);
      TH1D fMean("fMean", "fMean", kNPtBins, kPtBins);
      TH1D fAlpha("fAlpha", "fAlpha", kNPtBins, kPtBins);
      int nUsedPtBins = 22;

      for (int iPtBin = 1; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        double ptMin = fTOFrawYield.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fTOFrawYield.GetXaxis()->GetBinUpEdge(iPtBin);
        double centMin = kCentBinsLimitsDeuteron[iCent][0];
        double centMax = kCentBinsLimitsDeuteron[iCent][1];

        int pTbinsIndexMin = fTOFSignal->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fTOFSignal->GetYaxis()->FindBin(ptMax - 0.005);
        int centBinMin = kCentBinsDeuteron[iCent][0];
        int centBinMax = kCentBinsDeuteron[iCent][1];

        // project histogram
        std::cout << "Pt bins: min=" << ptMin << "; max=" << ptMax << "; indexMin=" << pTbinsIndexMin << "; indexMax=" << pTbinsIndexMax << "; centralityBinMin=" << centBinMin << "; centralityBinMax=" << centBinMax << std::endl;
        TH1D *tofSignalProjection = fTOFSignal->ProjectionZ(Form("f%sTOFSignal_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);

        // roofit variables
        RooRealVar tofSignal("tofSignal", "#it{m}^{2}-#it{m}^{2}_{PDG}", kTOFSignalMin, kTOFSignalMax, "GeV^{2}/#it{c^{4}}");

        // roofit histogram
        RooDataHist data("data", "data", RooArgList(tofSignal), tofSignalProjection);

        // gaussian core fit
        tofSignalProjection->SetAxisRange(-0.5, 1.);
        double maxBin = tofSignalProjection->GetBinCenter(tofSignalProjection->GetMaximumBin());
        tofSignalProjection->SetAxisRange(kTOFSignalMin, kTOFSignalMax);
        TF1 gaus("gausCore", "gaus(0)+expo(3)", kTOFSignalMin, kTOFSignalMax);
        gaus.SetParLimits(0, 100., 1.e6);
        gaus.SetParLimits(1, 0.0, 0.9);
        if (ptMin < 1.19)
          gaus.SetParameter(1, 0.8);
        else
          gaus.SetParameter(1, 0.2);
        gaus.SetParLimits(2, 0.1, 0.3);
        gaus.SetParLimits(3, 0., 100.);
        gaus.SetParLimits(4, -1.0, -0.1);
        tofSignalProjection->Fit("gausCore", "QMRL", "", maxBin - 0.8, maxBin + 0.2);

        // build composite model
        RooRealVar mean("#mu", "mean", gaus.GetParameter(1), 0.0, gaus.GetParameter(1) + 0.5, "GeV^{2}/#it{c^{4}}");
        RooRealVar *sigma = new RooRealVar("#sigma", "sigma", gaus.GetParameter(2), gaus.GetParameter(2) - 0.05, gaus.GetParameter(2) + 0.2, "GeV^{2}/#it{c^{4}}");
        RooRealVar *alpha1 = new RooRealVar("#alpha_{L}", "alpha1", -2.0, -0.5);
        RooRealVar *alpha = new RooRealVar("#alpha_{R}", "alpha", 0.0, 10.);
        RooRealVar a1("a1","a1",1.,1.,2.);
        RooRealVar a2("a2","a2",1., 1.,2.);
        RooRealVar n1("n1","n1", 2.,10.);
        RooRealVar n2("n2","n2", 2.,10.);
        if (ptMin < 1.41)
        {
          alpha->setRange(0.5, 2.0);
          alpha->setVal(0.6);
        }
        else
        {
          alpha->setRange(1.0, 3.0);
          alpha->setVal(1.1);
          if (ptMin > 2.79 && ptMin < 2.81)
          { // special case #1
            alpha->setVal(1.5);
          }
          alpha1->setRange(-2.0, -1.0);
          if (ptMin > 2.39 && ptMin < 2.41 && iMatt == 1)
          { // special case #1
            alpha1->setVal(-1.9);
          }
        }
        RooAbsPdf *signal;
        if (iMatt == 0)
        {
          signal = new RooGausExp("signal", "signal", tofSignal, mean, *sigma, *alpha);
        }
        else
        {
          signal = new RooGausDExp("signal", "signal", tofSignal, mean, *sigma, *alpha1, *alpha);
        }
        if (useDSCB && iMatt == 1){
          signal = new RooDSCBShape("signal", "signal", tofSignal, mean, *sigma, a1, n1, a2, n2);
        }
        RooAbsPdf *background1;
        RooAbsPdf *background2;
        RooAbsPdf *background;
        RooRealVar *slope;
        RooRealVar *slope1;
        RooRealVar *slope2;
        RooRealVar *nBackground1;
        RooRealVar *nBackground2;
        if (bkg_shape == 1)
        { // expo
          if ((ptMin < 1.41) || (ptMin > 1.59 && ptMin < 1.61 && iCent == 1))
          {
            slope1 = new RooRealVar("#tau_{mismatch}", "slope1", -0.2, -1., -0.1);
            background = (RooAbsPdf *)new RooExponential("background", "background", tofSignal, *slope1);
          }
          else
          {
            slope1 = new RooRealVar("#tau_{mismatch}", "slope1", -2.1, -5., -2.0);
            slope2 = new RooRealVar("#tau_{p}", "slope2", -0.2, -1.0, -0.1);
            background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
            background2 = (RooAbsPdf *)new RooExponential("background2", "background2", tofSignal, *slope2);
            nBackground1 = new RooRealVar("#it{f}_{Bkg,p}", "nBackground1", 0., 1.);

            background = new RooAddPdf("background", "background", RooArgList(*background1, *background2), RooArgList(*nBackground1));
          }
        }
        else
        { // TODO: sum of expo + straight line
          std::cout << "No background shape with bkg_shape = 0!" << std::endl;
          return;
        }
        RooRealVar nSignal("N_{sig}", "nSignal", 1., 1000000.);
        RooRealVar nBackground("N_{bkg}", "nBackground", 1., 10000000.);
        RooAddPdf *model = new RooAddPdf("model", "model", RooArgList(*signal, *background), RooArgList(nSignal, nBackground));

        if (extractSignal)
        {
          // fit model
          RooFitResult *r;
          for (int i = 0; i < 2; ++i)
          {
            r = model->fitTo(data, RooFit::Save());
          }
          std::cout << "fit status: " << r->status() << ";" << std::endl;
          std::cout << "covariance quality: " << r->covQual() << std::endl;

          // signal range
          double signal_min = mean.getVal() - kNSigma * sigma->getVal(), signal_max = mean.getVal() + kNSigma * sigma->getVal();
          tofSignal.setRange("signalRange", signal_min, signal_max);

          // background integral
          auto components = model->getComponents();
          auto bkgPdf = (RooAbsPdf *)components->find("background");
          double bkgIntegral = (((RooAbsPdf *)model->pdfList().at(1))->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
          double bkgIntegral_val = nBackground.getVal() * bkgIntegral;
          // std::cout << bkgIntegral << std::endl;

          double rawYield, rawYieldError, counts;
          if (binCounting)
          {
            // total counts
            counts = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", signal_min, signal_max));
            rawYield = counts - bkgIntegral_val; // signal=counts-error
            rawYieldError = TMath::Sqrt(counts); // counts=signal+bkg => correct error from poisson statistics
            if (binCountingNoFit)
            {
              rawYield = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", -0.7815 * kNSigma, 0.7815 * kNSigma)); // only for He4 in signal loss studies
            }
          }
          else
          {
            rawYield = nSignal.getVal();
            rawYieldError = nSignal.getError();
          }

          // fill raw yield histogram
          fTOFrawYield.SetBinContent(fTOFrawYield.FindBin((ptMax + ptMin) / 2.), rawYield);
          fTOFrawYield.SetBinError(fTOFrawYield.FindBin((ptMax + ptMin) / 2.), rawYieldError);

          // fill significance histogram
          if (binCounting)
            fSignificance.SetBinContent(fSignificance.FindBin((ptMax + ptMin) / 2.), rawYield / rawYieldError);
          else
            fSignificance.SetBinContent(fSignificance.FindBin((ptMax + ptMin) / 2.), rawYield / TMath::Sqrt(rawYield + bkgIntegral_val));
          fSignificance.SetBinError(fSignificance.FindBin((ptMax + ptMin) / 2.), 0.);

          // fill sigma histogram
          fSigma.SetBinContent(fSigma.FindBin((ptMax + ptMin) / 2.), sigma->getVal());
          fSigma.SetBinError(fSigma.FindBin((ptMax + ptMin) / 2.), sigma->getError());

          // fill mean histogram
          fMean.SetBinContent(fMean.FindBin((ptMax + ptMin) / 2.), mean.getVal());
          fMean.SetBinError(fMean.FindBin((ptMax + ptMin) / 2.), mean.getError());

          // fill alpha histogram
          fAlpha.SetBinContent(fAlpha.FindBin((ptMax + ptMin) / 2.), alpha->getVal());
          fAlpha.SetBinError(fAlpha.FindBin((ptMax + ptMin) / 2.), alpha->getError());
        }

        // frame
        int nBins = 36;
        TString plotTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]);
        RooPlot *xframe = tofSignal.frame(RooFit::Bins(nBins), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigma_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], ptMin, ptMax)));
        data.plotOn(xframe, RooFit::Name("dataNsigma"));

        if (extractSignal)
        {
          model->plotOn(xframe, RooFit::Components("background"), RooFit::Name("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
          model->plotOn(xframe, RooFit::Components("signal"), RooFit::Name("signal"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
          model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue));
          xframe->remove("model",false);
          xframe->remove("background",false);
          xframe->remove("signal",false);
          model->plotOn(xframe, RooFit::Components("background"), RooFit::Name("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
          model->plotOn(xframe, RooFit::Components("signal"), RooFit::Name("signal"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
          model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue));
          double x_min = 0.62;
          if(ptMax < 1.1) x_min = 0.12;
          model->paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.4f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(x_min, 0.99, 0.88));
          xframe->getAttText()->SetTextSize(0.03);
          xframe->getAttLine()->SetLineWidth(0);
          xframe->getAttFill()->SetFillColor(kWhite);
          xframe->getAttFill()->SetFillStyle(4000);
          xframe->GetYaxis()->SetMaxDigits(2);
        }

        tofSignalProjection->Write();
        xframe->Write();

        // save to png
        TCanvas canv;
        canv.SetName(plotTitle);
        xframe->Draw("");
        // canv.SetLogy();
        canv.Print(Form("%s/signal_extraction/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], ptMin, ptMax));

        tofSignalProjection->SetMarkerStyle(20);
        tofSignalProjection->SetMarkerSize(0.8);
        tofSignalProjection->Draw("pe");
        canv.Print(Form("%s/signal_extraction_preliminary/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], ptMin, ptMax));

      } // end of loop on pt bins

      // set raw yield histogram style and write to file
      fTOFrawYield.SetTitle(TString::Format("%s raw yield, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fTOFrawYield.SetName(TString::Format("f%sTOFrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fTOFrawYield.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fTOFrawYield.GetYaxis()->SetTitle("#it{N_{raw}}");
      fTOFrawYield.SetMarkerStyle(20);
      fTOFrawYield.SetMarkerSize(0.8);
      fTOFrawYield.SetOption("PE");
      fTOFrawYield.Write();

      // set significance histogram style and write to file
      fSignificance.SetTitle(TString::Format("%s significance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fSignificance.SetName(TString::Format("f%sSignificance_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fSignificance.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSignificance.GetYaxis()->SetTitle("S/#sqrt{N}");
      fSignificance.SetMarkerStyle(20);
      fSignificance.SetMarkerSize(0.8);
      fSignificance.SetOption("PE");
      fSignificance.Write();

      // set sigma histogram style and write to file
      fSigma.SetTitle(TString::Format("%s sigma, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fSigma.SetName(TString::Format("f%sSigma_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fSigma.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSigma.GetYaxis()->SetTitle("#sigma (GeV^2/#it{c}^{4})");
      fSigma.SetMarkerStyle(20);
      fSigma.SetMarkerSize(0.8);
      fSigma.SetOption("PE");
      fSigma.Write();

      // set mean histogram style and write to file
      fMean.SetTitle(TString::Format("%s mean, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fMean.SetName(TString::Format("f%sMean_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fMean.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fMean.GetYaxis()->SetTitle("#mu (GeV^2/#it{c}^{4})");
      fMean.SetMarkerStyle(20);
      fMean.SetMarkerSize(0.8);
      fMean.SetOption("PE");
      fMean.Write();

      // set sigma histogram style and write to file
      fAlpha.SetTitle(TString::Format("%s alpha, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fAlpha.SetName(TString::Format("f%sAlpha_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fAlpha.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fAlpha.GetYaxis()->SetTitle("#alpha");
      fAlpha.SetMarkerStyle(20);
      fAlpha.SetMarkerSize(0.8);
      fAlpha.SetOption("PE");
      fAlpha.Write();
    } // end of loop on centrality bin
  }   // end of loop on antimatter/matter
  outFile->Close();
} // end of macro Signal.cpp
