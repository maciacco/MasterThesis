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

using namespace utils;
using namespace deuteron;

const double kNSigma = 3; // define interval for bin counting

void SignalBinned(const char *cutSettings = "", const bool binCounting = false, const int bkg_shape = 1, const char *inFileDat = "AnalysisResults", const char *outFileName = "SignalDeuteron", const char *outFileOption = "recreate", const bool extractSignal = true, const bool binCountingNoFit = false)
{
  // make signal extraction plots directory
  system(Form("mkdir %s/signal_extraction", kPlotDir));

  // killing RooFit output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel = kError; // Suppressing warning outputs
  gStyle->SetOptStat(0);

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
  std::string listName = Form("mpuccio_deuterons_%s",cutSettings);
  TTList *list = (TTList *)dataFile->Get(listName.data());

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  { // loop on antimatter/matter
    // Get histograms from file
    std::string histName = Form("f%sTOFsignal",kAntimatterMatter[iMatt]);
    TH3F *fTOFSignal = (TH3F *)list->Get(histName.data());
    if (!fTOFSignal)
    {
      std::cout << "Hstogram not found!" << std::endl; // check data TFile opening
      return;
    }

    // make plot subdirectory
    system(Form("mkdir %s/signal_extraction/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));

    dirOutFile->cd();
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    { // loop over centrality classes
      double pTbins[25] = {0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.4f, 1.6f, 1.8f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.2f, 4.4f, 5.0f, 6.0f, 8.0f};

      TH1D fTOFrawYield("fRawYield", "fRawYield", kNPtBins, pTbins);          // declare raw yields histogram
      TH1D fSignificance("fSignificance", "fSignificance", kNPtBins, pTbins); // declare significance
      TH1D fSigma("fSigma", "fSigma", kNPtBins, pTbins);
      TH1D fMean("fMean", "fMean", kNPtBins, pTbins);
      int nUsedPtBins = 22;

      for (int iPtBin = 1; iPtBin < nUsedPtBins+1; ++iPtBin)
      { // loop on pT bins
        double ptMin = fTOFrawYield.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fTOFrawYield.GetXaxis()->GetBinUpEdge(iPtBin);
        double centMin = kCentBinsLimitsDeuteron[iCent][0];
        double centMax = kCentBinsLimitsDeuteron[iCent][1];

        int pTbinsIndexMin = fTOFSignal->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fTOFSignal->GetYaxis()->FindBin(ptMax-0.005);

        // project histogram
        std::cout<<"Pt bins index: min="<<pTbinsIndexMin<<"; max="<<pTbinsIndexMax<<std::endl;
        TH1D *tofSignalProjection = fTOFSignal->ProjectionZ(Form("f%sTOFSignal_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsLimitsDeuteron[iCent][0],  kCentBinsLimitsDeuteron[iCent][1], pTbinsIndexMin, pTbinsIndexMax);

        // roofit variables
        RooRealVar tofSignal("tofSignal", "#it{m}^{2}-#it{m}^{2}_{PDG}", kTOFSignalMin, kTOFSignalMax, "GeV^{2}/#it{c^{4}}");

        // roofit histogram
        RooDataHist data("data", "data", RooArgList(tofSignal), tofSignalProjection);

        // gaussian core fit
        tofSignalProjection->SetAxisRange(-0.5,1.);
        double maxBin = tofSignalProjection->GetBinCenter(tofSignalProjection->GetMaximumBin());
        TF1 gaus("gausCore","gaus",maxBin-0.2f,maxBin+0.2f);
        tofSignalProjection->Fit("gausCore","RL","",maxBin-0.2f,maxBin+0.2f);
        tofSignalProjection->SetAxisRange(-2.,3.);

        // build composite model
        RooRealVar mean("#mu", "mean", gaus.GetParameter(1), gaus.GetParameter(1)-0.4, gaus.GetParameter(1)+0.4, "GeV^{2}/#it{c^{4}}");
        if(ptMin>0.85f){
          mean.setConstant(kTRUE);
        }

        RooRealVar *sigma = new RooRealVar("#sigma", "sigma", gaus.GetParameter(2), 0.1f, gaus.GetParameter(2)+0.2f, "GeV^{2}/#it{c^{4}}");
        if(ptMin>0.85f){
          sigma->setMin(0.09f);
          sigma->setMax(gaus.GetParameter(2)+0.1f);
        }

        RooRealVar *tau = new RooRealVar("#tau", "tau", 0.6f, 2.0f);
        if(ptMin>1.15f){
          tau->setMin(0.9f);
          tau->setVal(1.1f);
        }
        RooGausExp signal("signal", "signal", tofSignal, mean, *sigma, *tau);
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
          TF1 expoLeft("expoLeft","expo",-2.f,-1.4f);
          tofSignalProjection->Fit("expoLeft");
          slope1 = new RooRealVar("slope_{1}", "slope1", expoLeft.GetParameter(1), -7.f, -1.f);

          slope2 = new RooRealVar("slope_{2}", "slope2", -0.2f, -1.f, -0.1f);

          background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
          background2 = (RooAbsPdf *)new RooExponential("background2", "background2", tofSignal, *slope2);
          nBackground1 = new RooRealVar("nBackground1", "nBackground1", 0.f, 1.f);
          if(ptMin<1.65f){nBackground1->setVal(0.f); nBackground1->setConstant(kTRUE); slope1->setConstant(kTRUE);}
          background = new RooAddPdf("background","background",RooArgList(*background1,*background2),RooArgList(*nBackground1));
        }
        else
        {
          std::cout << "No background shape with bkg_shape = 0!" << std::endl;
          return;
        }
        RooRealVar nSignal("N_{sig}", "nSignal", 1., 10000000.);
        RooRealVar nBackground("N_{bkg}", "nBackground", 10., 10000000.);
        RooAddPdf model("model", "model", RooArgList(signal, *background), RooArgList(nSignal, nBackground));

        if (extractSignal)
        {
          // fit model
          RooFitResult *r = model.fitTo(data, RooFit::Save());

          // signal range
          double signal_min = mean.getVal() - kNSigma * sigma->getVal(), signal_max = mean.getVal() + kNSigma * sigma->getVal();
          tofSignal.setRange("signalRange", signal_min, signal_max);

          // background integral
          auto components = model.getComponents();
          auto bkgPdf = (RooAbsPdf *)components->find("background");
          double bkgIntegral = (((RooAbsPdf *)model.pdfList().at(1))->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
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

          // fill sigma histogram
          fMean.SetBinContent(fMean.FindBin((ptMax + ptMin) / 2.), mean.getVal());
          fMean.SetBinError(fMean.FindBin((ptMax + ptMin) / 2.), mean.getError());
        }

        // frame
        int nBins = 36;
        TString plotTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]);
        RooPlot *xframe = tofSignal.frame(RooFit::Bins(nBins), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigma_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], ptMin, ptMax)));
        data.plotOn(xframe, RooFit::Name("dataNsigma"));

        if (extractSignal)
        {
          model.plotOn(xframe, RooFit::Components("background"), RooFit::Name("background"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
          model.plotOn(xframe, RooFit::Components("signal"), RooFit::Name("signal"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
          model.plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue));
          model.paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.4f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(0.68, 0.96, 0.96));
        }

        tofSignalProjection->Write();
        xframe->Write();

        // save to png
        TCanvas canv;
        canv.SetName(plotTitle);
        xframe->Draw("");
        canv.Print(Form("%s/signal_extraction/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], ptMin, ptMax));

      } // end of loop on pt bins

      // set raw yield histogram style and write to file
      fTOFrawYield.SetTitle(TString::Format("%s raw yield, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fTOFrawYield.SetName(TString::Format("f%sTPCrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
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

      // set sigma histogram style and write to file
      fMean.SetTitle(TString::Format("%s mean, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fMean.SetName(TString::Format("f%sMean_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]));
      fMean.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSigma.GetYaxis()->SetTitle("#mu (GeV^2/#it{c}^{4})");
      fMean.SetMarkerStyle(20);
      fMean.SetMarkerSize(0.8);
      fMean.SetOption("PE");
      fMean.Write();
    } // end of loop on centrality bin
  }   // end of loop on antimatter/matter
  outFile->Close();
} // end of macro Signal.cpp
