// SignalUnbinned.cpp

#include <stdlib.h>
#include <string>

#include <Riostream.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TString.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPad.h>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooArgList.h>
#include <RooAbsPdf.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooBinning.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooCBShape.h>
#include <RooHist.h>
#include <RooFitResult.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"
#include "../utils/RooGausExp.h"
#include "../utils/RooDSCBShape.h"
#include "../utils/RooGausDExp.h"

const double extend_roi = 0;
const bool fit_peak = false;

const char* start_time[] = {"T0FillNsigma","NoT0FillNsigma","nSigma"};

using namespace utils;
using namespace proton;

const double kNSigma = 3; // define interval for bin counting

const char* sidebands_fit[] = {"leftSideband,rightSideband","rightSideband"};

void SignalBinned(const char *cutSettings = "", const double roi_nsigma = 8., const bool binCounting = false, const int bkg_shape = 1, const char *inFileDat = "AnalysisResults", const char *outFileName = "SignalProton", const char *outFileOption = "recreate", const bool extractSignal = true, const bool binCountingNoFit = false)
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

  int iNsigma = 0;
  if (roi_nsigma > 7.8 && roi_nsigma < 8.1) iNsigma = 1;
  else if (roi_nsigma > 8.1) iNsigma = 2; 
  TFile *outFile = TFile::Open(TString::Format("%s/%s.root", kOutDir, outFileName), outFileOption); // output file
  if (outFile->GetDirectory(Form("%s_%d_%d_%d", cutSettings, binCounting, bkg_shape,iNsigma)))
    return;
  TDirectory *dirOutFile = outFile->mkdir(Form("%s_%d_%d_%d", cutSettings, binCounting, bkg_shape, iNsigma));
  TFile *dataFile = TFile::Open(TString::Format("%s/%s_kINT7_10_30_50_90.root", kDataDir, inFileDat)); // open data TFile _LHC18qr_lowPtProton_ITSrecalibrated.root _largeNsigma_cutDCAxyChi2TPC

  if (!dataFile)
  {
    if (kVerbose) std::cout << "File not found!" << std::endl; // check data TFile opening
    return;
  }

  // get TTList
  std::string listName = Form("nuclei_proton_%s", cutSettings);
  TTList *list = (TTList *)dataFile->Get(listName.data());

  // merge antimatter + matter histograms
  std::string histNameTPCA = Form("f%sTPCcounts", kAntimatterMatter[0]);
  std::string histNameTPCM = Form("f%sTPCcounts", kAntimatterMatter[1]);
  std::string histNameA = Form("f%sTOF%s", kAntimatterMatter[0],start_time[1]);
  std::string histNameM = Form("f%sTOF%s", kAntimatterMatter[1],start_time[1]);
  TH3F *fTPCSignalA = (TH3F *)list->Get(histNameTPCA.data());
  TH3F *fTOFSignalA = (TH3F *)list->Get(histNameA.data());
  TH3F *fTPCSignalAll_ = (TH3F *)fTPCSignalA->Clone(fTPCSignalA->GetName());
  TH3F *fTOFSignalAll = (TH3F *)fTOFSignalA->Clone(fTOFSignalA->GetName());
  TH3F *fTPCSignalM = (TH3F *)list->Get(histNameTPCM.data());
  TH3F *fTOFSignalM = (TH3F *)list->Get(histNameM.data());
  fTPCSignalAll_->Add(fTPCSignalM);
  fTOFSignalAll->Add(fTOFSignalM);
  TH3F *fTPCSignalAll[2];
  fTPCSignalAll[0]=new TH3F(*fTPCSignalAll_);
  fTPCSignalAll[1]=new TH3F(*fTPCSignalAll_);

  /////////////////////////////////////////////////////////////////////////////////////
  // FIT PROTON PEAK - SAVE PARAMETERS
  /////////////////////////////////////////////////////////////////////////////////////
  double fitParameterMean[kNCentClasses][kNPtBins];
  double fitParameterSigma[kNCentClasses][kNPtBins];
  double fitParameterAlphaL[kNCentClasses][kNPtBins];
  double fitParameterAlphaR[kNCentClasses][kNPtBins];
  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  // FIT PROTON PEAK - SAVE PARAMETERS - TPC
  /////////////////////////////////////////////////////////////////////////////////////
  double fitParameterMeanTPC[kNCentClasses][kNPtBins];
  double fitParameterSigmaTPC[kNCentClasses][kNPtBins];
  double fitParameterAlphaLTPC[kNCentClasses][kNPtBins];
  double fitParameterAlphaRTPC[kNCentClasses][kNPtBins];
  double fitParameterMeanTPCBkg[kNCentClasses][kNPtBins];
  double fitParameterSigmaTPCBkg[kNCentClasses][kNPtBins];
  double fitParameterTauTPCBkg[kNCentClasses][kNPtBins];
  double fitParameterTauTPCBkgExp[kNCentClasses][kNPtBins];
  double fitNSigTPC[kNCentClasses][kNPtBins];
  double fitNBkgTPC[kNCentClasses][kNPtBins];
  /////////////////////////////////////////////////////////////////////////////////////

  for (int iMatt = 1; iMatt > -1; --iMatt)
  { // loop on antimatter/matter
    // Get histograms from file
    std::string histNameTPC = Form("f%sTPCcounts", kAntimatterMatter[iMatt]);
    std::string histName = Form("f%sTOF%s", kAntimatterMatter[iMatt],start_time[1]);
    std::string histNameT0Fill = Form("f%sTOF%s", kAntimatterMatter[iMatt],start_time[0]);
    TH3F *fTOFSignal1 = (TH3F *)list->Get(histName.data());
    if (!fTOFSignal1)
    {
      if (kVerbose) std::cout << "Hstogram not found!" << std::endl; // check data TFile opening
      return;
    }

    TH3F *fTPCSignal = (TH3F *)list->Get(histNameTPC.data());
    TH3F *fTOFSignal = (TH3F *)fTOFSignal1->Clone(fTOFSignal1->GetName());
    TH3F *fTOFSignalT0Fill = (TH3F *)list->Get(histNameT0Fill.data());

    // make plot subdirectory
    system(Form("mkdir %s/signal_extraction/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));
    system(Form("mkdir %s/signal_extraction_preliminary/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));

    dirOutFile->cd();
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {                                                                          // loop over centrality classes
      double roi_nsigma_up = roi_nsigma;
      double roi_nsigma_down = roi_nsigma;
      TH1D fTOFrawYield("fRawYield", "fRawYield", kNPtBins, kPtBins);          // declare raw yields histogram
      TH1D fTPCrawYield("fTPCRawYield", "fTPCRawYield", kNPtBins, kPtBins);          // declare raw yields histogram
      fTPCrawYield.GetXaxis()->SetRangeUser(0.3,1.);
      TH1D fSignificance("fSignificance", "fSignificance", kNPtBins, kPtBins); // declare significance
      TH1D fSigma("fSigma", "fSigma", kNPtBins, kPtBins);
      TH1D fMean("fMean", "fMean", kNPtBins, kPtBins);
      TH1D fAlphaL("fAlphaL", "fAlphaL", kNPtBins, kPtBins);
      TH1D fAlphaR("fAlphaR", "fAlphaR", kNPtBins, kPtBins);
      int nUsedPtBins = 42; // up to 3.00 GeV/c

      for (int iPtBin = 5; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        double ptMin = fTOFrawYield.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fTOFrawYield.GetXaxis()->GetBinUpEdge(iPtBin);
        double centMin = kCentBinsLimitsProton[iCent][0];
        double centMax = kCentBinsLimitsProton[iCent][1];

        int pTbinsIndexMin = fTOFSignal->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fTOFSignal->GetYaxis()->FindBin(ptMax - 0.005);
        int centBinMin = fTOFSignal->GetXaxis()->FindBin(kCentBinsLimitsProton[iCent][0]);      //kCentBinsProton[iCent][0];
        int centBinMax = fTOFSignal->GetXaxis()->FindBin(kCentBinsLimitsProton[iCent][1] - 1.); //kCentBinsProton[iCent][1];

        // project histogram
        if (kVerbose) std::cout << "Pt bins: min=" << ptMin << "; max=" << ptMax << "; indexMin=" << pTbinsIndexMin << "; indexMax=" << pTbinsIndexMax << "; centralityBinMin=" << centBinMin << "; centralityBinMax=" << centBinMax << std::endl;
        TH1D *tpcSignalProjection = fTPCSignal->ProjectionZ(Form("f%sTPCSignal_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTPCSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTPCSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        TH1D *tofSignalProjection = fTOFSignal->ProjectionZ(Form("f%sTOFSignal_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        tofSignalProjection->Rebin(2);
        TH1D *tofSignalProjectionT0Fill = fTOFSignalT0Fill->ProjectionZ(Form("f%sTOFSignalT0Fill_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        tofSignalProjectionT0Fill->Rebin(2);
        tofSignalProjectionT0Fill->Write();

        // project histogram (antimatter + matter)
        TH1D *tpcSignalProjectionAll = fTPCSignalAll[iMatt]->ProjectionZ(Form("f%sTPCSignalAll_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTPCSignalAll[iMatt]->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTPCSignalAll[iMatt]->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        TH1D *tofSignalProjectionAll = fTOFSignalAll->ProjectionZ(Form("f%sTOFSignalAll_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignalAll->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignalAll->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        tofSignalProjectionAll->Rebin(2);

        // limits
        double maxNsigma=20;
        double nSigmaLeft = -20;
        double nSigmaRight = -15.;
        if (ptMin > 0.79)
        {
          if (ptMin > 1.81)
          {
            nSigmaLeft = -17.;
            nSigmaRight = -12.;
          };
          if (ptMin > 1.99)
          {
            nSigmaLeft = -15.;
            nSigmaRight = -10.;
          };
          if (ptMin > 2.18)
          {
            nSigmaLeft = -13.;
            nSigmaRight = -8.;
          };
          if (ptMin > 2.4)
          {
            nSigmaLeft = -10.;
            nSigmaRight = -5.;
          };
          tofSignalProjectionAll->GetXaxis()->SetRangeUser(nSigmaLeft, nSigmaRight);
          tofSignalProjection->GetXaxis()->SetRangeUser(nSigmaLeft, nSigmaRight);
          double maximum = tofSignalProjectionAll->GetBinCenter(tofSignalProjectionAll->GetMaximumBin());
          nSigmaLeft = maximum + 2.5;
          if (ptMin > 2.1) nSigmaLeft = maximum + 2.;
          if (ptMin > 2.49) nSigmaLeft = maximum + 1.;
          nSigmaRight = nSigmaLeft + 3.;
          if (kVerbose) std::cout << "nSigmaLeft = " << nSigmaLeft << std::endl;
        }
        else {
          //continue;
          nSigmaLeft = -25;
          maxNsigma=25;
        }
        if (iCent==4){
          if (ptMin<2.29)maxNsigma=40;
          else maxNsigma=20;
        }
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(nSigmaLeft, kTOFnSigmaMax);
        tofSignalProjection->GetXaxis()->SetRangeUser(nSigmaLeft, kTOFnSigmaMax);

        // DEFINE SIGNAL REGION
        TF1 signalRegionFit("signalRegionFit", "gaus", -20., 20.);
        signalRegionFit.SetParLimits(0, 1., 1.e7);
        signalRegionFit.SetParLimits(1, 0.2, 0.4);
        signalRegionFit.SetParLimits(2, 1.2, 1.6);
        signalRegionFit.SetLineColor(kBlue);
        if(ptMin>0.99){
          tofSignalProjectionAll->GetXaxis()->SetRangeUser(-.5, .5);
          double maximum_signal = tofSignalProjectionAll->GetBinCenter(tofSignalProjectionAll->GetMaximumBin());
          tofSignalProjectionAll->GetXaxis()->SetRangeUser(-20., 20.);
          tofSignalProjectionAll->Fit("signalRegionFit", "QRL+", "", maximum_signal - 1., maximum_signal + 1.);
        }
        else if(ptMin<0.99){
          tofSignalProjection->GetXaxis()->SetRangeUser(-.5, .5);
          double maximum_signal = tofSignalProjection->GetBinCenter(tofSignalProjection->GetMaximumBin());
          tofSignalProjection->GetXaxis()->SetRangeUser(-20., 20.);
          tofSignalProjection->Fit("signalRegionFit", "QRL+", "", maximum_signal - .7, maximum_signal + .7);
        }
        double mean_tmp = signalRegionFit.GetParameter(1);
        double rms_tmp = signalRegionFit.GetParameter(2);

        // roofit data
        RooRealVar tofSignal("tofSignal", "n#sigma_{p}", nSigmaLeft, maxNsigma, "a.u.");

        // roofit histogram
        RooDataHist data("data", "data", RooArgList(tofSignal), tofSignalProjection);
        RooDataHist dataT0Fill("dataT0Fill", "dataT0Fill", RooArgList(tofSignal), tofSignalProjectionT0Fill);

        // roofit histogram
        RooDataHist dataAll("dataAll", "dataAll", RooArgList(tofSignal), tofSignalProjectionAll);
        if (kVerbose) std::cout << "Number of entries (roofit) = " << dataAll.sumEntries() << std::endl;
        if (kVerbose) std::cout << "Number of entries (root) = " << tofSignalProjectionAll->GetEntries() << std::endl;

        // build composite model
        RooRealVar mean("#mu", "mean", -1., 1., "a.u.");
        RooRealVar *sigma = new RooRealVar("#sigma", "sigma", 1.0, 1.4, "a.u.");
        RooRealVar *alphaL = new RooRealVar("#alpha_{L}", "alphaL", -1.5, -0.8);
        RooRealVar *alphaR = new RooRealVar("#alpha_{R}", "alphaR", 0.8, 1.5);

        RooAbsPdf *signal = new RooGausDExp("signal", "signal", tofSignal, mean, *sigma, *alphaL, *alphaR);

        RooAddPdf *background0;
        RooAbsPdf *background1;
        RooAbsPdf *background2;
        RooRealVar *slope1;
        RooRealVar *slope2;
        RooRealVar *nBackground1;
        RooRealVar *nBackground2;
        RooAddPdf *modelPeak;
        RooFitResult *r;

        RooAddPdf *model;
        RooRealVar nSignal("N_{sig}", "nSignal", 1., 1.e8);

        slope1 = new RooRealVar("#tau_{1}", "slope1", -10., 10.);
        slope2 = new RooRealVar("#tau_{2}", "slope2",-0.9, -10.,0.);
        nBackground1 = new RooRealVar("#it{N}_{Bkg,1}", "nBackground1", 0., 1.e7);
        
        if (bkg_shape == 1)
        { // expo

          if ((ptMin < 1.39 && ptMin > 0.74 && iCent<4)||(ptMin<1.39 && ptMin > 0.74&&iCent==4))
          {
            if (ptMin<0.99)slope1 = new RooRealVar("#tau_{1}", "slope1", -10., 0.);
            background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
            nBackground2 = new RooRealVar("#it{N}_{Bkg,2}", "nBackground2", 0., 0., 1.e9);
            nBackground2->setConstant();
            model = new RooAddPdf("model", "model", RooArgList(*background1), RooArgList(*nBackground1));
          }
          else
          {
            background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
            background0 = new RooAddPdf("background0","background0",RooArgList(*background1),RooArgList(*nBackground1));
            background2 = (RooAbsPdf *)new RooExponential("background2", "background2", tofSignal, *slope2);
            nBackground2 = new RooRealVar("#it{N}_{Bkg,2}", "nBackground2", 1., 1.e9);
            model = new RooAddPdf("model", "model", RooArgList(*background1, *background2), RooArgList(*nBackground1, *nBackground2));
          }
        }
        else
          if (kVerbose) std::cout << "!!!!!" << std::endl; // TODO: UPDATE FOLLOWING BLOCK5);

        int covQ = -999;

        if (extractSignal)
        {
          roi_nsigma_down = roi_nsigma;
          roi_nsigma_up = roi_nsigma;
          if (ptMin < 1.) {roi_nsigma_down=roi_nsigma+3; // default = 9sigma
            roi_nsigma_up=roi_nsigma+6; // 12 sigma
          }
          
          if (ptMin > 1.51) roi_nsigma_down=roi_nsigma-2; // default = 6sigma
          if (ptMin > 1.99) {
            roi_nsigma_down=roi_nsigma-3; // default = 5sigma
            roi_nsigma_up=roi_nsigma+1.;
          }
          if (ptMin > 2.49) {
            roi_nsigma_down=roi_nsigma-4.5; // default = 4sigma
            roi_nsigma_up=roi_nsigma+2.;
          }
          if (ptMin > 2.69) {
            roi_nsigma_down=roi_nsigma-4.5; // default = 3sigma
            roi_nsigma_up=roi_nsigma+3.;
          }
          tofSignal.setRange("leftSideband", nSigmaLeft, mean_tmp - roi_nsigma_down * rms_tmp);
          if (iCent==4&&ptMin<2.49){
            tofSignal.setRange("leftSideband", nSigmaLeft, mean_tmp - (roi_nsigma_down+1.) * rms_tmp);
          }
          tofSignal.setRange("rightSideband", mean_tmp + roi_nsigma_up * rms_tmp, maxNsigma);
          if (iCent==4) {
            if (ptMin<2.29)tofSignal.setRange("rightSideband", 20., maxNsigma);
            else tofSignal.setRange("rightSideband", 14., maxNsigma);
          }
          // fit TOF signal distribution

          if (ptMin>1.99){
            for(int I=0;I<2;++I)background1->fitTo(dataAll, RooFit::Range("rightSideband"));
            slope1->setConstant();
            for (int I=0;I<2;++I)background0->fitTo(data, RooFit::Range("rightSideband"));
            nBackground1->setConstant();
          }
          else {
            for (int I=0;I<2;++I)model->fitTo(dataAll, RooFit::Range(sidebands_fit[(int)((iCent/4)&&((ptMin<1.39 && ptMin > 0.74)))]));
          }
          model->fitTo(data, RooFit::Save(), RooFit::Range(sidebands_fit[(int)((iCent/4)&&((ptMin<1.39 && ptMin > 0.74)))]));
          r = model->fitTo(data, RooFit::Save(), RooFit::Range(sidebands_fit[(int)((iCent/4)&&((ptMin<1.39 && ptMin > 0.74)))]));

          if (kVerbose) std::cout << "fit status: " << r->status() << ";" << std::endl;
          if (kVerbose) std::cout << "covariance quality: " << r->covQual() << std::endl;
          covQ = r->covQual();
          
          // signal range
          if (ptMin<0.99){
            roi_nsigma_down=roi_nsigma;
            roi_nsigma_up=roi_nsigma;
          }
          tofSignal.setRange("signalRange", mean_tmp - (roi_nsigma_down+extend_roi) * rms_tmp, mean_tmp + (roi_nsigma_up+extend_roi) * rms_tmp);
          tofSignal.setRange("aFitRange", mean_tmp - (roi_nsigma_down+extend_roi) * rms_tmp, mean_tmp + (roi_nsigma_up+extend_roi) * rms_tmp);
        }

        // frame
        int nBins = tofSignalProjection->GetNbinsX();
        TString plotTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]);
        RooPlot *xframe = tofSignal.frame(RooFit::Bins(nBins), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigma_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax)));
        data.plotOn(xframe, RooFit::Name("dataNsigma"));
        tofSignal.setRange("full", kTOFnSigmaMin, kTOFnSigmaMax);
        tofSignal.setRange("model", nSigmaLeft, maxNsigma);
        if (extractSignal)
        {
          model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange(sidebands_fit[(int)((iCent/4)&&((ptMin<1.39 && ptMin > 0.74)))]), RooFit::Range(sidebands_fit[(int)((iCent/4)&&((ptMin<1.39 && ptMin > 0.74)))]));
          model->paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(0.58812,0.911028,0.861955));
          xframe->getAttLine()->SetLineWidth(0);
          xframe->getAttText()->SetTextFont(44);
          xframe->getAttText()->SetTextSize(16);
          xframe->remove("model",false);
          model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange(sidebands_fit[(int)((iCent/4)&&((ptMin<1.39 && ptMin > 0.74)))]), RooFit::Range("model"));
          xframe->GetYaxis()->SetMaxDigits(2);

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // FIT PROTON PEAK
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          if (fit_peak){
            nBackground1->setConstant(true);
            nBackground2->setConstant(true);
            slope1->setConstant(true);
            slope2->setConstant(true);
            if (iMatt == 0){
              mean.setVal(fitParameterMean[iCent][iPtBin]);
              mean.setConstant(true);
              sigma->setVal(fitParameterSigma[iCent][iPtBin]);
              sigma->setConstant(true);
              alphaL->setVal(fitParameterAlphaL[iCent][iPtBin]);
              alphaL->setConstant(true);
              alphaR->setVal(fitParameterAlphaR[iCent][iPtBin]);
              alphaR->setConstant(true);
            }
            modelPeak = new RooAddPdf("model2", "model2", RooArgList(*background1, *background2, *signal), RooArgList(*nBackground1, *nBackground2, nSignal));
            if (iMatt == 1){
              modelPeak->fitTo(data, RooFit::Range("signalRange"));
              modelPeak->fitTo(data, RooFit::Range("signalRange"));
            }
            else{
              modelPeak->fitTo(data, RooFit::Range("aFitRange"));
              modelPeak->fitTo(data, RooFit::Range("aFitRange"));
            }
            if (iMatt == 1){
              fitParameterMean[iCent][iPtBin]=mean.getVal();
              fitParameterSigma[iCent][iPtBin]=sigma->getVal();
              fitParameterAlphaL[iCent][iPtBin]=alphaL->getVal();
              fitParameterAlphaR[iCent][iPtBin]=alphaR->getVal();
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (iMatt ==1)
              modelPeak->plotOn(xframe, RooFit::Name("modelPeak"), RooFit::LineColor(kRed), RooFit::NormRange("signalRange"), RooFit::Range("signalRange"));
            else
              modelPeak->plotOn(xframe, RooFit::Name("modelPeak"), RooFit::LineColor(kRed), RooFit::NormRange("aFitRange"), RooFit::Range("signalRange"));
          }

          // background integral
          std::cout << "ptMin = " << ptMin << std::endl;
          double bkgIntegral_1 = (((RooAbsPdf *)(model->pdfList().at(0)))->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
          double bkgIntegral_2 = 0;
          if (((!(ptMin < 1.39 && ptMin > 0.74)) && iCent<4)!= (((ptMin>1.39 || ptMin < 0.74)&&iCent==4))) bkgIntegral_2 = (((RooAbsPdf *)(model->pdfList().at(1)))->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
          double bkgIntegral_val = nBackground1->getVal() *bkgIntegral_1 + nBackground2->getVal() * bkgIntegral_2;

          double rawYield, rawYieldError, counts, countsT0Fill;
          if (binCounting)
          {
            // total counts
            counts = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", mean_tmp - (roi_nsigma_down+extend_roi) * rms_tmp, mean_tmp + (roi_nsigma_up+extend_roi) * rms_tmp));
            countsT0Fill = 0.;
            if (ptMin<2.&&iCent==4)countsT0Fill=dataT0Fill.sumEntries(Form("tofSignal>%f && tofSignal<%f", mean_tmp - (4.) * rms_tmp, mean_tmp + (4.) * rms_tmp));
            else if (ptMin<2.2&&iCent==4)countsT0Fill=dataT0Fill.sumEntries(Form("tofSignal>%f && tofSignal<%f", mean_tmp - (2.) * rms_tmp, mean_tmp + (4.) * rms_tmp));
            rawYield = counts - bkgIntegral_val + countsT0Fill; // signal=counts-error
            if (kVerbose) std::cout << "Counts: " << rawYield << ", Background: "<< bkgIntegral_val << std::endl;
            rawYieldError = TMath::Sqrt(counts + countsT0Fill); // counts=signal+bkg => correct error from poisson statistics
            if (binCountingNoFit)
            {
              rawYield = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", -0.7815 * kNSigma, 0.7815 * kNSigma)); // only for He4 in signal loss studies
            }
            if (r->status()!=0 || r->covQual()<2.5 || ptMin<0.49) rawYield=0;
          }
          else
          {
            rawYield = nSignal.getVal();
            rawYieldError = nSignal.getError();
          }

          // tpcSignalProjection->SetMarkerStyle(20);
          // tpcSignalProjection->SetMarkerSize(1.);
          // tpcSignalProjection->SetDrawOption("pe");
          // tpcSignalProjection->SetMarkerColor(kBlack);
          // tpcSignalProjection->GetXaxis()->SetLabelSize(0.04);
          // tpcSignalProjection->GetXaxis()->SetTitleSize(0.04);
          // tpcSignalProjection->GetXaxis()->SetTitle("n#sigma_{p}");
          // tpcSignalProjection->GetYaxis()->SetTitle(Form("Entries / %.2f",tpcSignalProjection->GetBinWidth(1)));


          // tpcSignalProjectionAll->GetXaxis()->SetRangeUser(-.3,.3);
          // double max_bkg_sig = tpcSignalProjectionAll->GetBinCenter(tpcSignalProjectionAll->GetMaximumBin());
          // tpcSignalProjectionAll->Fit("gaus","RLM","",max_bkg_sig-.5,max_bkg_sig+.5);
          // double mean_sig_all = tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1);
          // double sigma_sig_all = tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2);
          // std::cout<<"mu = "<<mean_sig_all<<"; sigma = "<<sigma_sig_all<<std::endl;
          // tpcSignalProjectionAll->GetXaxis()->SetRangeUser(-15.,-3.);
          // double max_bkg = tpcSignalProjectionAll->GetBinCenter(tpcSignalProjectionAll->GetMaximumBin());
          // tpcSignalProjectionAll->Fit("gaus","RLM","",max_bkg-.8,max_bkg+.8);
          // tpcSignalProjectionAll->GetXaxis()->SetRangeUser(-24.,24.);
          // tpcSignalProjectionAll->Write(Form("%s_All",tpcSignalProjection->GetName()));

          // tpcSignalProjectionAll->Rebin(2);
          // tpcSignalProjection->Rebin(2);

          // double left_fit_range_factor = +1.;
          // if (ptMin>0.59) left_fit_range_factor = -.5;
          // double left_range = -999;
          // for (int iB=1;iB<tpcSignalProjectionAll->GetNbinsX();++iB){
          //   if (tpcSignalProjectionAll->GetBinCenter(iB)>tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2)){
          //     left_range=tpcSignalProjectionAll->GetXaxis()->GetBinLowEdge(iB);
          //     break;
          //   }
          // }
          // RooRealVar nSigmaTPC("nSigmaTPC","n#sigma_{TPC,p}",left_range,6,"a.u.");
          // RooDataHist hNSigmaTPC("hNSigmaTPC","hNSigmaTPC",RooArgList(nSigmaTPC),tpcSignalProjection);
          // RooDataHist hNSigmaTPCAll("hNSigmaTPCAll","hNSigmaTPCAll",RooArgList(nSigmaTPC),tpcSignalProjectionAll);
          // RooRealVar meanTPC("#mu","#mu",-.3,.3);
          // RooRealVar sigmaTPC("#sigma","#sigma",0.8,1.5);
          // RooRealVar alpha1TPC("#alpha_{1}","#alpha_{1}",-3.,-.8);
          // //if (ptMin>0.89)alpha1TPC.setMax(-0.8);
          // if (ptMin < 0.56) alpha1TPC.setRange(-6.,-8.);
          // RooRealVar alpha2TPC("#alpha_{2}","#alpha_{2}",0.8,3.);
          // RooRealVar meanTPCBkg("#mu_{Bkg}","#mu_{Bkg}",tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)-2.,tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+2.);
          // RooRealVar sigmaTPCBkg("#sigma_{Bkg}","#sigma_{Bkg}",0.7,1.8);
          // RooRealVar tauTPCBkg1("#tau_{Bkg,1}","#tau_{Bkg,1}",-4.,-.1);
          // RooRealVar tauTPCBkg("#tau_{Bkg}","#tau_{Bkg}",0.8,2.5);
          // RooGausDExp signalTPC("signalTPC","signalTPC",nSigmaTPC,meanTPC,sigmaTPC,alpha1TPC,alpha2TPC);
          // RooGausExp bkgTPC("bkgTPC","bkgTPC",nSigmaTPC,meanTPCBkg,sigmaTPCBkg,tauTPCBkg);
          // RooRealVar nSignalTPC("n_{Signal}","n_{Signal}",0.,1.e10);
          // RooRealVar nBkgTPC("n_{Bkg}","n_{Bkg}",0.,1.e12);
          // RooAddPdf modelTPC("modelTPC","modelTPC",RooArgList(signalTPC,bkgTPC),RooArgList(nSignalTPC,nBkgTPC));
          // //if (ptMin > 0.44) left_fit_range_factor = -.5;
          // double lowX = left_range;
          // nSigmaTPC.setRange("fitRange",left_range,3.);
          // if (ptMin<0.54){
          //   nSigmaTPC.setRange("fitRange",-6.,3.);
          //   lowX=-6;
          // }
          // for (int i=0;i<2;++i)modelTPC.fitTo(hNSigmaTPCAll,RooFit::Range("fitRange"));
          // /* meanTPC.setConstant();
          // sigmaTPC.setConstant();
          // alpha2TPC.setConstant();
          // alpha1TPC.setConstant(); */
          // // tauTPCBkg.setConstant();
          // modelTPC.fitTo(hNSigmaTPC,RooFit::Range("fitRange"));
          // auto *r=modelTPC.fitTo(hNSigmaTPC,RooFit::Range("fitRange"),RooFit::Save());
          // RooPlot *frameTPC = nSigmaTPC.frame(RooFit::Name(tpcSignalProjection->GetName()),RooFit::Title(Form("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1])));
          // hNSigmaTPC.plotOn(frameTPC,RooFit::Name("nsigmaTPC"));
          // double lowTPCrange=3;
          // double tpc_range_limit = (roi_nsigma-8.)*1;
          // double tpc_tmp_mean=0., tpc_tmp_rms=0.;
          // nSigmaTPC.setRange("signalRange",meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal());
          // modelTPC.plotOn(frameTPC,RooFit::Components("bkgTPC"),RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
          // modelTPC.plotOn(frameTPC,RooFit::Components("signalTPC"),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
          // modelTPC.plotOn(frameTPC,RooFit::Name("model"));
          // modelTPC.paramOn(frameTPC,RooFit::Parameters(RooArgSet(meanTPC,sigmaTPC)),RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", frameTPC->chiSquare("model", "nsigmaTPC"))), RooFit::Layout(0.58812,0.811028,0.561955));
          // //frameTPC->GetXaxis()->SetRangeUser(tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2),5);
          // frameTPC->GetXaxis()->SetRangeUser(lowX,5);
          // frameTPC->getAttText()->SetTextFont(44);
          // frameTPC->getAttText()->SetTextSize(18);
          // frameTPC->getAttLine()->SetLineWidth(0);
          // frameTPC->Write();

          // double background=0.;
          // if (ptMin > 0.49){
          //   background = ((RooAbsPdf *)bkgTPC.createIntegral(RooArgSet(nSigmaTPC), RooFit::NormSet(RooArgSet(nSigmaTPC)), RooFit::Range("signalRange")))->getVal();
          //   background *= nBkgTPC.getVal();
          //   std::cout<< "background = " << background << std::endl;
          // }
          // TCanvas c(Form("%s",tpcSignalProjection->GetName()),Form("%s",tpcSignalProjection->GetName()));
          // c.cd();
          // TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.4, 1.0, 1.0, 0);
          // pad1->SetFillColor(0);
          // pad1->SetFrameBorderMode(0);
          // pad1->SetBottomMargin(0.06);
          // pad1->SetLogy();
          // pad1->Draw();
          // TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.4, 0);
          // pad2->SetFillColor(0);
          // pad2->SetFrameBorderMode(0);
          // pad2->SetFillColor(0);
          // pad2->SetFrameBorderMode(0);
          // pad2->SetTopMargin(0.10);
          // pad2->SetBottomMargin(0.25);
          // pad2->Draw();
          // pad1->cd();
          // tpcSignalProjection->GetXaxis()->SetRangeUser(lowX,5);
          // double rawYieldTPC = hNSigmaTPC.sumEntries(Form("nSigmaTPC>%f && nSigmaTPC<%f",meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal()));
          // if(ptMin>0.49)rawYieldTPC-=background;
          // double rawYieldErrorTPC = TMath::Sqrt(rawYieldTPC);
          // pad1->SetLogy();
          // frameTPC->SetMinimum(1);
          // frameTPC->SetMaximum(tpcSignalProjection->GetMaximum()*1.5);
          // frameTPC->GetXaxis()->SetRangeUser(lowX,5);
          // frameTPC->Draw();
          // pad2->cd();
          // auto xframeTPC1 = (RooPlot *)frameTPC->Clone();
          // const double *a=tpcSignalProjection->GetXaxis()->GetXbins()->GetArray();
          // TH1D h_diff("h_diff","h_diff",tpcSignalProjection->GetNbinsX(),a);
          // TF1 *f_sig=signalTPC.asTF(nSigmaTPC,RooArgList(meanTPC,sigmaTPC,alpha1TPC,alpha2TPC),nSigmaTPC);
          // auto sig_int = ((RooAbsPdf *)signalTPC.createIntegral(RooArgSet(nSigmaTPC), RooFit::NormSet(RooArgSet(nSigmaTPC)), RooFit::Range("plotRange")))->getVal();
          // TF1 *f_bkg=bkgTPC.asTF(nSigmaTPC,RooArgList(meanTPCBkg,sigmaTPCBkg,tauTPCBkg),nSigmaTPC);
          // auto bkg_int = ((RooAbsPdf *)bkgTPC.createIntegral(RooArgSet(nSigmaTPC), RooFit::NormSet(RooArgSet(nSigmaTPC)), RooFit::Range("plotRange")))->getVal();
          // std::cout<<bkg_int<<std::endl;
          // f_sig->Write();
          // f_bkg->Write();
          // TH1D h_f("h_f","h_f",tpcSignalProjection->GetNbinsX(),a);
          // for (int iB=1; iB<tpcSignalProjection->GetNbinsX()+1; ++iB){
          //   if (tpcSignalProjection->GetBinCenter(iB)>3.) continue;
          //   double eval_sig=f_sig->Eval(tpcSignalProjection->GetBinCenter(iB));
          //   double eval_bkg=f_bkg->Eval(tpcSignalProjection->GetBinCenter(iB));
          //   h_diff.SetBinContent(iB,tpcSignalProjection->GetBinContent(iB)-(nSignalTPC.getVal()*eval_sig*tpcSignalProjection->GetBinWidth(1)/sig_int)-(nBkgTPC.getVal()*eval_bkg*tpcSignalProjection->GetBinWidth(1)/bkg_int));
          //   h_diff.SetBinError(iB,tpcSignalProjection->GetBinError(iB));
          //   h_f.SetBinContent(iB,(nSignalTPC.getVal()*eval_sig*tpcSignalProjection->GetBinWidth(1)/sig_int)+(nBkgTPC.getVal()*eval_bkg*tpcSignalProjection->GetBinWidth(1)/bkg_int));
          // }
          // h_diff.Divide(tpcSignalProjection);
          // h_f.Write();
          // for (int iB=1; iB<tpcSignalProjection->GetNbinsX()+1; ++iB){
          //   if (tpcSignalProjection->GetBinCenter(iB)>3.) continue;
          //   h_diff.SetBinError(iB,tpcSignalProjection->GetBinError(iB)/tpcSignalProjection->GetBinContent(iB));
          // }
          // h_diff.GetXaxis()->SetRangeUser(lowX,5);
          // h_diff.SetMarkerStyle(20);
          // h_diff.SetMarkerSize(1);
          // h_diff.GetYaxis()->SetRangeUser(-1,1);
          // h_diff.Draw("pe");
          // TLine hLine(lowX,0,5,0);
          // hLine.SetLineStyle(kDashed);
          // hLine.Draw("same");
          // c.Write(frameTPC->GetName());
          // c.Print(Form("%s/signal_extraction/%s_%s_%d_%d/TPC_cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));

          tpcSignalProjectionAll->Rebin(1);
          tpcSignalProjection->Rebin(1);
          tpcSignalProjection->SetMarkerStyle(20);
          tpcSignalProjection->SetMarkerSize(1.);
          tpcSignalProjection->SetDrawOption("pe");
          tpcSignalProjection->SetMarkerColor(kBlack);
          tpcSignalProjection->GetXaxis()->SetLabelSize(0.04);
          tpcSignalProjection->GetXaxis()->SetTitleSize(0.04);
          tpcSignalProjection->GetXaxis()->SetTitle("n#sigma_{p}");
          tpcSignalProjection->GetYaxis()->SetTitle(Form("Entries / %.2f",tpcSignalProjection->GetBinWidth(1)));

          tpcSignalProjectionAll->GetXaxis()->SetRangeUser(-15.,-3.);
          double max_bkg = tpcSignalProjectionAll->GetBinCenter(tpcSignalProjectionAll->GetMaximumBin());
          for (int i=0;i<2;++i)tpcSignalProjectionAll->Fit("gaus","RLM+","",max_bkg-2.5,max_bkg+2.5);
          tpcSignalProjectionAll->GetXaxis()->SetRangeUser(-24.,24.);
          double left_fit_range_factor = 2.5;
          if (ptMin > 0.54 && ptMin<0.84) left_fit_range_factor = 1.3; // fits with noitspid do not converge
          else if (ptMin > 0.84) left_fit_range_factor = -1.5;
          double left_range = -999;
          for (int iB=1;iB<tpcSignalProjectionAll->GetNbinsX();++iB){
            if (tpcSignalProjectionAll->GetBinCenter(iB)>(tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2))){
              left_range=tpcSignalProjectionAll->GetXaxis()->GetBinLowEdge(iB);
              break;
            }
          }
          tpcSignalProjectionAll->Write();
          RooRealVar nSigmaTPC("nSigmaTPC","n#sigma_{TPC,p}",left_range,24.,"a.u.");
          RooDataHist hNSigmaTPC("hNSigmaTPC","hNSigmaTPC",RooArgList(nSigmaTPC),tpcSignalProjection);
          RooDataHist hNSigmaTPCAll("hNSigmaTPCAll","hNSigmaTPCAll",RooArgList(nSigmaTPC),tpcSignalProjectionAll);
          RooRealVar meanTPC("#mu","#mu",-2.,2.);
          RooRealVar sigmaTPC("#sigma","#sigma",0.8,2.);
          RooRealVar alpha1TPC("#alpha_{1}","#alpha_{1}",-3.5,-0.8);
          RooRealVar alpha2TPC("#alpha_{2}","#alpha_{2}",0.8,3.5);
          RooRealVar meanTPCBkg("#mu_{Bkg}","#mu_{Bkg}",-20.,-6.);
          if (ptMin>0.7) meanTPCBkg.setMax(-4.);
          if (ptMin>0.71) {
            meanTPCBkg.setMin(-10.);
            meanTPCBkg.setMax(-4.5);
            alpha1TPC.setMin(-3.);
            alpha2TPC.setMax(3.);
          }
          if (ptMin>0.89) {
            meanTPCBkg.setMin(-7.);
            meanTPCBkg.setMax(-2.5);
            alpha1TPC.setMin(-3.);
            alpha1TPC.setMax(-1.);
            alpha2TPC.setMax(3.);
          }
          RooRealVar sigmaTPCBkg("#sigma_{Bkg}","#sigma_{Bkg}",0.5,2.);
          RooRealVar tauTPCBkg("#tau_{Bkg}","#tau_{Bkg}",0.8,2.5);
          RooRealVar tauTPCBkgExp("#tau_{Bkg}","#tau_{Bkg}",-10.,0.);
          RooGausDExp signalTPC("signalTPC","signalTPC",nSigmaTPC,meanTPC,sigmaTPC,alpha1TPC,alpha2TPC);
          RooAbsPdf *bkgTPC=nullptr;
          if (ptMin>0.84) bkgTPC=new RooGausExp("bkgTPC","bkgTPC",nSigmaTPC,meanTPCBkg,sigmaTPCBkg,tauTPCBkg);
          else bkgTPC=new RooExponential("bkgTPC","bkgTPC",nSigmaTPC,tauTPCBkgExp);
          RooRealVar nSignalTPC("n_{Signal}","n_{Signal}",0.,1.e9);
          RooRealVar nBkgTPC("n_{Bkg}","n_{Bkg}",0.,1.e10);
          RooAddPdf modelTPC("modelTPC","modelTPC",RooArgList(signalTPC,*bkgTPC),RooArgList(nSignalTPC,nBkgTPC));
          nSigmaTPC.setRange("fitRange",tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2),2.);
          if (ptMin<0.96){
            if (iMatt == 0){
              meanTPC.setVal(fitParameterMeanTPC[iCent][iPtBin]);
              meanTPC.setConstant(true);
              sigmaTPC.setVal(fitParameterSigmaTPC[iCent][iPtBin]);
              sigmaTPC.setConstant(true);
              alpha1TPC.setVal(fitParameterAlphaLTPC[iCent][iPtBin]);
              alpha1TPC.setConstant(true);
              alpha2TPC.setVal(fitParameterAlphaRTPC[iCent][iPtBin]);
              alpha2TPC.setConstant(true);
              meanTPCBkg.setVal(fitParameterMeanTPCBkg[iCent][iPtBin]);
              sigmaTPCBkg.setVal(fitParameterSigmaTPCBkg[iCent][iPtBin]);
              tauTPCBkg.setVal(fitParameterTauTPCBkg[iCent][iPtBin]);
              tauTPCBkgExp.setVal(fitParameterTauTPCBkgExp[iCent][iPtBin]);
              nSignalTPC.setVal(fitNSigTPC[iCent][iPtBin]);
              nBkgTPC.setVal(fitNBkgTPC[iCent][iPtBin]);
            }
            if (iMatt == 1){
              for (int i=0;i<4;++i)auto r=modelTPC.fitTo(hNSigmaTPCAll,RooFit::Range("fitRange"));
              fitParameterMeanTPC[iCent][iPtBin]=meanTPC.getVal();
              fitParameterSigmaTPC[iCent][iPtBin]=sigmaTPC.getVal();
              fitParameterAlphaLTPC[iCent][iPtBin]=alpha1TPC.getVal();
              fitParameterAlphaRTPC[iCent][iPtBin]=alpha2TPC.getVal();
              fitParameterMeanTPCBkg[iCent][iPtBin]=meanTPCBkg.getVal();
              fitParameterSigmaTPCBkg[iCent][iPtBin]=sigmaTPCBkg.getVal();
              fitParameterTauTPCBkg[iCent][iPtBin]=tauTPCBkg.getVal();
              fitParameterTauTPCBkgExp[iCent][iPtBin]=tauTPCBkgExp.getVal();
              fitNSigTPC[iCent][iPtBin]=nSignalTPC.getVal();
              fitNBkgTPC[iCent][iPtBin]=nBkgTPC.getVal();
              meanTPC.setConstant(true);
              sigmaTPC.setConstant(true);
              alpha1TPC.setConstant(true);
              alpha2TPC.setConstant(true);
            }
          }
          RooPlot *frameTPC_all = nSigmaTPC.frame(RooFit::Name(tpcSignalProjection->GetName()),RooFit::Title(Form("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1])));
          hNSigmaTPCAll.plotOn(frameTPC_all,RooFit::Name("nsigmaTPC"));
          double lowTPCrange=3;
          double tpc_range_limit = (roi_nsigma-8.);
          double tpc_tmp_mean=0., tpc_tmp_rms=0.;
          nSigmaTPC.setRange("signalRange",meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal());
          modelTPC.plotOn(frameTPC_all,RooFit::Components("bkgTPC"),RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
          modelTPC.plotOn(frameTPC_all,RooFit::Components("signalTPC"),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
          modelTPC.plotOn(frameTPC_all,RooFit::Name("model"));
          modelTPC.paramOn(frameTPC_all,RooFit::Parameters(RooArgSet(meanTPC,sigmaTPC)),RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", frameTPC_all->chiSquare("model", "nsigmaTPC"))), RooFit::Layout(0.58812,0.811028,0.561955));
          frameTPC_all->GetXaxis()->SetRangeUser(tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2)+0.5,5);
          frameTPC_all->getAttText()->SetTextFont(44);
          frameTPC_all->getAttText()->SetTextSize(18);
          frameTPC_all->getAttLine()->SetLineWidth(0);
          frameTPC_all->Write();
          for (int i=0;i<2;++i)modelTPC.fitTo(hNSigmaTPC,RooFit::Range("fitRange"));
          nSigmaTPC.setRange("signalRange",meanTPC.getVal()-5.*sigmaTPC.getVal(),meanTPC.getVal()+5.*sigmaTPC.getVal());
          RooPlot *frameTPC = nSigmaTPC.frame(RooFit::Name(tpcSignalProjection->GetName()),RooFit::Title(Form("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1])));
          hNSigmaTPC.plotOn(frameTPC,RooFit::Name("nsigmaTPC"));
          lowTPCrange=3;
          tpc_range_limit = (roi_nsigma-8.);
          tpc_tmp_mean=0.;
          tpc_tmp_rms=0.;
          nSigmaTPC.setRange("signalRange",meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal());
          modelTPC.plotOn(frameTPC,RooFit::Components("bkgTPC"),RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
          modelTPC.plotOn(frameTPC,RooFit::Components("signalTPC"),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
          modelTPC.plotOn(frameTPC,RooFit::Name("model"));
          modelTPC.paramOn(frameTPC,RooFit::Parameters(RooArgSet(meanTPC,sigmaTPC)),RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", frameTPC->chiSquare("model", "nsigmaTPC"))), RooFit::Layout(0.58812,0.811028,0.561955));
          frameTPC->GetXaxis()->SetRangeUser(tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2)+0.5,5);
          frameTPC->getAttText()->SetTextFont(44);
          frameTPC->getAttText()->SetTextSize(18);
          frameTPC->getAttLine()->SetLineWidth(0);
          frameTPC->Write();

          double background=0.;
          if (ptMin > 0.49){
            background = ((RooAbsPdf *)bkgTPC->createIntegral(RooArgSet(nSigmaTPC), RooFit::NormSet(RooArgSet(nSigmaTPC)), RooFit::Range("signalRange")))->getVal();
            background *= nBkgTPC.getVal();
            std::cout<<"background="<<background;
            // TF1* fbkg = nullptr;
            // if (ptMin>0.84)fbkg = bkgTPC->asTF(nSigmaTPC,RooArgList(meanTPCBkg,sigmaTPCBkg,tauTPCBkg));
            // else fbkg = bkgTPC->asTF(nSigmaTPC,RooArgList(tauTPCBkgExp));
            // background = fbkg->Integral(meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal());
          }
          TCanvas c(Form("%s",tpcSignalProjection->GetName()),Form("%s",tpcSignalProjection->GetName()));
          tpcSignalProjection->GetXaxis()->SetRangeUser(tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2)+0.5,5);
          double rawYieldTPC = hNSigmaTPC.sumEntries(Form("nSigmaTPC>%f && nSigmaTPC<%f",meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal()));
          //tpcSignalProjection->GetXaxis()->SetRangeUser(meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal());
          //double rawYieldTPC = tpcSignalProjection->Integral("width");
          double rawYieldErrorTPC = TMath::Sqrt(rawYieldTPC);
          if(ptMin>0.49)rawYieldTPC-=background;
          c.SetLogy();
          frameTPC->SetMinimum(1);
          frameTPC->SetMaximum(tpcSignalProjection->GetMaximum()*1.2);
          frameTPC->Draw();
          c.Print(Form("%s/signal_extraction/%s_%s_%d_%d/TPC_cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));

          // fill raw yield histogram
          fTPCrawYield.SetBinContent(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), rawYieldTPC);
          fTPCrawYield.SetBinError(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), rawYieldErrorTPC);
          if ((ptMin<0.49 || ptMin>0.99 || r->covQual()<2.5 || r->status()!=0)&&(ptMin<0.59||ptMin>0.71)){
            fTPCrawYield.SetBinContent(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), 0.);
            fTPCrawYield.SetBinError(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), 0.);
          }
          
          // tpcSignalProjectionAll->Rebin(1);
          // tpcSignalProjection->Rebin(1);
          // tpcSignalProjection->SetMarkerStyle(20);
          // tpcSignalProjection->SetMarkerSize(1.);
          // tpcSignalProjection->SetDrawOption("pe");
          // tpcSignalProjection->SetMarkerColor(kBlack);
          // tpcSignalProjection->GetXaxis()->SetLabelSize(0.04);
          // tpcSignalProjection->GetXaxis()->SetTitleSize(0.04);
          // tpcSignalProjection->GetXaxis()->SetTitle("n#sigma_{p}");
          // tpcSignalProjection->GetYaxis()->SetTitle(Form("Entries / %.2f",tpcSignalProjection->GetBinWidth(1)));

          // tpcSignalProjectionAll->GetXaxis()->SetRangeUser(-15.,-3.);
          // double max_bkg = tpcSignalProjectionAll->GetBinCenter(tpcSignalProjectionAll->GetMaximumBin());
          // tpcSignalProjectionAll->Fit("gaus","RLM+","",max_bkg-2.5,max_bkg+2.5);
          // tpcSignalProjectionAll->GetXaxis()->SetRangeUser(-24.,24.);
          // double left_fit_range_factor = 2.5;
          // if (ptMin > 0.54 && ptMin<0.84) left_fit_range_factor = 1.5;
          // else if (ptMin > 0.84) left_fit_range_factor = -1.5;
          // double left_range = -999;
          // for (int iB=1;iB<tpcSignalProjectionAll->GetNbinsX();++iB){
          //   if (tpcSignalProjectionAll->GetBinCenter(iB)>tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2)){
          //     left_range=tpcSignalProjectionAll->GetXaxis()->GetBinLowEdge(iB);
          //     break;
          //   }
          // }
          // RooRealVar nSigmaTPC("nSigmaTPC","n#sigma_{TPC,p}",left_range,24.,"a.u.");
          // RooDataHist hNSigmaTPC("hNSigmaTPC","hNSigmaTPC",RooArgList(nSigmaTPC),tpcSignalProjection);
          // RooDataHist hNSigmaTPCAll("hNSigmaTPCAll","hNSigmaTPCAll",RooArgList(nSigmaTPC),tpcSignalProjectionAll);
          // RooRealVar meanTPC("#mu","#mu",-2.,2.);
          // RooRealVar sigmaTPC("#sigma","#sigma",0.8,3.);
          // RooRealVar alpha1TPC("#alpha_{1}","#alpha_{1}",-3.,-0.8);
          // RooRealVar alpha2TPC("#alpha_{2}","#alpha_{2}",0.8,3.);
          // RooRealVar meanTPCBkg("#mu_{Bkg}","#mu_{Bkg}",-20.,-6.);
          // if (ptMin>0.7) meanTPCBkg.setMax(-5.);
          // if (ptMin>0.8) {
          //   meanTPCBkg.setMin(-8.);
          //   meanTPCBkg.setMax(-3.5);
          // }
          // RooRealVar sigmaTPCBkg("#sigma_{Bkg}","#sigma_{Bkg}",0.5,2.);
          // RooRealVar tauTPCBkg("#tau_{Bkg}","#tau_{Bkg}",0.5,2.5);
          // RooGausDExp signalTPC("signalTPC","signalTPC",nSigmaTPC,meanTPC,sigmaTPC,alpha1TPC,alpha2TPC);
          // RooGausExp bkgTPC("bkgTPC","bkgTPC",nSigmaTPC,meanTPCBkg,sigmaTPCBkg,tauTPCBkg);
          // RooRealVar nSignalTPC("n_{Signal}","n_{Signal}",0.,1.e9);
          // RooRealVar nBkgTPC("n_{Bkg}","n_{Bkg}",0.,1.e9);
          // RooAddPdf modelTPC("modelTPC","modelTPC",RooArgList(signalTPC,bkgTPC),RooArgList(nSignalTPC,nBkgTPC));
          // nSigmaTPC.setRange("fitRange",tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2),2.);
          // for (int i=0;i<2;++i)modelTPC.fitTo(hNSigmaTPCAll,RooFit::Range("fitRange"));
          // for (int i=0;i<2;++i)modelTPC.fitTo(hNSigmaTPC,RooFit::Range("fitRange"));
          // nSigmaTPC.setRange("signalRange",meanTPC.getVal()-5.*sigmaTPC.getVal(),meanTPC.getVal()+5.*sigmaTPC.getVal());
          // RooPlot *frameTPC = nSigmaTPC.frame(RooFit::Name(tpcSignalProjection->GetName()),RooFit::Title(Form("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1])));
          // hNSigmaTPC.plotOn(frameTPC,RooFit::Name("nsigmaTPC"));
          // double lowTPCrange=3;
          // double tpc_range_limit = (roi_nsigma-8.);
          // double tpc_tmp_mean=0., tpc_tmp_rms=0.;
          // nSigmaTPC.setRange("signalRange",meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal());
          // modelTPC.plotOn(frameTPC,RooFit::Components("bkgTPC"),RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
          // modelTPC.plotOn(frameTPC,RooFit::Components("signalTPC"),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
          // modelTPC.plotOn(frameTPC,RooFit::Name("model"));
          // modelTPC.paramOn(frameTPC,RooFit::Parameters(RooArgSet(meanTPC,sigmaTPC)),RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", frameTPC->chiSquare("model", "nsigmaTPC"))), RooFit::Layout(0.58812,0.811028,0.561955));
          // frameTPC->GetXaxis()->SetRangeUser(tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2)+0.5,5);
          // frameTPC->getAttText()->SetTextFont(44);
          // frameTPC->getAttText()->SetTextSize(18);
          // frameTPC->getAttLine()->SetLineWidth(0);
          // frameTPC->Write();

          // double background=0.;
          // if (ptMin > 0.49){
          //   background = ((RooAbsPdf *)bkgTPC.createIntegral(RooArgSet(nSigmaTPC), RooFit::NormSet(RooArgSet(nSigmaTPC)), RooFit::Range("signalRange")))->getVal();
          //   background *= nBkgTPC.getVal();
          // }
          // TCanvas c(Form("%s",tpcSignalProjection->GetName()),Form("%s",tpcSignalProjection->GetName()));
          // tpcSignalProjection->GetXaxis()->SetRangeUser(tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(1)+left_fit_range_factor*tpcSignalProjectionAll->GetFunction("gaus")->GetParameter(2)+0.5,5);
          // double rawYieldTPC = hNSigmaTPC.sumEntries(Form("nSigmaTPC>%f && nSigmaTPC<%f",meanTPC.getVal()-(lowTPCrange+tpc_range_limit)*sigmaTPC.getVal(),meanTPC.getVal()+(5.+tpc_range_limit)*sigmaTPC.getVal()));
          // if(ptMin>0.49)rawYieldTPC-=background;
          // double rawYieldErrorTPC = TMath::Sqrt(rawYieldTPC);
          // c.SetLogy();
          // frameTPC->SetMinimum(1);
          // frameTPC->SetMaximum(tpcSignalProjection->GetMaximum()*1.2);
          // frameTPC->Draw();
          // c.Print(Form("%s/signal_extraction/%s_%s_%d_%d/TPC_cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));

          // // fill raw yield histogram
          // fTPCrawYield.SetBinContent(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), rawYieldTPC);
          // fTPCrawYield.SetBinError(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), rawYieldErrorTPC);
          // if (ptMin>0.99){
          //   fTPCrawYield.SetBinContent(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), 0.);
          //   fTPCrawYield.SetBinError(fTPCrawYield.FindBin((ptMax + ptMin) / 2.), 0.);
          // }

          // fill raw yield histogram
          fTOFrawYield.SetBinContent(fTOFrawYield.FindBin((ptMax + ptMin) / 2.), rawYield);
          fTOFrawYield.SetBinError(fTOFrawYield.FindBin((ptMax + ptMin) / 2.), rawYieldError);

          // fill significance histogram
          if (binCounting)
            fSignificance.SetBinContent(fSignificance.FindBin((ptMax + ptMin) / 2.), rawYield / rawYieldError);
          else
            //fSignificance.SetBinContent(fSignificance.FindBin((ptMax + ptMin) / 2.), rawYield / TMath::Sqrt(rawYield + bkgIntegral_val));
            fSignificance.SetBinError(fSignificance.FindBin((ptMax + ptMin) / 2.), 0.);

          // fill sigma histogram
          fSigma.SetBinContent(fSigma.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParameter(2));
          fSigma.SetBinError(fSigma.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParError(2));

          // fill mean histogram
          fMean.SetBinContent(fMean.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParameter(1));
          fMean.SetBinError(fMean.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParError(1));

          // fill alphaR histogram
          fAlphaR.SetBinContent(fAlphaR.FindBin((ptMax + ptMin) / 2.), alphaR->getVal());
          fAlphaR.SetBinError(fAlphaR.FindBin((ptMax + ptMin) / 2.), alphaR->getError());

          // fill alphaL histogram
          fAlphaL.SetBinContent(fAlphaL.FindBin((ptMax + ptMin) / 2.), alphaL->getVal());
          fAlphaL.SetBinError(fAlphaL.FindBin((ptMax + ptMin) / 2.), alphaL->getError());
        
        }
        // save to pdf
        TCanvas canv("canv", "canv");
        canv.cd();
        TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.4, 1.0, 1.0, 0);
        pad1->SetFillColor(0);
        pad1->SetFrameBorderMode(0);
        pad1->SetBottomMargin(0.06);
        pad1->SetLogy();
        pad1->Draw();
        TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.4, 0);
        pad2->SetFillColor(0);
        pad2->SetFrameBorderMode(0);
        pad2->SetFillColor(0);
        pad2->SetFrameBorderMode(0);
        pad2->SetTopMargin(0.10);
        pad2->SetBottomMargin(0.25);
        pad2->Draw();
        pad1->cd();
        xframe->GetXaxis()->SetLabelOffset(0.005);
        xframe->GetXaxis()->SetTitle("");
        xframe->GetXaxis()->SetTitleOffset(-0.005);
        xframe->Draw("");
        xframe->GetXaxis()->SetLabelSize(0.06);
        xframe->GetYaxis()->SetLabelSize(0.06);
        xframe->GetYaxis()->SetTitleSize(0.07);
        xframe->GetYaxis()->SetTitleOffset(0.72);
        tofSignalProjection->GetXaxis()->SetRangeUser(-0.5, 0.5);
        double peakMinimum = tofSignalProjection->GetBinContent(tofSignalProjection->GetMinimumBin());
        double peakMaximum = tofSignalProjection->GetBinContent(tofSignalProjection->GetMaximumBin());
        xframe->GetYaxis()->SetRangeUser(10.,peakMaximum+1e3);
        if (kVerbose) std::cout << "peakMaximum=" << peakMaximum << std::endl;
        TLine lsx(mean_tmp - (roi_nsigma_down+extend_roi) * rms_tmp, 0, mean_tmp - (roi_nsigma_down+extend_roi) * rms_tmp, peakMaximum);
        lsx.SetLineStyle(kDashed);
        lsx.Draw("same");
        TLine ldx(mean_tmp + (roi_nsigma_up+extend_roi) * rms_tmp, 0, mean_tmp + (roi_nsigma_up+extend_roi) * rms_tmp, peakMaximum*0.03);
        //if(ptMin>1.19)ldx.SetY2(peakMaximum*0.6);
        ldx.SetLineStyle(kDashed);
        ldx.Draw("same");
        canv.SetName(plotTitle);
        TLatex covarianceQuality;
        covarianceQuality.SetNDC();
        covarianceQuality.SetX(0.55);
        covarianceQuality.SetY(0.4);
        covarianceQuality.SetTitle(Form("Covariance quality = %d", covQ));
        pad1->Update();

        RooHist *hresid;
        RooPlot *frame2;

        if (extractSignal)
        {
          auto xframe1 = (RooPlot *)xframe->Clone();
          hresid = xframe1->residHist();
          frame2 = tofSignal.frame(RooFit::Title(""));
          frame2->addPlotable(hresid, "P");
          frame2->GetYaxis()->SetMaxDigits(2);
        }

        tofSignalProjection->Write();
        tpcSignalProjection->Write();
        xframe->Write();

        pad2->cd();
        frame2->Draw("");
        frame2->GetXaxis()->SetLabelOffset(0.005);
        frame2->GetXaxis()->SetTitle("n#sigma_{p} (a.u.)");
        frame2->GetXaxis()->SetTitleOffset(0.9);
        frame2->GetYaxis()->SetNdivisions(5);
        frame2->GetXaxis()->SetTitleSize(0.11);
        frame2->GetXaxis()->SetLabelSize(0.095);
        frame2->GetYaxis()->SetLabelSize(0.09);
        //frame2->GetYaxis()->SetRangeUser(-1.e4, 2.e4);
        frame2->GetYaxis()->SetTitleSize(0.095);
        frame2->GetYaxis()->SetTitleOffset(0.52);
        frame2->GetYaxis()->SetTitle(Form("#frac{ Ev. - Bg. }{ %.1f a.u. }", tofSignalProjection->GetXaxis()->GetBinWidth(2)));
        frame2->SetTitle("   ");
        //covarianceQuality.Draw("same");
        //canv.SetLogy();
        TLine lsx1(mean_tmp - roi_nsigma_down * rms_tmp, -1.e4, mean_tmp - roi_nsigma_down * rms_tmp, 2.e4);
        lsx1.SetLineStyle(kDashed);
        //lsx1.Draw("same");
        TLine ldx1(mean_tmp + roi_nsigma_up * rms_tmp, -1.e4, mean_tmp + roi_nsigma_up * rms_tmp, 2.e4);
        ldx1.SetLineStyle(kDashed);
        //ldx1.Draw("same");
        TLine zero(nSigmaLeft, 0, maxNsigma, 0);
        zero.SetLineStyle(kDashed);
        zero.Draw("same");
        pad2->Update();
        canv.Write(Form("%s_%s_%d_%d_cent_%.0f_%.0f_pt_%.2f_%.2f", kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));
        canv.Print(Form("%s/signal_extraction/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));
        tofSignalProjection->SetMarkerStyle(20);
        tofSignalProjection->SetMarkerSize(0.8);
        TCanvas canv1("c1", "c1");
        tofSignalProjection->Draw("pe");
        canv1.Print(Form("%s/signal_extraction_preliminary/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));
      } // end of loop on pt bins

      // set raw yield histogram style and write to file
      fTPCrawYield.SetTitle(TString::Format("%s raw yield, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fTPCrawYield.SetName(TString::Format("f%sTPCrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fTPCrawYield.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fTPCrawYield.GetYaxis()->SetTitle("#it{N_{raw}}");
      fTPCrawYield.SetMarkerStyle(20);
      fTPCrawYield.SetMarkerSize(0.8);
      fTPCrawYield.SetOption("PE");
      fTPCrawYield.Write();

      // set raw yield histogram style and write to file
      fTOFrawYield.SetTitle(TString::Format("%s raw yield, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fTOFrawYield.SetName(TString::Format("f%sTOFrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fTOFrawYield.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fTOFrawYield.GetYaxis()->SetTitle("#it{N_{raw}}");
      fTOFrawYield.SetMarkerStyle(20);
      fTOFrawYield.SetMarkerSize(0.8);
      fTOFrawYield.SetOption("PE");
      fTOFrawYield.Write();

      // set significance histogram style and write to file
      fSignificance.SetTitle(TString::Format("%s significance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fSignificance.SetName(TString::Format("f%sSignificance_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fSignificance.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSignificance.GetYaxis()->SetTitle("S/#sqrt{N}");
      fSignificance.SetMarkerStyle(20);
      fSignificance.SetMarkerSize(0.8);
      fSignificance.SetOption("PE");
      fSignificance.Write();

      // set sigma histogram style and write to file
      fSigma.SetTitle(TString::Format("%s sigma, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fSigma.SetName(TString::Format("f%sSigma_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fSigma.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSigma.GetYaxis()->SetTitle("#sigma (GeV^2/#it{c}^{4})");
      fSigma.SetMarkerStyle(20);
      fSigma.SetMarkerSize(0.8);
      fSigma.SetOption("PE");
      fSigma.Write();

      // set mean histogram style and write to file
      fMean.SetTitle(TString::Format("%s mean, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fMean.SetName(TString::Format("f%sMean_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fMean.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fMean.GetYaxis()->SetTitle("#mu (GeV^2/#it{c}^{4})");
      fMean.SetMarkerStyle(20);
      fMean.SetMarkerSize(0.8);
      fMean.SetOption("PE");
      fMean.Write();

      // set sigma histogram style and write to file
      fAlphaR.SetTitle(TString::Format("%s alphaR, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fAlphaR.SetName(TString::Format("f%sAlphaR_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fAlphaR.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fAlphaR.GetYaxis()->SetTitle("#alpha_{R}");
      fAlphaR.SetMarkerStyle(20);
      fAlphaR.SetMarkerSize(0.8);
      fAlphaR.SetOption("PE");
      fAlphaR.Write();

      // set sigma histogram style and write to file
      fAlphaL.SetTitle(TString::Format("%s alphaL, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fAlphaL.SetName(TString::Format("f%sAlphaL_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]));
      fAlphaL.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fAlphaL.GetYaxis()->SetTitle("#alpha_{L}");
      fAlphaL.SetMarkerStyle(20);
      fAlphaL.SetMarkerSize(0.8);
      fAlphaL.SetOption("PE");
      fAlphaL.Write();

      //delete model, mean, sigma, alphaL, alphaR;
    } // end of loop on centrality bin
  }   // end of loop on antimatter/matter
  outFile->Close();
} // end of macro Signal.cpp
