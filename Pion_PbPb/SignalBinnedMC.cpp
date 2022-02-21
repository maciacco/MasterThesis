// SignalUnbinnedMC.cpp

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

#include "../utils/Utils.h"
#include "../utils/Config.h"
#include "../utils/RooGausExp.h"
#include "../utils/RooDSCBShape.h"
#include "../utils/RooGausDExp.h"

using namespace utils;
using namespace pion;

double roi_nsigma_up = 8.;
double roi_nsigma_down = 8.;

double computeExponentialNormalisation(double, double, double);

const double kNSigma = 3; // define interval for bin counting

void SignalBinnedMC(const char *cutSettings = "", const double roi_min_limit_input = 1.5, const double roi_max_limit_input = 11., const double mismatch_min_limit_input = 8.5, const double mismatch_max_limit_input = 13.5, const bool binCounting = false, const int bkg_shape = 1, const char *inFileDat = "mc_20g7_likeData_largeNsigma", const char *outFileName = "SignalProton", const char *outFileOption = "recreate", const bool extractSignal = true, const bool useDSCB = false, const bool binCountingNoFit = false)
{
  //double roi_max_limit=12.;

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
  //TFile *mcFile_21l5 = TFile::Open(TString::Format("%s/%s_largeNsigma.root", kDataDir, inFileDat)); // open data TFile
  TFile *mcFile_21l5 = TFile::Open(TString::Format("%s/%s.root", kDataDir, inFileDat)); // open data TFile
  TFile *mcFile_20g7 = TFile::Open(TString::Format("%s/%s.root", kDataDir, "mc_20g7_likeData_largeNsigma")); // open data TFile

  //TFile *dataFile2 = TFile::Open(TString::Format("%s/%s.root", kDataDir, inFileDat)); // open data TFile
  if (!mcFile_21l5)
  {
    std::cout << "File not found!" << std::endl; // check data TFile opening
    return;
  }

  // get TTList
  std::string listName_21l5 = Form("nuclei_pion_mcFalse_%s", cutSettings);
  TTList *list_21l5 = (TTList *)mcFile_21l5->Get(listName_21l5.data());
  // std::string listName_20g7 = Form("nuclei_proton_%s", cutSettings);
  // TTList *list_20g7 = (TTList *)mcFile_20g7->Get(listName_20g7.data());
  //TTList *list2 = (TTList *)dataFile2->Get(listName.data());

  // merge antimatter + matter histograms
  std::string histNameA = Form("f%sTOFnSigma", kAntimatterMatter[0]);
  std::string histNameM = Form("f%sTOFnSigma", kAntimatterMatter[1]);
  TH3F *fTOFSignalA = (TH3F *)list_21l5->Get(histNameA.data());
  //TH3F *fTOFSignalA_20g7 = (TH3F *)list_20g7->Get(histNameA.data());
  /* if (ADD20g7)
    fTOFSignalA->Add(fTOFSignalA_20g7); */
  // TH3F *fTOFSignalA2 = (TH3F *)list2->Get(histNameA.data());
  TH3F *fTOFSignalAll = (TH3F *)fTOFSignalA->Clone(fTOFSignalA->GetName());
  TH3F *fTOFSignalM = (TH3F *)list_21l5->Get(histNameM.data());
  //TH3F *fTOFSignalM_20g7 = (TH3F *)list_20g7->Get(histNameM.data());
  /* if (ADD20g7)
    fTOFSignalM->Add(fTOFSignalM_20g7); */
  //TH3F *fTOFSignalM2 = (TH3F *)list2->Get(histNameM.data());
  //fTOFSignalAll->Add(fTOFSignalA2);
/*   if (ADD20g7)
    fTOFSignalAll->Add(fTOFSignalM); */
  //fTOFSignalAll->Add(fTOFSignalM2);

  for (int iMatt = 1; iMatt > -1; --iMatt)
  { // loop on antimatter/matter
    // Get histograms from file
    std::string histName = Form("f%sTOFnSigma", kAntimatterMatter[iMatt]);
    TH3F *fTOFSignal1 = (TH3F *)list_21l5->Get(histName.data());
    fTOFSignal1->SetName("A");
    // TH3F *fTOFSignal2 = (TH3F *)list_20g7->Get(histName.data());
    // fTOFSignal1->SetName("B");
    if (!fTOFSignal1)
    {
      std::cout << "Hstogram not found!" << std::endl; // check data TFile opening
      return;
    }

    TH3F *fTOFSignal = (TH3F *)fTOFSignal1->Clone(histName.data());
    // if (ADD20g7)
    //   fTOFSignal->Add(fTOFSignal2);

    // make plot subdirectory
    system(Form("mkdir %s/signal_extraction/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));
    system(Form("mkdir %s/signal_extraction_preliminary/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));

    dirOutFile->cd();

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {                                                                          // loop over centrality classes
      TH1D fTOFrawYield("fRawYield", "fRawYield", kNPtBins, kPtBins);          // declare raw yields histogram
      TH1D fSignificance("fSignificance", "fSignificance", kNPtBins, kPtBins); // declare significance
      TH1D fSigma("fSigma", "fSigma", kNPtBins, kPtBins);
      TH1D fMean("fMean", "fMean", kNPtBins, kPtBins);
      TH1D fAlphaL("fAlphaL", "fAlphaL", kNPtBins, kPtBins);
      TH1D fAlphaR("fAlphaR", "fAlphaR", kNPtBins, kPtBins);
      int nUsedPtBins = 27; // up to 2.00 GeV/c
      //int nUsedPtBins = 39;

      for (int iPtBin = 7; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        double ptMin = fTOFrawYield.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fTOFrawYield.GetXaxis()->GetBinUpEdge(iPtBin);
        double centMin = kCentBinsLimitsPion[iCent][0];
        double centMax = kCentBinsLimitsPion[iCent][1];

        int pTbinsIndexMin = fTOFSignal->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fTOFSignal->GetYaxis()->FindBin(ptMax - 0.005);
        int centBinMin = fTOFSignal->GetXaxis()->FindBin(kCentBinsLimitsPion[iCent][0]);      //kCentBinsPion[iCent][0];
        int centBinMax = fTOFSignal->GetXaxis()->FindBin(kCentBinsLimitsPion[iCent][1] - 1.); //kCentBinsPion[iCent][1];

        // project histogram
        std::cout << "Pt bins: min=" << ptMin << "; max=" << ptMax << "; indexMin=" << pTbinsIndexMin << "; indexMax=" << pTbinsIndexMax << "; centralityBinMin=" << centBinMin << "; centralityBinMax=" << centBinMax << std::endl;
        TH1D *tofSignalProjection = fTOFSignal->ProjectionZ(Form("f%sTOFSignal_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        // tofSignalProjection->Rebin(4);

        // project histogram (antimatter + matter)
        TH1D *tofSignalProjectionAll = fTOFSignalAll->ProjectionZ(Form("f%sTOFSignalAll_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignalAll->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignalAll->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        // tofSignalProjectionAll->Rebin(1);

        // limits
        double nSigmaLeft = -10;
        double nSigmaRight = -5.;
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(nSigmaLeft, kTOFnSigmaMax);
        tofSignalProjection->GetXaxis()->SetRangeUser(nSigmaLeft, kTOFnSigmaMax);

        // DEFINE SIGNAL REGION
        TF1 signalRegionFit("signalRegionFit", "gaus", -20., 20.);
        signalRegionFit.SetParLimits(0, 1., 1.e7);
        signalRegionFit.SetParLimits(1, -0.4, 0.4);
        signalRegionFit.SetParLimits(2, 0.8, 1.6);
        signalRegionFit.SetLineColor(kBlue);
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-.5, .5);
        double maximum_signal = tofSignalProjectionAll->GetBinCenter(tofSignalProjectionAll->GetMaximumBin());
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-20., 20.);
        tofSignalProjectionAll->Fit("signalRegionFit", "QRL+", "", maximum_signal - 1., maximum_signal + 1.);
        double mean_tmp = signalRegionFit.GetParameter(1);
        double rms_tmp = signalRegionFit.GetParameter(2);
        double roi_max_limit = mean_tmp+(roi_max_limit_input-2.)*rms_tmp; // smaller signal window in MC than in data
        if (ptMin>1.09)
          roi_max_limit = mean_tmp+(roi_max_limit_input-3.)*rms_tmp;
        if (ptMin>1.49)
          roi_max_limit = mean_tmp+(roi_max_limit_input-5.)*rms_tmp;

        // DEFINE K SIGNAL REGION -> MISMATCH FIT REGION
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-50., 50.);
        TF1 signalRegionFitK("signalRegionFitK", "gaus", -50., 50.);
        signalRegionFitK.SetParLimits(0, 1., 1.e7);
        signalRegionFitK.SetParLimits(1, 0., 20.);
        signalRegionFitK.SetParLimits(2, 0.8, 3.);
        signalRegionFitK.SetLineColor(kBlue);
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(5., 30.);
        if (ptMin > 1.54)tofSignalProjectionAll->GetXaxis()->SetRangeUser(5., 15.);
        double maximum_signal_K = tofSignalProjectionAll->GetBinCenter(tofSignalProjectionAll->GetMaximumBin());
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-20., 20.);
        tofSignalProjectionAll->Fit("signalRegionFitK", "QRL+", "", maximum_signal_K - 1., maximum_signal_K + 1.);
        double mean_tmp_K = signalRegionFitK.GetParameter(1);
        double rms_tmp_K = signalRegionFitK.GetParameter(2);
        double mismatch_left_limit = mean_tmp_K+mismatch_min_limit_input*rms_tmp_K;
        double mismatch_right_limit = mean_tmp_K+mismatch_max_limit_input*rms_tmp_K;

        // roofit data
        double maxNsigma=17.;
        if (ptMin>0.99)
        maxNsigma=18.;
        if (ptMin>1.04)
        maxNsigma=16.5;
        if (ptMin>1.14)
        maxNsigma=15.;
        if (ptMin>1.19)
        maxNsigma=13.;
        if (ptMin>1.24){
          maxNsigma=mean_tmp_K-1.*rms_tmp_K;
        }
        if (ptMin>1.54){
          maxNsigma=mean_tmp_K-.5*rms_tmp_K;
        }
        RooRealVar tofSignal("tofSignal", "n#sigma_{p}", nSigmaLeft, maxNsigma, "a.u.");
        RooRealVar *tofSignal_full=new RooRealVar("tofSignal_full", "n#sigma_{p}full_", nSigmaLeft, 20., "a.u.");
        if (ptMin>1.09)tofSignal_full=new RooRealVar("tofSignal_full", "n#sigma_{p}full_", nSigmaLeft, mismatch_right_limit, "a.u.");
        
        // roofit histogram
        RooDataHist data("data", "data", RooArgList(tofSignal), tofSignalProjection);
        RooDataHist data_full("data_full", "data_full", RooArgList(*tofSignal_full), tofSignalProjection);

        // roofit histogram
        RooDataHist dataAll("dataAll", "dataAll", RooArgList(tofSignal), tofSignalProjectionAll);
        std::cout << "Number of entries (roofit) = " << dataAll.sumEntries() << std::endl;
        std::cout << "Number of entries (root) = " << tofSignalProjectionAll->GetEntries() << std::endl;

        RooAbsPdf *background0;
        RooAbsPdf *background1;
        RooAbsPdf *background2;
        RooRealVar *slope1;
        RooRealVar *slope2;
        RooRealVar *nBackground1;
        RooRealVar *nBackground2;
        RooRealVar *nBackground3;
        RooAddPdf *model;
        RooAbsPdf *modelAll;
        RooAbsPdf *modelMismatch;
        RooAbsPdf *modelTail;
        RooFitResult *r;

        slope1 = new RooRealVar("#tau_{1}", "slope1", -0.1, -.7, 0.);
        nBackground1 = new RooRealVar("#it{N}_{Bkg,1}", "nBackground1", 0.,1.);
        if (bkg_shape == 1 && ptMin < 0.99)
        { // expo
          background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
          modelAll = (RooAbsPdf*)new RooExponential("modelAll", "modelAll", tofSignal, *slope1);
        }
        else if (bkg_shape == 1 && ptMin > 0.99 && ptMin < 1.09)
        { // double expo
          slope2 = new RooRealVar("#tau_{2}", "slope2", 0., 2.);
          background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
          background2 = (RooAbsPdf *)new RooExponential("background2", "background2", tofSignal, *slope2);
          modelAll = (RooAbsPdf*)new RooAddPdf("modelAll", "modelAll", RooArgList(*background1,*background2), RooArgList(*nBackground1));
        }
        else if (bkg_shape == 1 && ptMin > 1.09)
        { // double expo with fixed mismatch
          slope1 = new RooRealVar("#tau_{1}", "slope1", -0.1, -.7, 0.7);
          slope2 = new RooRealVar("#tau_{2}", "slope2", 0., 2.);
          background0 = (RooAbsPdf *)new RooExponential("background1", "background1", *tofSignal_full, *slope1);
          background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
          background2 = (RooAbsPdf *)new RooExponential("background2", "background2", tofSignal, *slope2);
          nBackground3 = new RooRealVar("#it{N}_{Bkg,3}", "nBackground3", 0.,1.e8);
          nBackground1 = new RooRealVar("#it{N}_{Bkg,1}", "nBackground1", 0., 1.e10);
          modelMismatch = (RooAddPdf*)new RooAddPdf("model", "model", RooArgList(*background0), RooArgList(*nBackground1));
        }
        else
          std::cout << "!!!!!" << std::endl; // TODO: UPDATE FOLLOWING BLOCK

        int covQ = -999;
        double intersectionBinCenter=-20.;
        intersectionBinCenter=mean_tmp-roi_min_limit_input*rms_tmp;
        double signalRightLimit=roi_max_limit;/* -1.;
        if(ptMin>0.89)
        signalRightLimit=roi_max_limit-2.;
        if(ptMin>1.09)
        signalRightLimit=roi_max_limit-3.;
        if(ptMin>1.14)
        signalRightLimit=roi_max_limit-5.; */
        
        if (extractSignal)
        {

          if (ptMin < 1.09){
            tofSignal.setRange("rightSideband", signalRightLimit, maxNsigma);
            
            modelAll->fitTo(dataAll, RooFit::Range("rightSideband"));
            r = modelAll->fitTo(dataAll, RooFit::Save(), RooFit::Range("rightSideband"));
            
            slope1->setConstant();
            if (ptMin > 0.99)
              slope2->setConstant();
            nBackground1->setConstant();

            if (bkg_shape == 1 && ptMin < 0.99)
            { // expo
              nBackground2 = new RooRealVar("#it{N}_{Bkg}", "nBackground2", 1., 1.e9);
              model = new RooAddPdf("model", "model", RooArgList(*modelAll), RooArgList(*nBackground2));
            }
            else if (bkg_shape == 1 && ptMin > 0.99)
            { // double expo
              nBackground2 = new RooRealVar("#it{N}_{Bkg}", "nBackground2", 1., 1.e9);
              model = new RooAddPdf("model", "model", RooArgList(*modelAll), RooArgList(*nBackground2));
            }
            model->fitTo(data, RooFit::Range("rightSideband"));
            r = model->fitTo(data, RooFit::Save(), RooFit::Range("rightSideband"));
          }
          else{
            // pion sideband
            tofSignal.setRange("rightSideband", mean_tmp_K-3.*rms_tmp_K, maxNsigma);

            // K sideband
            tofSignal_full->setRange("rightSidebandK", mismatch_left_limit, mismatch_right_limit);

            // full range
            tofSignal_full->setRange("full", nSigmaLeft,mismatch_right_limit);

            // narrower range
            tofSignal_full->setRange("small", nSigmaLeft,maxNsigma);

            for (int I=0;I<2;++I)modelMismatch->fitTo(data_full, RooFit::Save(), RooFit::Range("rightSidebandK"));
            double mismatch_integral_=((RooAbsPdf *)modelMismatch->createIntegral(RooArgSet(*tofSignal_full), RooFit::NormSet(RooArgSet(*tofSignal_full)), RooFit::Range("full")))->getVal();
            double mismatch_integral_small=((RooAbsPdf *)modelMismatch->createIntegral(RooArgSet(*tofSignal_full), RooFit::NormSet(RooArgSet(*tofSignal_full)), RooFit::Range("small")))->getVal();
            std::cout << "mism_integral_full = " << mismatch_integral_ << "; mism_integral_small = " << mismatch_integral_small << std::endl;
            slope1->setConstant();
            nBackground2 = new RooRealVar("#it{N}_{Bkg}", "nBackground", 0., 1.e10);
            nBackground2->setVal(nBackground1->getVal()*mismatch_integral_small);
            nBackground2->setConstant();
            model = new RooAddPdf("model", "model", RooArgList(*background1,*background2), RooArgList(*nBackground2,*nBackground3));
            for (int I=0;I<2;++I)r = model->fitTo(data, RooFit::Save(), RooFit::Range("rightSideband"));
            /* slope1->setConstant(false);
            nBackground2->setConstant(false); */
          }
          
          // find signal region
          double leftSignalLimit=-10;
          bool intersection=false;
          int iB=1;
          
          int binShiftIndex=tofSignalProjection->FindBin(leftSignalLimit);
          /* while (!intersection){
            double pdfValue=0;
            double binContent=data.weight(iB);
            if (ptMin<0.99){
              TF1* bkgFunction=modelAll->asTF(tofSignal,RooArgList(*slope1));
              double normBkgFunction=nBackground2->getVal() * tofSignalProjection->GetBinWidth(1) / bkgFunction->Integral(-10.,maxNsigma);
              pdfValue=normBkgFunction*bkgFunction->Eval(tofSignalProjection->GetBinCenter(iB+binShiftIndex));
            }
            else if (ptMin>0.99) {
              TF1* bkgFunction1=background1->asTF(tofSignal,RooArgList(*slope1));
              TF1* bkgFunction2=background2->asTF(tofSignal,RooArgList(*slope2));
              double val1=bkgFunction1->Eval(tofSignalProjection->GetBinCenter(iB+binShiftIndex));
              double val2=bkgFunction2->Eval(tofSignalProjection->GetBinCenter(iB+binShiftIndex));
              double frac=nBackground1->getVal();
              double integral_1 = bkgFunction1->Integral(-10.,maxNsigma);
              double integral_2 = bkgFunction2->Integral(-10.,maxNsigma);
              double val_tot = val1/integral_1+frac*val2/integral_2;
              double normBkgFunction=nBackground2->getVal() * tofSignalProjection->GetBinWidth(1);
              pdfValue=normBkgFunction*val_tot;
            }
            if(binContent>pdfValue)
            {
              intersection=true;
              std::cout << "bin = " << iB << "; weight = " << binContent << "; pdf = " << pdfValue << std::endl;
            }
            else iB++;
          } */
          //intersectionBinCenter=tofSignalProjection->GetBinCenter(iB+binShiftIndex);
          std::cout << "intersection bin center = " << intersectionBinCenter << std::endl;
          tofSignal.setRange("signalRange", intersectionBinCenter, signalRightLimit);

          std::cout << "fit status: " << r->status() << ";" << std::endl;
          std::cout << "covariance quality: " << r->covQual() << std::endl;
          if (r->status()!=0) continue;
          covQ = r->covQual();
        }

        // frame
        int nBins = tofSignalProjection->GetNbinsX();
        TString plotTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]);
        RooPlot *xframe_full = tofSignal_full->frame(RooFit::Bins(nBins), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigmaFull_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], ptMin, ptMax)));
        RooPlot *xframe = tofSignal.frame(RooFit::Bins(nBins), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigma_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], ptMin, ptMax)));
        data.plotOn(xframe, RooFit::Name("dataNsigma"));
        data_full.plotOn(xframe_full, RooFit::Name("dataNsigma"));
        tofSignal.setRange("full", nSigmaLeft, 20);
        tofSignal.setRange("model", nSigmaLeft, maxNsigma);
        if (extractSignal)
        {
          if (ptMin<1.09)
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("rightSideband"));
          else
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("rightSideband"));
          model->paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(0.58812,0.911028,0.861955));
          xframe->getAttLine()->SetLineWidth(0);
          xframe->getAttText()->SetTextFont(44);
          xframe->getAttText()->SetTextSize(16);
          xframe->remove("model",false);
          if (ptMin<1.09)
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("model"));
          else{
            modelMismatch->plotOn(xframe_full, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSidebandK"), RooFit::Range("full"));
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("model"));
            //modelMismatch->plotOn(xframe, RooFit::Name("modelMismatch"), RooFit::LineColor(kGreen), RooFit::NormRange("rightSidebandK"), RooFit::Range("model"));
            //modelMismatch->paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", xframe->chiSquare("modelMismatch", "dataNsigma"))), RooFit::Layout(0.58812,0.911028,0.861955));
          }
          xframe->GetYaxis()->SetMaxDigits(2);

          // background integral
          double bkgIntegral = 0;
          double bkgIntegral_val = 0;
          if (ptMin < 1.09){
            bkgIntegral = ((RooAbsPdf *)model->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
            bkgIntegral_val = (nBackground2->getVal()) * bkgIntegral;
          }
          else{
            bkgIntegral = ((RooAbsPdf *)model->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
            bkgIntegral_val = (nBackground2->getVal()+nBackground1->getVal()) * bkgIntegral;
          }

          double rawYield, rawYieldError, counts;
          if (binCounting)
          {
            // total counts
            counts = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", intersectionBinCenter, signalRightLimit));
            rawYield = counts - bkgIntegral_val; // signal=counts-error
            std::cout << "Counts: " << counts << ", Background: "<< bkgIntegral_val << std::endl;
            rawYieldError = TMath::Sqrt(counts); // counts=signal+bkg => correct error from poisson statistics
            if (binCountingNoFit)
            {
              rawYield = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", -0.7815 * kNSigma, 0.7815 * kNSigma)); // only for He4 in signal loss studies
            }
          }
          else
          {
            std::cout << "No fit to signal peak!" << std::endl;
          }

          // fill raw yield histogram
          fTOFrawYield.SetBinContent(fTOFrawYield.FindBin((ptMax + ptMin) / 2.), rawYield);
          fTOFrawYield.SetBinError(fTOFrawYield.FindBin((ptMax + ptMin) / 2.), rawYieldError);

          // fill significance histogram
          if (binCounting)
            fSignificance.SetBinContent(fSignificance.FindBin((ptMax + ptMin) / 2.), rawYield / rawYieldError);
          else
            fSignificance.SetBinError(fSignificance.FindBin((ptMax + ptMin) / 2.), 0.);

          // fill sigma histogram
          fSigma.SetBinContent(fSigma.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParameter(2));
          fSigma.SetBinError(fSigma.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParError(2));

          // fill mean histogram
          fMean.SetBinContent(fMean.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParameter(1));
          fMean.SetBinError(fMean.FindBin((ptMax + ptMin) / 2.), signalRegionFit.GetParError(1));
        
        }
        // save to png
        TCanvas canv("canv", "canv");
        canv.cd();
        TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.4, 1.0, 1.0, 0);
        pad1->SetFillColor(0);
        pad1->SetFrameBorderMode(0);
        pad1->SetBottomMargin(0.06);
        //pad1->SetLogy();
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
        //xframe->GetYaxis()->SetRangeUser(-1.e4, 2.e4);
        xframe->GetYaxis()->SetTitleSize(0.07);
        xframe->GetYaxis()->SetTitleOffset(0.72);
        tofSignalProjection->GetXaxis()->SetRangeUser(-0.5, 0.5);
        double peakMaximum = tofSignalProjection->GetBinContent(tofSignalProjection->GetMaximumBin());
        std::cout << "peakMaximum=" << peakMaximum << std::endl;
        TLine lsx(intersectionBinCenter, 0, intersectionBinCenter, peakMaximum);
        lsx.SetLineStyle(kDashed);
        lsx.Draw("same");
        TLine ldx(signalRightLimit, 0, signalRightLimit, peakMaximum*0.75);
        if(ptMin>1.19)ldx.SetY2(peakMaximum*0.6);
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
        xframe_full->Write();
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
        TLine lsx1(intersectionBinCenter, -1.e4, intersectionBinCenter, 2.e4);
        lsx1.SetLineStyle(kDashed);
        //lsx1.Draw("same");
        TLine ldx1(signalRightLimit, -1.e4, signalRightLimit, 2.e4);
        ldx1.SetLineStyle(kDashed);
        //ldx1.Draw("same");
        TLine zero(nSigmaLeft, 0, maxNsigma, 0);
        zero.SetLineStyle(kDashed);
        zero.Draw("same");
        pad2->Update();
        canv.Write(Form("%s_%s_%d_%d_cent_%.0f_%.0f_pt_%.2f_%.2f", kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], ptMin, ptMax));
        canv.Print(Form("%s/signal_extraction/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], ptMin, ptMax));
        tofSignalProjection->SetMarkerStyle(20);
        tofSignalProjection->SetMarkerSize(0.8);
        TCanvas canv1("c1", "c1");
        tofSignalProjection->Draw("pe");
        canv1.Print(Form("%s/signal_extraction_preliminary/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.pdf", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], ptMin, ptMax));
      } // end of loop on pt bins

      // set raw yield histogram style and write to file
      fTOFrawYield.SetTitle(TString::Format("%s raw yield, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fTOFrawYield.SetName(TString::Format("f%sTOFrawYield_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fTOFrawYield.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fTOFrawYield.GetYaxis()->SetTitle("#it{N_{raw}}");
      fTOFrawYield.SetMarkerStyle(20);
      fTOFrawYield.SetMarkerSize(0.8);
      fTOFrawYield.SetOption("PE");
      fTOFrawYield.Write();

      // set significance histogram style and write to file
      fSignificance.SetTitle(TString::Format("%s significance, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fSignificance.SetName(TString::Format("f%sSignificance_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fSignificance.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSignificance.GetYaxis()->SetTitle("S/#sqrt{N}");
      fSignificance.SetMarkerStyle(20);
      fSignificance.SetMarkerSize(0.8);
      fSignificance.SetOption("PE");
      fSignificance.Write();

      // set sigma histogram style and write to file
      fSigma.SetTitle(TString::Format("%s sigma, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fSigma.SetName(TString::Format("f%sSigma_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fSigma.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fSigma.GetYaxis()->SetTitle("#sigma (GeV^2/#it{c}^{4})");
      fSigma.SetMarkerStyle(20);
      fSigma.SetMarkerSize(0.8);
      fSigma.SetOption("PE");
      fSigma.Write();

      // set mean histogram style and write to file
      fMean.SetTitle(TString::Format("%s mean, %.0f-%.0f%%", kAntimatterMatterLabel[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fMean.SetName(TString::Format("f%sMean_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]));
      fMean.GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fMean.GetYaxis()->SetTitle("#mu (GeV^2/#it{c}^{4})");
      fMean.SetMarkerStyle(20);
      fMean.SetMarkerSize(0.8);
      fMean.SetOption("PE");
      fMean.Write();

      //delete model, mean, sigma, alphaL, alphaR;
    } // end of loop on centrality bin
  }   // end of loop on antimatter/matter
  outFile->Close();
} // end of macro Signal.cpp