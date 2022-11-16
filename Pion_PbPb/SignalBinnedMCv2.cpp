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

void SignalBinnedMCv2(const char *cutSettings = "", const double roi_min_limit_input = 1.5, const double roi_max_limit_input = 11., const double mismatch_min_limit_input = 8.5, const double mismatch_max_limit_input = 13.5, const bool binCounting = false, const int bkg_shape = 1, const char *inFileDat = "mc_20g7_likeData_largeNsigma", const char *outFileName = "SignalProton", const char *outFileOption = "recreate", const bool extractSignal = true, const bool useDSCB = false, const bool binCountingNoFit = false)
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
  TDirectory *dirOutFile = outFile->mkdir(Form("%s_%d_%d", cutSettings, binCounting, bkg_shape));
  TFile *mcFile_21l5 = TFile::Open(TString::Format("%s/%s.root", kDataDir, inFileDat)); // open data TFile
  TFile *mcFile_20g7 = TFile::Open(TString::Format("%s/%s.root", kDataDir, "mc_20g7_likeData_largeNsigma")); // open data TFile

  if (!mcFile_21l5)
  {
    if (kVerbose) std::cout << "File not found!" << std::endl; // check data TFile opening
    return;
  }

  // get TTList
  std::string listName_21l5 = Form("nuclei_pion_mcFalse_%s", cutSettings);
  TTList *list_21l5 = (TTList *)mcFile_21l5->Get(listName_21l5.data());

  // merge antimatter + matter histograms
  std::string histNameA = Form("f%sTOFnSigma", kAntimatterMatter[0]);
  std::string histNameM = Form("f%sTOFnSigma", kAntimatterMatter[1]);
  TH3F *fTOFSignalA = (TH3F *)list_21l5->Get(histNameA.data());
  TH3F *fTOFSignalAll = (TH3F *)fTOFSignalA->Clone(fTOFSignalA->GetName());
  TH3F *fTOFSignalM = (TH3F *)list_21l5->Get(histNameM.data());
  fTOFSignalAll->Add(fTOFSignalM);

  for (int iMatt = 1; iMatt > -1; --iMatt)
  { // loop on antimatter/matter
    // Get histograms from file
    std::string histName = Form("f%sTOFnSigma", kAntimatterMatter[iMatt]);
    TH3F *fTOFSignal1 = (TH3F *)list_21l5->Get(histName.data());
    fTOFSignal1->SetName("A");

    if (!fTOFSignal1)
    {
      if (kVerbose) std::cout << "Hstogram not found!" << std::endl; // check data TFile opening
      return;
    }

    TH3F *fTOFSignal = (TH3F *)fTOFSignal1->Clone(histName.data());

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

      for (int iPtBin = 7; iPtBin < nUsedPtBins + 1 -1; ++iPtBin)
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
        if (kVerbose) std::cout << "Pt bins: min=" << ptMin << "; max=" << ptMax << "; indexMin=" << pTbinsIndexMin << "; indexMax=" << pTbinsIndexMax << "; centralityBinMin=" << centBinMin << "; centralityBinMax=" << centBinMax << std::endl;
        TH1D *tofSignalProjection = fTOFSignal->ProjectionZ(Form("f%sTOFSignal_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        tofSignalProjection->Rebin(2);

        // project histogram (antimatter + matter)
        TH1D *tofSignalProjectionAll = fTOFSignalAll->ProjectionZ(Form("f%sTOFSignalAll_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignalAll->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignalAll->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        tofSignalProjectionAll->Rebin(2);

        // limits
        double nSigmaLeft = -10;
        double nSigmaRight = -5.;
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(nSigmaLeft, kTOFnSigmaMax);
        tofSignalProjection->GetXaxis()->SetRangeUser(nSigmaLeft, kTOFnSigmaMax);

        // DEFINE SIGNAL REGION
        TF1 signalRegionFit("signalRegionFit", "gaus", -20., 20.);
        // signalRegionFit.SetParLimits(0, 1., 1.e7);
        // signalRegionFit.SetParLimits(1, -0.4, 0.4);
        // signalRegionFit.SetParLimits(2, 0.8, 1.6);
        signalRegionFit.SetLineColor(kBlue);
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-.5, .5);
        double maximum_signal = tofSignalProjectionAll->GetBinCenter(tofSignalProjectionAll->GetMaximumBin());
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-20., 20.);
        tofSignalProjectionAll->Fit("signalRegionFit", "QRL+", "", maximum_signal - 1., maximum_signal + 1.);
        tofSignalProjectionAll->Fit("signalRegionFit", "QRL+", "", signalRegionFit.GetParameter(1) - 1., signalRegionFit.GetParameter(1) + 1.);
        double mean_tmp = signalRegionFit.GetParameter(1);
        double rms_tmp = signalRegionFit.GetParameter(2);
        double roi_max_limit = mean_tmp+(roi_max_limit_input-2.)*rms_tmp; // smaller signal window in MC than in data
        if (ptMin>1.29)
          roi_max_limit = mean_tmp+(roi_max_limit_input-3.)*rms_tmp;
        if (ptMin>1.39)
          roi_max_limit = mean_tmp+(roi_max_limit_input-4.)*rms_tmp;

        // DEFINE K SIGNAL REGION -> MISMATCH FIT REGION
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-50., 50.);
        TF1 signalRegionFitK("signalRegionFitK", "gaus", -50., 50.);
        // signalRegionFitK.SetParLimits(0, 1., 1.e7);
        // signalRegionFitK.SetParLimits(1, 0., 20.);
        // signalRegionFitK.SetParLimits(2, 0.8, 3.);
        signalRegionFitK.SetLineColor(kBlue);
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(5., 50.);
        if (ptMin > 1.2)tofSignalProjectionAll->GetXaxis()->SetRangeUser(5., 15.);
        double maximum_signal_K = tofSignalProjectionAll->GetBinCenter(tofSignalProjectionAll->GetMaximumBin());
        tofSignalProjectionAll->Fit("signalRegionFitK", "QRL+", "", maximum_signal_K - 1., maximum_signal_K + 1.);
        tofSignalProjectionAll->Fit("signalRegionFitK", "QRL+", "", signalRegionFitK.GetParameter(1) - 2., signalRegionFitK.GetParameter(1) + 2.);
        double mean_tmp_K = signalRegionFitK.GetParameter(1);
        double rms_tmp_K = signalRegionFitK.GetParameter(2);
        double mismatch_left_limit = mean_tmp_K+(mismatch_min_limit_input+2.)*rms_tmp_K;
        double mismatch_right_limit = mean_tmp_K+(mismatch_max_limit_input+4.)*rms_tmp_K;
        if (iCent==4&&ptMin>1.29){
          mismatch_left_limit = mean_tmp_K+(mismatch_min_limit_input+3.)*rms_tmp_K;
          mismatch_right_limit = mean_tmp_K+(mismatch_max_limit_input+3.)*rms_tmp_K;
        }
        else if (iCent==2 || (iCent==4&&ptMin>0.89)){
          mismatch_left_limit = mean_tmp_K+(mismatch_min_limit_input+7.)*rms_tmp_K;
          if (iCent==4)mismatch_left_limit = mean_tmp_K+(mismatch_min_limit_input+6.)*rms_tmp_K;
          mismatch_right_limit = mean_tmp_K+(mismatch_max_limit_input+6.)*rms_tmp_K;
        }
        else if (iCent==4&&ptMin>0.79&&ptMin<0.89){
          mismatch_left_limit = mean_tmp_K+(mismatch_min_limit_input+4.)*rms_tmp_K;
          mismatch_right_limit = mean_tmp_K+(mismatch_max_limit_input+6.)*rms_tmp_K;
        }
        if (ptMin>1.39){
          mismatch_right_limit = mean_tmp_K+(mismatch_max_limit_input+1)*rms_tmp_K;
          if (iCent==4 && ptMin>1.49){
            mismatch_left_limit = 17;
            mismatch_right_limit = 21;
          }
          else if (iCent==4){
            mismatch_left_limit = 19;
            mismatch_right_limit = 24;
          }
        }
        if (mismatch_left_limit>50) mismatch_left_limit = 43.;
        if (mismatch_right_limit>50) mismatch_right_limit = 50.;
        std::cout<<"pT = "<<ptMin<<"mismatch range (right to K) = "<<mismatch_left_limit<<" - "<<mismatch_right_limit<<std::endl;
        tofSignalProjectionAll->GetXaxis()->SetRangeUser(-20., 20.);

        // roofit data
        double maxNsigma=20.;
        if ((ptMin>0.89&&iCent<4)||(ptMin>0.79&&iCent==4)){
          maxNsigma=mean_tmp_K-1.*rms_tmp_K;
        }
        RooRealVar tofSignal("tofSignal", "n#sigma_{p}", nSigmaLeft, maxNsigma, "a.u.");
        RooRealVar *tofSignal_full=new RooRealVar("tofSignal_full", "n#sigma_{p}full_", nSigmaLeft, 20., "a.u.");
        if ((ptMin>0.89&&iCent<4)||(ptMin>0.79&&iCent==4))tofSignal_full=new RooRealVar("tofSignal_full", "n#sigma_{p}full_", nSigmaLeft, mismatch_right_limit, "a.u.");
        RooRealVar *tofSignal_full_all=new RooRealVar("tofSignal_full_all", "n#sigma_{p}full_all_", nSigmaLeft, 20., "a.u.");
        if ((ptMin>0.89&&iCent<4)||(ptMin>0.79&&iCent==4))tofSignal_full_all=new RooRealVar("tofSignal_full_all", "n#sigma_{p}full_all_", nSigmaLeft, mismatch_right_limit, "a.u.");
        
        // roofit histogram
        RooDataHist data("data", "data", RooArgList(tofSignal), tofSignalProjection);
        RooDataHist data_full("data_full", "data_full", RooArgList(*tofSignal_full), tofSignalProjection);
        RooDataHist data_full_all("data_full_all", "data_full_all", RooArgList(*tofSignal_full), tofSignalProjectionAll);

        // roofit histogram
        RooDataHist dataAll("dataAll", "dataAll", RooArgList(tofSignal), tofSignalProjectionAll);
        if (kVerbose) std::cout << "Number of entries (roofit) = " << dataAll.sumEntries() << std::endl;
        if (kVerbose) std::cout << "Number of entries (root) = " << tofSignalProjectionAll->GetEntries() << std::endl;

        RooAbsPdf *background0;
        RooAbsPdf *background1;
        RooAbsPdf *background2;
        RooAbsPdf *background3;
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
        nBackground1 = new RooRealVar("#it{N}_{Bkg,1}", "nBackground1", 0.,10.);
        if (bkg_shape == 1 && ((iCent<4&&ptMin < 0.89)||(iCent==4&&ptMin<0.79)))
        { // expo
          background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
          modelAll = (RooAbsPdf*)new RooExponential("modelAll", "modelAll", tofSignal, *slope1);
        }
        else if (bkg_shape == 1 && !((iCent<4&&ptMin < 0.89)||(iCent==4&&ptMin<0.79)))
        { // double expo with fixed mismatch
          slope1 = new RooRealVar("#tau_{1}", "slope1", -0.1, -.7, 0.);
          slope2 = new RooRealVar("#tau_{2}", "slope2", 0., 10.);
          background0 = (RooAbsPdf *)new RooExponential("background1", "background1", *tofSignal_full, *slope1);
          background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
          background2 = (RooAbsPdf *)new RooExponential("background2", "background2", tofSignal, *slope2);
          background3 = (RooAbsPdf *)new RooExponential("background3", "background3", *tofSignal_full, *slope2);
          nBackground3 = new RooRealVar("#it{N}_{Bkg,3}", "nBackground3", 0.,1.e10);
          nBackground1 = new RooRealVar("#it{N}_{Bkg,1}", "nBackground1", 1., 0., 1.e5);
          modelMismatch = (RooAddPdf*)new RooAddPdf("model", "model", RooArgList(*background0), RooArgList(*nBackground1));
        }
        else
          if (kVerbose) std::cout << "!!!!!" << std::endl; // TODO: UPDATE FOLLOWING BLOCK

        int covQ = -999;
        double intersectionBinCenter=-20.;
        intersectionBinCenter=mean_tmp-(roi_min_limit_input+1.)*rms_tmp;
        // if (ptMin>0.9)intersectionBinCenter=mean_tmp-(roi_min_limit_input+0.5)*rms_tmp;
        // if (ptMin>1.09)intersectionBinCenter=mean_tmp-(roi_min_limit_input)*rms_tmp;
        // if (ptMin>1.29)intersectionBinCenter=mean_tmp-(roi_min_limit_input-0.5)*rms_tmp;
        double signalRightLimit=roi_max_limit;
        if (extractSignal)
        {

          if (((iCent<4&&ptMin < 0.89)||(iCent==4&&ptMin<0.79))){
            tofSignal.setRange("rightSideband", signalRightLimit, maxNsigma);
            
            modelAll->fitTo(dataAll, RooFit::Range("rightSideband"));
            r = modelAll->fitTo(dataAll, RooFit::Save(), RooFit::Range("rightSideband"));
            
            slope1->setConstant();
            nBackground1->setConstant();

            if (bkg_shape == 1 && ((iCent<4&&ptMin < 0.89)||(iCent==4&&ptMin<0.79)))
            { // expo
              nBackground2 = new RooRealVar("#it{N}_{Bkg}", "nBackground2", 0., 1.e9);
              model = new RooAddPdf("model", "model", RooArgList(*modelAll), RooArgList(*nBackground2));
            }
            model->fitTo(data, RooFit::Range("rightSideband"));
            r = model->fitTo(data, RooFit::Save(), RooFit::Range("rightSideband"));
          }
          else{
            // pion sideband
            tofSignal.setRange("rightSideband", roi_max_limit, maxNsigma);

            // K sideband
            tofSignal_full->setRange("rightSidebandK", mismatch_left_limit, mismatch_right_limit);

            // full range
            tofSignal_full->setRange("full", nSigmaLeft,mismatch_right_limit);

            // narrower range
            tofSignal_full->setRange("small", nSigmaLeft,maxNsigma);

            for (int I=0;I<2;++I)modelMismatch->fitTo(data_full_all, RooFit::Save(), RooFit::Range("rightSidebandK"));
            slope1->setConstant();
            for (int I=0;I<2;++I)modelMismatch->fitTo(data_full, RooFit::Save(), RooFit::Range("rightSidebandK"));
            nBackground1->setConstant();
            double mismatch_integral_=((RooAbsPdf *)modelMismatch->createIntegral(RooArgSet(*tofSignal_full), RooFit::NormSet(RooArgSet(*tofSignal_full)), RooFit::Range("full")))->getVal();
            double mismatch_integral_small=((RooAbsPdf *)modelMismatch->createIntegral(RooArgSet(*tofSignal_full), RooFit::NormSet(RooArgSet(*tofSignal_full)), RooFit::Range("small")))->getVal();
            if (kVerbose) std::cout << "mism_integral_full = " << mismatch_integral_ << "; mism_integral_small = " << mismatch_integral_small << std::endl;
            nBackground2 = new RooRealVar("#it{N}_{Bkg}", "nBackground", 0., 1.e7);
            nBackground2->setVal(nBackground1->getVal()*mismatch_integral_small);
            nBackground2->setConstant();
            model = new RooAddPdf("model", "model", RooArgList(*background1,*background2), RooArgList(*nBackground2,*nBackground3));
            for (int I=0;I<2;++I)r = model->fitTo(data, RooFit::Save(), RooFit::Range("rightSideband"));
            nBackground2->setConstant(false);
          }
          
          // find signal region
          double leftSignalLimit=-10;
          bool intersection=false;
          int iB=1;
          int binShiftIndex=tofSignalProjection->FindBin(leftSignalLimit);
          if (kVerbose) std::cout << "intersection bin center = " << intersectionBinCenter << std::endl;

          if (kVerbose) std::cout << "fit status: " << r->status() << ";" << std::endl;
          if (kVerbose) std::cout << "covariance quality: " << r->covQual() << std::endl;
          covQ = r->covQual();
          if (covQ<2.5 || r->status()!= 0) continue;
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
          if (((iCent<4&&ptMin < 0.89)||(iCent==4&&ptMin<0.79)))
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("rightSideband"));
          else
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("rightSideband"));
          model->paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(0.58812,0.911028,0.861955));
          xframe->getAttLine()->SetLineWidth(0);
          xframe->getAttText()->SetTextFont(44);
          xframe->getAttText()->SetTextSize(16);
          xframe->remove("model",false);
          if (((iCent<4&&ptMin < 0.89)||(iCent==4&&ptMin<0.79)))
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("model"));
          else{
            modelMismatch->plotOn(xframe_full, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSidebandK"), RooFit::Range("full"));
            model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("rightSideband"), RooFit::Range("model"));
          }
          xframe->GetYaxis()->SetMaxDigits(2);


          auto xframe1 = (RooPlot *)xframe->Clone();
          RooHist *hresid = xframe1->residHist();
          int binGtZero = -1.;
          int ibin=1;
          while(binGtZero<0){
            if (hresid->GetPointY(ibin)>0.){
              binGtZero=ibin;
            }
            else ibin++;
          }
          std::cout<<"============>>>>>> bin gt zero = "<<binGtZero<<std::endl;
          if (intersectionBinCenter<hresid->GetPointX(binGtZero)) intersectionBinCenter = hresid->GetPointX(binGtZero);
          tofSignal.setRange("signalRange", intersectionBinCenter, signalRightLimit);


          // background integral
          double bkgIntegral = 0;
          double bkgIntegral_1 = 0;
          double bkgIntegral_val = 0;
          if (((iCent<4&&ptMin < 0.89)||(iCent==4&&ptMin<0.79))){
            bkgIntegral = ((RooAbsPdf *)model->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
            bkgIntegral_val = (nBackground2->getVal()) * bkgIntegral;
          }
          else{
            bkgIntegral = (((RooAbsPdf *)(model->pdfList().at(0)))->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
            bkgIntegral_1 = (((RooAbsPdf *)(model->pdfList().at(1)))->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
            bkgIntegral_val = nBackground2->getVal() * bkgIntegral_1+nBackground1->getVal() * bkgIntegral;
          }

          double rawYield, rawYieldError, counts;
          if (binCounting)
          {
            // total counts
            counts = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", intersectionBinCenter, signalRightLimit));
            rawYield = counts - bkgIntegral_val; // signal=counts-error
            if (kVerbose) std::cout << "Counts: " << counts << ", Background: "<< bkgIntegral_val << std::endl;
            rawYieldError = TMath::Sqrt(counts); // counts=signal+bkg => correct error from poisson statistics
            if (binCountingNoFit)
            {
              rawYield = data.sumEntries(Form("tofSignal>%f && tofSignal<%f", -0.7815 * kNSigma, 0.7815 * kNSigma)); // only for He4 in signal loss studies
            }
          }
          else
          {
            if (kVerbose) std::cout << "No fit to signal peak!" << std::endl;
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
        if (kVerbose) std::cout << "peakMaximum=" << peakMaximum << std::endl;
        TLine lsx(intersectionBinCenter, 0, intersectionBinCenter, peakMaximum);
        lsx.SetLineStyle(kDashed);
        lsx.Draw("same");
        TLine ldx(signalRightLimit, 0, signalRightLimit, peakMaximum*0.75);
        //TLine ldx(mean_tmp_K-2.*rms_tmp_K, 0, mean_tmp_K-2.*rms_tmp_K, peakMaximum*0.75);
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

        tofSignalProjectionAll->Write();
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