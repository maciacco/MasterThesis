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

#include "../utils/Utils.h"
#include "../utils/Config.h"
#include "../utils/RooGausExp.h"
#include "../utils/RooDSCBShape.h"
#include "../utils/RooGausDExp.h"

using namespace utils;
using namespace proton;

double computeExponentialNormalisation(double,double,double);

const double kNSigma = 3; // define interval for bin counting

void SignalBinned(const char *cutSettings = "", const bool binCounting = false, const int bkg_shape = 1, const char *inFileDat = "AnalysisResults", const char *outFileName = "SignalProton", const char *outFileOption = "recreate", const bool extractSignal = true, const bool useDSCB = true, const bool binCountingNoFit = false)
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
  TFile *dataFile = TFile::Open(TString::Format("%s/%s_largeNsigma.root", kDataDir, inFileDat)); // open data TFile
  //TFile *dataFile2 = TFile::Open(TString::Format("%s/%s.root", kDataDir, inFileDat)); // open data TFile
  if (!dataFile)
  {
    std::cout << "File not found!" << std::endl; // check data TFile opening
    return;
  }

  // get TTList
  std::string listName = Form("nuclei_proton_%s", cutSettings);
  TTList *list = (TTList *)dataFile->Get(listName.data());
  //TTList *list2 = (TTList *)dataFile2->Get(listName.data());

  // merge antimatter + matter histograms
  std::string histNameA = Form("f%sTOFnSigma", kAntimatterMatter[0]);
  std::string histNameM = Form("f%sTOFnSigma", kAntimatterMatter[1]);
  TH3F *fTOFSignalA = (TH3F *)list->Get(histNameA.data());
  // TH3F *fTOFSignalA2 = (TH3F *)list2->Get(histNameA.data());
  TH3F *fTOFSignalAll = (TH3F *)fTOFSignalA->Clone(fTOFSignalA->GetName());
  TH3F *fTOFSignalM = (TH3F *)list->Get(histNameM.data());
  //TH3F *fTOFSignalM2 = (TH3F *)list2->Get(histNameM.data());
  //fTOFSignalAll->Add(fTOFSignalA2);
  fTOFSignalAll->Add(fTOFSignalM);
  //fTOFSignalAll->Add(fTOFSignalM2);

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  { // loop on antimatter/matter
    // Get histograms from file
    std::string histName = Form("f%sTOFnSigma", kAntimatterMatter[iMatt]);
    TH3F *fTOFSignal1 = (TH3F *)list->Get(histName.data());
    //TH3F *fTOFSignal2 = (TH3F *)list2->Get(histName.data());
    if (!fTOFSignal1)
    {
      std::cout << "Hstogram not found!" << std::endl; // check data TFile opening
      return;
    }

    TH3F *fTOFSignal = (TH3F *)fTOFSignal1->Clone(fTOFSignal1->GetName());
    //fTOFSignal->Add(fTOFSignal2);
  
    // make plot subdirectory
    system(Form("mkdir %s/signal_extraction/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));
    system(Form("mkdir %s/signal_extraction_preliminary/%s_%s_%d_%d", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape));

    dirOutFile->cd();
    for (int iCent = 0; iCent < kNCentClasses-2; ++iCent)
    { // loop over centrality classes
      TH1D fTOFrawYield("fRawYield", "fRawYield", kNPtBins, kPtBins);          // declare raw yields histogram
      TH1D fSignificance("fSignificance", "fSignificance", kNPtBins, kPtBins); // declare significance
      TH1D fSigma("fSigma", "fSigma", kNPtBins, kPtBins);
      TH1D fMean("fMean", "fMean", kNPtBins, kPtBins);
      TH1D fAlphaL("fAlphaL", "fAlphaL", kNPtBins, kPtBins);
      TH1D fAlphaR("fAlphaR", "fAlphaR", kNPtBins, kPtBins);
      int nUsedPtBins = 27;
      //int nUsedPtBins = 39;

      for (int iPtBin = 5; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        double ptMin = fTOFrawYield.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fTOFrawYield.GetXaxis()->GetBinUpEdge(iPtBin);
        double centMin = kCentBinsLimitsProton[iCent][0];
        double centMax = kCentBinsLimitsProton[iCent][1];

        int pTbinsIndexMin = fTOFSignal->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fTOFSignal->GetYaxis()->FindBin(ptMax - 0.005);
        int centBinMin = fTOFSignal->GetXaxis()->FindBin(kCentBinsLimitsProton[iCent][0]);//kCentBinsProton[iCent][0];
        int centBinMax = fTOFSignal->GetXaxis()->FindBin(kCentBinsLimitsProton[iCent][1]-1.);//kCentBinsProton[iCent][1];

        // project histogram
        std::cout << "Pt bins: min=" << ptMin << "; max=" << ptMax << "; indexMin=" << pTbinsIndexMin << "; indexMax=" << pTbinsIndexMax << "; centralityBinMin=" << centBinMin << "; centralityBinMax=" << centBinMax << std::endl;
        TH1D *tofSignalProjection = fTOFSignal->ProjectionZ(Form("f%sTOFSignal_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignal->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignal->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        tofSignalProjection->Rebin(2);

        // project histogram (antimatter + matter)
        TH1D *tofSignalProjectionAll = fTOFSignalAll->ProjectionZ(Form("f%sTOFSignalAll_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], centMin, centMax, fTOFSignalAll->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fTOFSignalAll->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), centBinMin, centBinMax, pTbinsIndexMin, pTbinsIndexMax);
        tofSignalProjectionAll->Rebin(2);
        
        // limits
        double nSigmaLeft = -20;
        if (ptMin > 1.56) nSigmaLeft = -16;
        // = mean_tmp+1.5*rms_tmp;
        // UNCOMMENT TO PRODUCE VARIABLE FIT RANGE
        /* if (ptMin < 2.01) nSigmaLeft = -12.;
        else if (ptMin < 2.41) nSigmaLeft = -8.;
        else if (ptMin < 2.81) nSigmaLeft = -6.;
        else nSigmaLeft = -5.; */

        // PRELIMINARY FIT :: MISMATCH
        // extract mismatch decay constant (all)
        double minNsigma = 15., maxNsigma = 20.;
        if (ptMin > 3.01) minNsigma = 13, maxNsigma = 16;
        TF1 expMismatch("expMismatch","expo",nSigmaLeft,kTOFnSigmaMax);
        tofSignalProjectionAll->Fit("expMismatch","QL+","",minNsigma,maxNsigma);
        expMismatch.SetLineColor(kRed);
        double expMismatchDecayConstant = expMismatch.GetParameter(1);
        double expMismatchPreExpFactor = expMismatch.GetParameter(0);
        double normRooFit = 5.e7;
        normRooFit = expMismatch.Integral(nSigmaLeft, kTOFnSigmaMax)/tofSignalProjectionAll->GetBinWidth(2);
        std::cout << "fit parameter = " << expMismatch.GetParameter(1) << "; normalisation = " << normRooFit << std::endl;

        // extract mismatch normalisation (split)
        TF1 expMismatch2("expMismatch2","expo",nSigmaLeft,kTOFnSigmaMax);
        expMismatch2.SetParLimits(0, 0., 50.);
        expMismatch2.FixParameter(1,expMismatchDecayConstant);
        tofSignalProjection->Fit("expMismatch2","QL+","",minNsigma,maxNsigma);
        expMismatch2.SetLineColor(kRed);
        double expMismatchDecayConstant2 = expMismatch2.GetParameter(1);
        double expMismatchPreExpFactor2 = expMismatch2.GetParameter(0);
        double normRooFit2 = 0.;
        normRooFit2 = expMismatch2.Integral(nSigmaLeft, kTOFnSigmaMax); normRooFit2/=tofSignalProjectionAll->GetBinWidth(2);
        std::cout << "binWidth = " << tofSignalProjectionAll->GetBinWidth(2) << "; fit parameter = " << expMismatch2.GetParameter(1) << "; normalisation = " << normRooFit2 << std::endl;
        // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // PRELIMINARY FIT :: POWER LAW BACKGROUND (K)
        // fit K tail (all)
        double minNsigmaTail = -20., maxNsigmaTail = -10.;
        if (ptMin > 1.56) minNsigmaTail = -16, maxNsigmaTail = -10;
        TF1 powerLawKtail("powerLawKtail","[0]*TMath::Power(x+[1],[2])+expo(3)",nSigmaLeft,kTOFnSigmaMax);
        powerLawKtail.SetParLimits(0,1.,1.e8);
        powerLawKtail.SetParameter(0,1.e7);
        powerLawKtail.SetParLimits(1,0.,50.);
        powerLawKtail.SetParameter(1,20.);
        powerLawKtail.SetParLimits(2,-10.,0.);
        powerLawKtail.FixParameter(3,expMismatch.GetParameter(0));
        powerLawKtail.FixParameter(4,expMismatch.GetParameter(1));
        powerLawKtail.SetLineColor(kGreen);
        tofSignalProjectionAll->Fit("powerLawKtail","QL+","",minNsigmaTail,maxNsigmaTail);
        double normRooFitTail = 5.e7;
        normRooFitTail = powerLawKtail.Integral(nSigmaLeft, kTOFnSigmaMax); normRooFitTail/=tofSignalProjectionAll->GetBinWidth(2);
        normRooFitTail-=normRooFit;

        // fit K tail (split)
        TF1 powerLawKtail2("powerLawKtail2","[0]*TMath::Power(x+[1],[2])+expo(3)",-20.,20.);
        powerLawKtail2.SetParLimits(0,1.,1.e8);
        powerLawKtail2.SetParameter(0,1.e7);
        powerLawKtail2.SetParLimits(1,0.,50.);
        powerLawKtail2.SetParameter(1,20.);
        powerLawKtail2.SetParLimits(2,-10.,0.);
        powerLawKtail2.SetParameter(2,-3.);
        powerLawKtail2.FixParameter(3,expMismatch2.GetParameter(0));
        powerLawKtail2.FixParameter(4,expMismatch2.GetParameter(1));
        powerLawKtail2.SetLineColor(kGreen);
        tofSignalProjection->Fit("powerLawKtail2","QL+","",-20.,-10.);

        // gaussian core fit
        /* tofSignalProjection->SetAxisRange(-0.5, 1.0);
        double maxProton = tofSignalProjection->GetMaximum(); */
        /* tofSignalProjection->SetAxisRange(-10.0, -3.0);
        double maxMesons = tofSignalProjection->GetBinCenter(tofSignalProjection->GetMaximumBin());
        std::cout<<maxMesons<<std::endl;
        tofSignalProjection->SetAxisRange(maxMesons-1.0, maxMesons+1.0);
        double mean_tmp = tofSignalProjection->GetMean();
        double rms_tmp = tofSignalProjection->GetRMS();
        tofSignalProjection->SetAxisRange(kTOFnSigmaMin, kTOFnSigmaMax); */

        /* 
        if (ptMin < 1.99) extendNsigmaLeft = 4.;
        else if (ptMin < 2.99) extendNsigmaLeft = 2.;
        else if (ptMin < 3.01) extendNsigmaLeft = 1.;
        else if (ptMin > 3.29) extendNsigmaLeft = -1.; */
        RooRealVar tofSignal("tofSignal", "n#sigma_{p}", nSigmaLeft, 20., "a.u.");

        // roofit histogram
        RooDataHist data("data", "data", RooArgList(tofSignal), tofSignalProjection);

        // roofit histogram
        RooDataHist dataAll("dataAll", "dataAll", RooArgList(tofSignal), tofSignalProjectionAll);
        std::cout<<"Number of entries (roofit) = "<<dataAll.sumEntries()<<std::endl;
        std::cout<<"Number of entries (root) = "<<tofSignalProjectionAll->GetEntries()<<std::endl;

        // build composite model
        RooRealVar mean("#mu", "mean", -0.1, 0.5, "a.u.");
        if (iCent > 0)mean.setRange(0., 0.6);
        RooRealVar *sigma;
        if (ptMin < 2.4) sigma = new RooRealVar("#sigma", "sigma", 0.9, 1.5, "a.u.");
        else sigma = new RooRealVar("#sigma", "sigma", 0.8, 1.3, "a.u.");
        RooRealVar *alphaL;
        /* if (iCent == 2) alphaL = new RooRealVar("#alpha_{L}", "alphaL", -3.0, -1.0);
        else  */alphaL = new RooRealVar("#alpha_{L}", "alphaL", -1.62, -1.8, -0.9);
        //if (ptMin > 2.51) alphaL->setConstant();
        RooRealVar *alphaR;
        if (ptMin < 2.01) alphaR = new RooRealVar("#alpha_{R}", "alphaR", 0.7, 1.5);
        else alphaR = new RooRealVar("#alpha_{R}", "alphaR", 0.6, 1.);
        RooRealVar a1("a1","a1",0.,10.);
        RooRealVar a2("a2","a2",0.,10.);
        RooRealVar n1("n1","n1", 0.,100.);
        RooRealVar n2("n2","n2", 0.,100.);

        // background crystal ball
        /* RooRealVar a11("a11","a11",1.,1.,2.);
        RooRealVar a21("a21","a21",1., 1.,2.);
        RooRealVar n11("n11","n11", 2.,10.);
        RooRealVar n21("n21","n21", 2.,10.); */

        RooAbsPdf *signal; if (!useDSCB)signal = new RooGausDExp("signal", "signal", tofSignal, mean, *sigma, *alphaL,  *alphaR);
        /* if (ptMin > 3.39)
        signal = new RooGausExp("signal", "signal", tofSignal, mean, *sigma, *alphaR);
        */
        if (useDSCB){
          signal = new RooDSCBShape("signal", "signal", tofSignal, mean, *sigma, a1, n1, a2, n2);
        }
        RooAbsPdf *background1;
        RooAbsPdf *background2;
        /* RooAbsPdf *background3; */
        /* RooAbsPdf *background; */
        /* RooRealVar *parameter; */
        RooRealVar *slope1;
        /* RooRealVar *parameter_a;
        RooRealVar *parameter_b; */
        RooRealVar *slope2;
        RooRealVar *mean2;
        RooRealVar *sigma2;
        /* RooRealVar *alpha2L = new RooRealVar("#alpha_{2,L}", "alpha2L", -3.0, -1.1);
        RooRealVar *alpha2R;
        alpha2R = new RooRealVar("#alpha_{2,R}", "alpha2R", 1.0, 0.9, 2.0); */
        /* RooRealVar *mean3;
        RooRealVar *sigma3; */
        RooRealVar *nBackground1;
        RooRealVar *nBackground2;

        RooAddPdf *model;
        RooRealVar nSignal("N_{sig}", "nSignal", 1., 1.e10);

        slope1 = new RooRealVar("#tau_{mismatch}", "slope1", expMismatchDecayConstant, -1., -0.0001);
        std::cout << "slope1 value = " << slope1->getVal() <<std::endl;
        nBackground1 = new RooRealVar("#it{N}_{Bkg,1}", "nBackground1", normRooFit, 1., 1.e9);
        std::cout << "nBackground1 value = " << nBackground1->getVal() <<std::endl;
        slope1->setConstant();
        nBackground1->setConstant();

        RooRealVar par0("par0","par0",powerLawKtail.GetParameter(1),0.,50.);
        par0.setConstant();
        RooRealVar par1("par1","par1",powerLawKtail.GetParameter(2),-10.,0.);
        par1.setConstant();
        if (bkg_shape == 1)
        { // expo

          /* if (ptMin < 1.51) // 2.41 if narrower range is considered
          {
            background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
            model = new RooAddPdf("model", "model", RooArgList(*signal, *background1), RooArgList(nSignal, *nBackground1));
          }
          else  */if (ptMin < 3.01)
          {
            background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
            background2 = (RooAbsPdf *)new RooGenericPdf("background2", "TMath::Power(@0+@1,@2)",RooArgList(tofSignal,par0,par1));
            nBackground2 = new RooRealVar("#it{N}_{Bkg,2}", "nBackground2",normRooFitTail, 1., 1.e9);
            //nBackground2->setConstant();
            model = new RooAddPdf("model", "model", RooArgList(*signal, *background1, *background2), RooArgList(nSignal, *nBackground1, *nBackground2));
          }
          else  // if (ptMin < 3.79)
          {
            //mean2 = new RooRealVar("#mu_{2}", "mu2", -11., -3.5);
            mean2 = new RooRealVar("#mu_{2}", "mu2", -8., -4.);
            sigma2 = new RooRealVar("#sigma_{2}", "sigma2", 0.5, 2.0);
            background2 = (RooAbsPdf *)new RooGaussian("background2", "background2", tofSignal, *mean2, *sigma2);// , *alpha2R );
            background1 = (RooAbsPdf *)new RooExponential("background1", "background1", tofSignal, *slope1);
            nBackground2 = new RooRealVar("#it{N}_{Bkg,2}", "nBackground2", 1., 1.e9);
            model = new RooAddPdf("model", "model", RooArgList(*signal, *background1, *background2), RooArgList(nSignal, *nBackground1, *nBackground2));
          }
        }
        else std::cout<<"!!!!!"<<std::endl; // TODO: UPDATE FOLLOWING BLOCK
        /* { // sum of expo + pol2
          parameter_a = new RooRealVar("parameter_a", "parameter_a", -0.2, -0.01);
          parameter_b = new RooRealVar("parameter_b", "parameter_b", 0., 0.05);
          if (ptMin < 2.41)
          {
            background = (RooAbsPdf *)new RooPolynomial("background", "background", tofSignal, RooArgList(*parameter_a, *parameter_b));
          }
          else if (ptMin < 3.04)
          {
            background1 = (RooAbsPdf *)new RooPolynomial("background1", "background1", tofSignal, RooArgList(*parameter_a, *parameter_b));
            if (iCent == 0 || iCent == 1)
            slope2 = new RooRealVar("#tau_{mesons}", "slope2", -5., -0.5);
            else
            slope2 = new RooRealVar("#tau_{mesons}", "slope2", -5., -0.5);
            background2 = (RooAbsPdf *)new RooExponential("background2", "background2", tofSignal, *slope2);nBackground2 = new RooRealVar("#it{f}_{Bkg,2}", "nBackground2", 0., 1.0);
            background = new RooAddPdf("background", "background", RooArgList(*background1, *background2), RooArgList(*nBackground2));
          }
          else // if (ptMin < 3.79) 
          {
            background3 = (RooAbsPdf *)new RooPolynomial("background3", "background3", tofSignal, RooArgList(*parameter_a, *parameter_b));
            //mean2 = new RooRealVar("#mu_{2}", "mu2", -11., -3.5);
            mean2 = new RooRealVar("#mu_{2}", "mu2", -8., -4.);
            sigma2 = new RooRealVar("#sigma_{2}", "sigma2", 0.5, 2.0);
            background1 = (RooAbsPdf *)new RooGaussian("background1", "background1", tofSignal, *mean2, *sigma2);// , *alpha2R );
            nBackground2 = new RooRealVar("#it{f}_{Bkg,2}", "nBackground2", 0., 1.0);

            background = new RooAddPdf("background", "background", RooArgList(*background1, *background3), RooArgList(*nBackground2));
          }
        } */
        
        
        int covQ = -999;

        if (extractSignal)
        {
          // fit model
          RooFitResult *r;

          //std::cout<<"ciao"<<std::endl;
          // fit antimatter + matter TOF signal distribution first
          for (int i = 0; i < 2; ++i)
          {
            r = model->fitTo(dataAll, RooFit::Save());
          }

          // fix parameters of the signal distribution
          /* nSignal.setRange(1., 1.e9);
          nBackground.setRange(1., 1.e7); */
          /* if (iCent == 0 || iCent == 2)
          alpha2R->setRange(1.0, 3.0);
          else
          alpha2R->setRange(1.0, 3.);
          if (bkg_shape == 1)
          { // expo
            if (ptMin < 2.61)
            {
              slope1->setRange(-0.2, -0.01);
              if (iCent == 0 || iCent == 1)
              slope2->setRange(-10., -1.5);
              else
              slope2->setRange(-5., -1.5);
            }
            else
            {
              slope1->setRange(-0.5, -0.01);
              mean2->setRange(-6., -2.0);
              sigma2->setRange(0.5, 2.0);
              nBackground2->setRange(0., 1.0);
            }
          } */
          mean.setConstant();
          sigma->setConstant();
          /* alphaR->setConstant();
          alphaL->setConstant(); */

          // fix mismatch background shape
          
          /* slope1->setConstant(false);
          slope1->setVal(expMismatchDecayConstant2);
          slope1->setConstant();
          nBackground1->setConstant(false);
          nBackground1->setVal(normRooFit2);
          nBackground1->setConstant(); */

          // fit TOF signal distribution
          for (int i = 0; i < 2; ++i)
          {
            //r = model->fitTo(data, RooFit::Save());
          }
          std::cout << "fit status: " << r->status() << ";" << std::endl;
          std::cout << "covariance quality: " << r->covQual() << std::endl;
          covQ = r->covQual();

          // signal range
          double signal_min = mean.getVal() - kNSigma * sigma->getVal(), signal_max = mean.getVal() + kNSigma * sigma->getVal();
          tofSignal.setRange("signalRange", signal_min, signal_max);

          // background integral
          auto components = model->getComponents();
          auto bkgPdf = (RooAbsPdf *)components->find("background1");
          double bkgIntegral = (((RooAbsPdf *)model->pdfList().at(1))->createIntegral(RooArgSet(tofSignal), RooFit::NormSet(RooArgSet(tofSignal)), RooFit::Range("signalRange")))->getVal();
          double bkgIntegral_val = nBackground1->getVal() * bkgIntegral;
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

          // fill alphaR histogram
          fAlphaR.SetBinContent(fAlphaR.FindBin((ptMax + ptMin) / 2.), alphaR->getVal());
          fAlphaR.SetBinError(fAlphaR.FindBin((ptMax + ptMin) / 2.), alphaR->getError());

          // fill alphaL histogram
          fAlphaL.SetBinContent(fAlphaL.FindBin((ptMax + ptMin) / 2.), alphaL->getVal());
          fAlphaL.SetBinError(fAlphaL.FindBin((ptMax + ptMin) / 2.), alphaL->getError());
        }

        // frame
        int nBins = tofSignalProjection->GetNbinsX();
        TString plotTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", ptMin, ptMax, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]);
        RooPlot *xframe = tofSignal.frame(RooFit::Bins(nBins), RooFit::Title(plotTitle), RooFit::Name(Form("f%sNSigma_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax)));
        dataAll.plotOn(xframe, RooFit::Name("dataNsigma"));
        tofSignal.setRange("full", nSigmaLeft, kTOFnSigmaMax);
        if (extractSignal)
        {
          if (background2)
            model->plotOn(xframe, RooFit::Components("background2"), RooFit::Name("background2"), RooFit::LineStyle(kDashed), RooFit::LineColor(kOrange));
          model->plotOn(xframe, RooFit::Components("background1"), RooFit::Name("background1"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
          model->plotOn(xframe, RooFit::Components("signal"), RooFit::Name("signal"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
          model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue));
          xframe->remove("model",false);
          xframe->remove("background1",false);
          if (background2)
            xframe->remove("background2",false);
          xframe->remove("signal",false);
          if (background2)
            model->plotOn(xframe, RooFit::Components("background2"), RooFit::Name("background2"), RooFit::LineStyle(kDashed), RooFit::LineColor(kOrange));
          model->plotOn(xframe, RooFit::Components("background1"), RooFit::Name("background1"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));
          model->plotOn(xframe, RooFit::Components("signal"), RooFit::Name("signal"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
          model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kBlue), RooFit::NormRange("full"));
          double x_min = 0.62;
          if(ptMax < 1.1) x_min = 0.12;
          model->paramOn(xframe, RooFit::Label(TString::Format("#chi^{2}/NDF = %2.4f", xframe->chiSquare("model", "dataNsigma"))), RooFit::Layout(x_min, 0.99, 0.88));
          xframe->getAttText()->SetTextSize(0.03);
          xframe->getAttLine()->SetLineWidth(0);
          xframe->getAttFill()->SetFillColor(kWhite);
          xframe->getAttFill()->SetFillStyle(4000);
          xframe->GetYaxis()->SetMaxDigits(2);
        }

        tofSignalProjectionAll->Write();
        xframe->Write();

        // save to png
        TCanvas canv;
        canv.SetName(plotTitle);
        TLatex covarianceQuality;
        covarianceQuality.SetNDC();
        covarianceQuality.SetX(0.55);
        covarianceQuality.SetY(0.4);
        covarianceQuality.SetTitle(Form("Covariance quality = %d", covQ));
        xframe->Draw("");
        //covarianceQuality.Draw("same");
        canv.SetLogy();
        canv.Print(Form("%s/signal_extraction/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));

        tofSignalProjection->SetMarkerStyle(20);
        tofSignalProjection->SetMarkerSize(0.8);
        tofSignalProjection->Draw("pe");
        canv.Print(Form("%s/signal_extraction_preliminary/%s_%s_%d_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1], ptMin, ptMax));

      } // end of loop on pt bins

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

double computeExponentialNormalisation (double slope, double lowLimit, double upLimit){
  double exp1 = TMath::Exp(slope*upLimit);
  double exp2 = TMath::Exp(slope*lowLimit);
  double norm = -1/slope*(exp2-exp1);
  return norm;
}
