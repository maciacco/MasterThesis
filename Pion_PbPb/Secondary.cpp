// Secondary.cpp
// This macro computes the fraction of primaries

#include <stdlib.h>
#include <string>

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3F.h>
#include <TString.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TFractionFitter.h>
#include <Fit/Fitter.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <TFitResult.h>
#include <RooFitResult.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using namespace utils;
using namespace pion;

bool use_uniform = false;

const double fitRange = 1.25;

void Secondary(const char *cutSettings = "", const double DCAxyCut=0.12, const char *inFileDatName = "AnalysisResults", const char *inFileMCName = "mc", const char *outFileName = "PrimaryPion", const bool use_roofit = false, const bool useAntiPionsAsPrimaries = false)
{
  // killing RooFit output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel = kError; // Suppressing warning outputs

  const int MAX_ITER = 2147483647;

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetTextFont(44);
  // make signal extraction plots directory
  system(Form("mkdir %s/primary_fraction", kPlotDir));

  // open files
  //TFile *inFileDat = TFile::Open(Form("%s/%s.root", kDataDir, inFileDatName));
  TFile *inFileDat = TFile::Open(Form("%s/%s_largeNsigma_pion.root", kDataDir, inFileDatName));
  TFile *inFileMC21l5 = TFile::Open(Form("%s/%s.root", kDataDir, "AnalysisResults_LHC21l5_full_largeDCA"/* inFileMCName */));
  TFile *inFileMC20e3a_1 = TFile::Open(Form("%s/%s.root", kDataDir, "LHC20e3a"));
  //TFile *inFileMC20e3a_2 = TFile::Open(Form("%s/%s_20e3a_runlist2_20210929.root", kDataDir, "mc"));
  //TFile *inFileMC1 = TFile::Open(Form("%s/%s.root", kDataDir, inFileMCName));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    double noSecMaterialThreshold = 0.f; // DO NOT USE SECONDARIES FROM MATERIAL

    // make plot subdirectory
    system(Form("mkdir %s/primary_fraction/%s_%s", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList(s)
    std::string listName_true = Form("nuclei_pion_mcTrue_%s", cutSettings);
    std::string listName = Form("nuclei_pion_%s", cutSettings);
    TTList *listData = (TTList *)inFileDat->Get(listName.data());
    TTList *listMc21l5 = (TTList *)inFileMC21l5->Get(listName_true.data());
    TTList *listMc20e3a_1 = (TTList *)inFileMC20e3a_1->Get(listName_true.data());
    //TTList *listMc20e3a_2 = (TTList *)inFileMC20e3a_2->Get(listName.data());

    // get histograms from files
    TH3F *fDCAdat = (TH3F *)listData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAprim, *fDCAprim1, *fDCAprim2, *fDCAprim3;
    if (useAntiPionsAsPrimaries)
    {
      fDCAprim = (TH3F *)listData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[0]));
    }
    else
    {
      fDCAprim1 = (TH3F *)listMc20e3a_1->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
      //fDCAprim2 = (TH3F *)listMc20e3a_2->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[0]));
      //fDCAprim3 = (TH3F *)listMc3->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[0]));
      fDCAprim = (TH3F *)fDCAprim1->Clone(fDCAprim1->GetName());
      //fDCAprim->Add(fDCAprim2);
      //fDCAprim->Add(fDCAprim3);
    }
    TH3F *fDCAsec1 = (TH3F *)listMc20e3a_1->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
    //TH3F *fDCAsec2 = (TH3F *)listMc20e3a_2->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
    //TH3F *fDCAsec3 = (TH3F *)listMc3->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsec = (TH3F *)fDCAsec1->Clone(fDCAsec1->GetName());
    //fDCAsec->Add(fDCAsec2);
    //fDCAsec->Add(fDCAsec3);
    TH3F *fDCAsecWD1 = (TH3F *)listMc20e3a_1->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
    //TH3F *fDCAsecWD2 = (TH3F *)listMc20e3a_2->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
    //TH3F *fDCAsecWD3 = (TH3F *)listMc3->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsecWD = (TH3F *)fDCAsecWD1->Clone(fDCAsecWD1->GetName());
    //fDCAsecWD->Add(fDCAsecWD2);
    //fDCAsecWD->Add(fDCAsecWD3);

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      TH1D fPrimaryFrac(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), kNPtBins, kPtBins);
      TH1D fSecondaryFrac(Form("f%sSecFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), kNPtBins, kPtBins);

      int nUsedPtBins = 21;

      for (int iPtBin = 7; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        //fPrimaryFrac.SetBinContent(iPtBin, 1.);

        double ptMin = fPrimaryFrac.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fPrimaryFrac.GetXaxis()->GetBinUpEdge(iPtBin);
        int pTbinsIndexMin = fDCAdat->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fDCAdat->GetYaxis()->FindBin(ptMax - 0.005);
        outFile->cd();

        // project TH3 histogram
        TH1D *fDCAdatProj;
        TH1D *fDCAMcProjPrim;
        TH1D *fDCAMcProjSec;
        TH1D *fDCAMcProjSecWD;
        TString canvTitleTOF;
        TString canvNameTOF;
        TH1D fRatioDCAPrim("", "", kNDCABinsMediumOld, kDCABinsMediumOld);
        TH1D fRatioDCASec("", "", kNDCABinsMediumOld, kDCABinsMediumOld);

        TString projTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]));
        fDCAdatProj = fDCAdat->ProjectionZ(TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[iCent][0], kCentBinsPion[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAdatProj->SetTitle(projTitle);
        fDCAMcProjPrim = fDCAprim->ProjectionZ(TString::Format("f%sDCAPrimaryTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[iCent][0], kCentBinsPion[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjPrim->SetTitle(projTitle);
        fDCAMcProjSec = fDCAsec->ProjectionZ(TString::Format("f%sDCASecondaryTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[3][0], kCentBinsPion[3][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSec->SetTitle(projTitle);
        fDCAMcProjSecWD = fDCAsecWD->ProjectionZ(TString::Format("f%sDCASecondaryWeakTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[3][0], kCentBinsPion[3][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSecWD->SetTitle(projTitle);

        // rebin
        fDCAdatProj = (TH1D *)fDCAdatProj->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        fDCAMcProjPrim = (TH1D *)fDCAMcProjPrim->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        fDCAMcProjSec = (TH1D *)fDCAMcProjSec->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        fDCAMcProjSecWD = (TH1D *)fDCAMcProjSecWD->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        /* 
        fDCAdatProj = (TH1D*)fDCAdatProj->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
        fDCAMcProjPrim = (TH1D*)fDCAMcProjPrim->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
        fDCAMcProjSec = (TH1D*)fDCAMcProjSec->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
        fDCAMcProjSecWD = (TH1D*)fDCAMcProjSecWD->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
 */
        canvTitleTOF = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]));
        canvNameTOF = TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax));
        TCanvas canv(canvNameTOF, canvTitleTOF);

        fDCAdatProj->GetYaxis()->SetTitle("Entries");
        fDCAdatProj->Write();
        fDCAMcProjPrim->GetYaxis()->SetTitle("Entries");
        fDCAMcProjPrim->Write();
        fDCAMcProjSec->GetYaxis()->SetTitle("Entries");
        fDCAMcProjSec->Write();
        fDCAMcProjSecWD->GetYaxis()->SetTitle("Entries");
        fDCAMcProjSecWD->Write();

        // fit data distribution with TFractionFitter
        int size = 3;
        if (ptMin > noSecMaterialThreshold)
          size = 2;
        TObjArray *mc = new TObjArray(size); // MC histograms are put in this array
        mc->Add(fDCAMcProjPrim);
        mc->Add(fDCAMcProjSecWD);
        if (ptMin < noSecMaterialThreshold)
          mc->Add(fDCAMcProjSec);

        // check with constant
        if (use_uniform)
        {
          TH1D *hWeight = (TH1D *)fDCAdatProj->Clone();
          hWeight->Clear();
          for (int i = 1; i < hWeight->GetNbinsX() + 1; ++i)
            hWeight->SetBinContent(i, 1.);
          mc->Pop();
          mc->Add(hWeight);
        }

        if (ptMin > noSecMaterialThreshold)
          mc->Pop(); // no secondaries

        // check primary fraction with roofit
        double prim_frac_roofit = 1., prim_frac_roofit_err = 0.;
        // ROOFIT UNIMPLEMENTED
        if (use_roofit)
        {
          // std::cout << "No RooFit implementation yet!" << std::endl;
          // return;
          RooRealVar *dca = new RooRealVar("DCA_{xy}", "DCAxy", -1.3, 1.3, "cm");
          RooDataHist *data = new RooDataHist("data", "data", *dca, fDCAdatProj);
          RooDataHist *prim_tmp = new RooDataHist("prim_tmp", "prim_tmp", *dca, fDCAMcProjPrim);
          RooHistPdf *prim = new RooHistPdf("prim", "prim", *dca, *prim_tmp);
          /* RooDataHist *sec_tmp = new RooDataHist("sec_tmp", "sec_tmp", *dca, fDCAMcProjSec);
          RooHistPdf *sec = new RooHistPdf("sec", "sec", *dca, *sec_tmp); */
          RooDataHist *sec_wd_tmp = new RooDataHist("sec_wd_tmp", "sec_wd_tmp", *dca, fDCAMcProjSecWD);
          RooHistPdf *sec_wd = new RooHistPdf("sec_wd", "sec_wd", *dca, *sec_wd_tmp);
          RooRealVar *primfrac;
          /* if (iMatt == 1) primfrac = new RooRealVar("#it{f_{prim}}","primfrac",0.7,0.6,1.);
          else */ primfrac = new RooRealVar("#it{f_{prim}}","primfrac",0.,1.);
          RooRealVar *sec_wd_frac = new RooRealVar("#it{f_{sec_wd}}","sec_wd_frac",0.0,0.7);
          RooAddPdf *model;
          /* if (iMatt == 1 && ptMin < noSecMaterialThreshold)
            model = new RooAddPdf("model", "model", RooArgList(*prim, *sec_wd, *sec), RooArgList(*primfrac, *sec_wd_frac));
          else */
            model = new RooAddPdf("model", "model", RooArgList(*prim, *sec_wd), RooArgList(*primfrac));

          RooFitResult* res = model->fitTo(*data,RooFit::Save());

          // integrate primaries (-0.12,0.12)
          dca->setRange("intRange",-DCAxyCut,DCAxyCut);
          Double_t prim_integral = (prim->createIntegral(*dca,RooFit::NormSet(*dca),RooFit::Range("intRange")))->getVal();
          //std::cout << setprecision(8) << prim_integral << std::endl;
          // integrate model
          Double_t tot_integral = (model->createIntegral(*dca,RooFit::NormSet(*dca),RooFit::Range("intRange")))->getVal();
          Double_t ratio = primfrac->getVal()*prim_integral/tot_integral;
          Double_t ratio_err = primfrac->getError()*prim_integral/tot_integral;//TMath::Sqrt(ratio * (1. - ratio) / (tot_integral*fDCAdatProj->Integral()));
          std::cout << "Prim frac error = " << primfrac->getError()*prim_integral/tot_integral << std::endl;

          if (ratio < 0.7){continue;
            primfrac->setVal(0.85);
            primfrac->setMin(0.65);
            primfrac->setMax(0.9);
            res = model->fitTo(*data,RooFit::Save());

            dca->setRange("intRange",-DCAxyCut,DCAxyCut);
            prim_integral = (prim->createIntegral(*dca,RooFit::NormSet(*dca),RooFit::Range("intRange")))->getVal();
            //std::cout << setprecision(8) << prim_integral << std::endl;
            // integrate model
            tot_integral = (model->createIntegral(*dca,RooFit::NormSet(*dca),RooFit::Range("intRange")))->getVal();
            ratio = primfrac->getVal()*prim_integral/tot_integral;
            ratio_err = primfrac->getError()*prim_integral/tot_integral;
          }
          if (res->status() != 0) continue;
          RooPlot *xframe = dca->frame(RooFit::Title("dca"), RooFit::Name(fDCAdatProj->GetName()));
          data->plotOn(xframe,RooFit::Name("data"));
          model->plotOn(xframe, RooFit::Components("prim"), RooFit::LineColor(kBlue));
          model->plotOn(xframe, RooFit::Components("sec_wd"), RooFit::LineColor(kOrange));
          /* if (iMatt == 1 && ptMin < noSecMaterialThreshold)
            model->plotOn(xframe, RooFit::Components("sec"), RooFit::LineColor(kRed)); */
          model->plotOn(xframe, RooFit::Name("model"), RooFit::LineColor(kGreen));
          std::cout << "Chi2 = " << xframe->chiSquare("model", "data") << std::endl;
          xframe->Write();

          TCanvas c(xframe->GetName(),xframe->GetTitle());
          double chi2 = xframe->chiSquare("model", "data");
          xframe->Draw();
          TLatex chiSq(0.65, 1.e7, Form("#chi^{2}/NDF=%.2f", chi2));
          chiSq.SetNDC();
          chiSq.SetY(0.55);
          chiSq.SetTextSize(22);
          chiSq.Draw("same");
          c.Write();

          // save primary fraction and its uncertainty
          prim_frac_roofit = ratio;
          prim_frac_roofit_err = ratio_err;

          // print fit results
          std::cout << "f_prim_fit = " << primfrac->getVal() << std::endl;
          std::cout << "fraction = " << primfrac->getVal()*prim_integral/tot_integral << std::endl;
        }

        TFractionFitter *fit = new TFractionFitter(fDCAdatProj, mc, "Q"); // initialise
        ROOT::Fit::Fitter *fitter = fit->GetFitter();

        // compute wd template fraction
        double dataIntegralDCAcut = fDCAdatProj->Integral(fDCAdatProj->FindBin(-DCAxyCut), fDCAdatProj->FindBin(DCAxyCut-0.001));
        double dataIntegral = fDCAdatProj->Integral();

        fit->SetRangeX(fDCAdatProj->FindBin(-fitRange), fDCAdatProj->FindBin(fitRange-0.001));
        ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
        if (ptMin < noSecMaterialThreshold)
        { 
          fit->Constrain(2, 0., 0.06);
        }
        if (iCent < 2)
          fit->Constrain(1, 0., 0.9);
        else if (iCent==2 && ptMin < 0.9)
          fit->Constrain(0, 0., 1.);

        TVirtualFitter::SetMaxIterations(MAX_ITER);    
        /* TVirtualFitter::SetPrecision(1e-2);  */
        /* TFitResultPtr  */ Int_t status = fit->Fit(); // perform the fit
        if (status != 0)
        {
          fit = new TFractionFitter(fDCAdatProj, mc, "Q"); // initialise
          fitter = fit->GetFitter();
          ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);

          TVirtualFitter::SetMaxIterations(MAX_ITER);
         // TVirtualFitter::SetPrecision(0.01);

          fit->SetRangeX(fDCAdatProj->FindBin(-fitRange), fDCAdatProj->FindBin(fitRange-0.001));
          fitter->Config().ParSettings(0).SetValue(0.710);
          fitter->Config().ParSettings(0).SetLimits(0.600, 0.990);
          fitter->Config().ParSettings(0).Release();
          fitter->Config().ParSettings(0).SetStepSize(1.e-3);
          fitter->Config().ParSettings(1).SetValue(0.300);
          fitter->Config().ParSettings(1).SetLimits(0.010, 0.400);
          fitter->Config().ParSettings(1).Release();
          fitter->Config().ParSettings(1).SetStepSize(1.e-3);
          if (ptMin < noSecMaterialThreshold)
          {
            fitter->Config().ParSettings(2).SetValue(0.0050);
            fitter->Config().ParSettings(2).SetLimits(0.0010, 0.0250);
            fitter->Config().ParSettings(2).Release();
            fitter->Config().ParSettings(2).SetStepSize(1e-5);
          }
          ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
          status = fit->Fit();
        }
        //TFitResultPtr result = fit->Fit(); // perform the fit
        if (status == 0 /* || status == 4 */ /* || (result->CovMatrixStatus() == 1 || result->CovMatrixStatus() == 3) */ /* 0==0 */)
        { // check on fit status
          TH1F *result = (TH1F *)fit->GetPlot();
          TH1F *mc1 = (TH1F *)fit->GetMCPrediction(0);
          TH1F *mc2 = (TH1F *)fit->GetMCPrediction(1);
          TH1F *mc3;
          mc1->SetName(Form("%s_Prediction", fDCAMcProjPrim->GetName()));
          mc2->SetName(Form("%s_Prediction", fDCAMcProjSec->GetName()));
          double integralData = fDCAdatProj->Integral();
          double fracMc1, fracMc2, fracMc3, errFracMc1, errFracMc2, errFracMc3;
          fit->GetResult(0, fracMc1, errFracMc1);
          fit->GetResult(1, fracMc2, errFracMc2);
          mc1->Scale(fracMc1 * integralData / mc1->Integral(), "width");
          mc2->Scale(fracMc2 * integralData / mc2->Integral(), "width");

          // set color
          fDCAdatProj->SetMarkerStyle(20);
          fDCAdatProj->SetMarkerSize(0.8);
          fDCAdatProj->SetMarkerColor(kBlack);
          result->SetLineColor(kGreen + 2);
          result->SetLineWidth(3);
          mc1->SetLineColor(kBlue);
          mc2->SetLineColor(kOrange);

          TLegend leg(0.574499 + 0.05, 0.60552, 0.879628 + 0.05, 0.866667);
          leg.AddEntry(fDCAdatProj, "data");
          leg.AddEntry(result, "fit");
          leg.AddEntry(mc1, "primary pions");
          leg.AddEntry(mc2, "sec. weak pions");
          leg.SetTextSize(0.035);
          leg.SetBorderSize(0);

          // draw on canvas
          double chi2 = fit->GetChisquare();
          gStyle->SetOptStat(0);
          fDCAdatProj->Scale(1, "width");
          fDCAdatProj->GetXaxis()->SetRangeUser(-fitRange,fitRange);
          result->Scale(1, "width");
          fDCAdatProj->Draw("Ep");
          result->Draw("histosame");
          fDCAdatProj->Draw("Epsame");

          if (ptMin < noSecMaterialThreshold)
          {
            mc3 = (TH1F *)fit->GetMCPrediction(2);
            mc3->SetName(Form("%s_Prediction", fDCAMcProjSecWD->GetName()));
            fit->GetResult(2, fracMc3, errFracMc3);
            mc3->Scale(fracMc3 * integralData / mc3->Integral(), "width");
            mc3->SetLineColor(kRed);
            leg.AddEntry(mc3, "secondary pions");
            mc3->Draw("histosame");
          }

          mc2->Draw("histosame");
          mc1->Draw("histosame");

          leg.Draw("same");
          TLatex chiSq(0.65, 1.e7, Form("#chi^{2}/NDF=%.2f", chi2 / fit->GetNDF()));
          chiSq.SetNDC();
          chiSq.SetY(0.55);
          TLatex covStatus(-1., 0.1 * result->GetMaximum(), Form("f_{prim} = %.2f", fracMc1 /* "covStatus = %d", status->CovMatrixStatus() */));
          chiSq.SetTextSize(22);
          covStatus.SetTextSize(22);
          chiSq.Draw("same");
          covStatus.Draw("same");

          // print fit results
          std::cout << "f_prim_TFFfit = " << fracMc1 << " +/- " << errFracMc1 << std::endl;

          // compute fraction of primaries and material secondaries
          double intPrimDCAcutError = 0.;
          double intPrimDCAcut = mc1->IntegralAndError(result->FindBin(-DCAxyCut), result->FindBin(DCAxyCut-0.001), intPrimDCAcutError);
          //double intSecDCAcut = mc3->Integral(result->FindBin(-0.12), result->FindBin(0.115));
          double intResDCAcutError = 0.;
          double intResDCAcut = result->IntegralAndError(result->FindBin(-DCAxyCut), result->FindBin(DCAxyCut-0.001), intResDCAcutError);
          double primaryRatio = intPrimDCAcut / intResDCAcut;
          double primaryRatioError = errFracMc1/fracMc1*primaryRatio;//primaryRatio * TMath::Sqrt(intPrimDCAcutError * intPrimDCAcutError / intPrimDCAcut / intPrimDCAcut + intResDCAcutError * intResDCAcutError / intResDCAcut / intResDCAcut); // TMath::Sqrt(primaryRatio * (1.f - primaryRatio) / intResDCAcut);
          //double primaryRatioError = TMath::Sqrt(primaryRatio * (1.f - primaryRatio) / intResDCAcut);
          if (primaryRatio < 1.e-7)
            primaryRatioError = 1. / intResDCAcut;
          //double secondaryRatio = intSecDCAcut / intResDCAcut;
          //double secondaryRatioError = TMath::Sqrt(secondaryRatio * (1.f - secondaryRatio) / intResDCAcut);
          if (primaryRatio>1.)continue;
          fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatio);
          fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatioError);

          std::cout << "fraction = " << intPrimDCAcut / intResDCAcut << std::endl;
          std::cout << " * * * * * * * * * * * * * * * * * " << std::endl;
          // roofit check
          if (use_roofit)
          {
            fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(ptMin + 0.005f), prim_frac_roofit);
            fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(ptMin + 0.005f), prim_frac_roofit_err);
          }

          /* if (ptMin < noSecMaterialThreshold) {
            fSecondaryFrac.SetBinContent(fSecondaryFrac.FindBin(ptMin + 0.005f), secondaryRatio);
            fSecondaryFrac.SetBinError(fSecondaryFrac.FindBin(ptMin + 0.005f), secondaryRatioError);
          } */

          canv.SetLogy();
          canv.Write();

          // write fitted histograms
          mc1->Write();
          mc2->Write();

          // scale starting histograms
          fDCAMcProjPrim->Scale(fracMc1 * integralData / fDCAMcProjPrim->Integral(), "width");

          // compare template with predictions from fit
          fRatioDCAPrim.SetName(Form("f%sRatioDCAPrim_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)));
          fRatioDCAPrim.SetTitle(fRatioDCAPrim.GetName());
          fDCAMcProjSecWD->Scale(fracMc2 * integralData / fDCAMcProjSecWD->Integral(), "width");

          if (ptMin < noSecMaterialThreshold)
          {
            fDCAMcProjSec->Scale(fracMc3 * integralData / fDCAMcProjSec->Integral(), "width");
          }

          fRatioDCASec.SetName(Form("f%sRatioDCASec_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)));
          fRatioDCASec.SetTitle(fRatioDCASec.GetName());

          auto kNDCABins = kNDCABinsMediumOld;
          for (int iDCA = 1; iDCA < kNDCABins + 1; ++iDCA)
          {
            double primPrediction = mc1->GetBinContent(iDCA);
            double primMc = fDCAMcProjPrim->GetBinContent(iDCA);
            double primPrediction_err = mc1->GetBinError(iDCA);
            double primMc_err = fDCAMcProjPrim->GetBinError(iDCA);
            (primMc > 1.e-1) ? fRatioDCAPrim.SetBinContent(iDCA, primPrediction / primMc) : fRatioDCAPrim.SetBinContent(iDCA, 0);
            (primMc > 1.e-1) ? fRatioDCAPrim.SetBinError(iDCA, TMath::Sqrt(primPrediction_err * primPrediction_err / primPrediction / primPrediction + primMc_err * primMc_err / primMc / primMc)) : fRatioDCAPrim.SetBinError(iDCA, 0.);

            if (ptMin < 10/* noSecMaterialThreshold */)
            {
              double secPrediction = mc2->GetBinContent(iDCA);
              double secMc = fDCAMcProjSecWD->GetBinContent(iDCA);
              double secPrediction_err = mc2->GetBinError(iDCA);
              double secMc_err = fDCAMcProjSec->GetBinError(iDCA);
              (secMc > 1.e-1) ? fRatioDCASec.SetBinContent(iDCA, secPrediction / secMc) : fRatioDCASec.SetBinContent(iDCA, 0);
              (secMc > 1.e-1) ? fRatioDCASec.SetBinError(iDCA, TMath::Sqrt(secPrediction_err * secPrediction_err / secPrediction / secPrediction + secMc_err * secMc_err / secMc / secMc)) : fRatioDCASec.SetBinError(iDCA, 0.);
            }
          }
          fRatioDCAPrim.GetXaxis()->SetTitle(kAxisTitleDCA);
          fRatioDCAPrim.GetYaxis()->SetTitle("Prediction / MC");
          fRatioDCAPrim.GetYaxis()->SetRangeUser(0., 8.);
          fRatioDCAPrim.Write();

          // write histograms to file
          fDCAdatProj->GetXaxis()->SetTitleSize(0.05);
          //fDCAdatProj->GetYaxis()->SetRangeUser(1e3,1e9);
          fDCAdatProj->Write();
          fDCAMcProjPrim->Write();
          if (/* ptMin < noSecMaterialThreshold */ptMin < 10)
          {
            fRatioDCASec.GetXaxis()->SetTitle(kAxisTitleDCA);
            fRatioDCASec.GetXaxis()->SetTitleSize(0.05);
            fRatioDCASec.GetYaxis()->SetTitle("Prediction / MC");
            fRatioDCASec.GetYaxis()->SetRangeUser(0., 8.);
            fRatioDCASec.Write();
            fDCAMcProjSec->Write();
          }

          // save canvas plot
          canv.Print(Form("%s/primary_fraction/%s_%s/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1], ptMin, ptMax));
        }
      } // end of loop on centrality bin

      // primary fraction fit with fFitFunc function
      gStyle->SetStatX(0.85);
      gStyle->SetStatY(0.5);
      gStyle->SetStatFontSize(0.035);
      TF1 fFitFunc(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), "1/(1+[0]*exp([1]*x))", 0.5f, 3.0f);
      fFitFunc.SetParLimits(0, 0., 10000.);
      fFitFunc.SetParLimits(1, -100., 0.);
      fPrimaryFrac.Fit(&fFitFunc, "MRQ");
      auto r = fPrimaryFrac.Fit(&fFitFunc, "MRQS");
      auto cov_mat = r->GetCovarianceMatrix();
      //corr_mat.Print();
      //std::cout<<"0,0 = "<<corr_mat(0,0)<<std::endl;
      //corr_mat.Print();
      TH2D hCovMat(Form("f%sCovMat_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]),"CorrelationMatrix",2,0,2,2,0,2);
      for (int i=0;i<2;++i){
        for(int j=0;j<2;++j){
          //std::cout << cov_mat(0,1) << std::endl;
          hCovMat.SetBinContent(i+1,j+1,cov_mat(i,j));
        }
      }
      hCovMat.Write();
      fFitFunc.Write();

      fPrimaryFrac.SetMarkerStyle(20);
      fPrimaryFrac.SetMarkerSize(0.8);
      fPrimaryFrac.SetOption("pe");
      fPrimaryFrac.GetYaxis()->SetTitle("#it{f}_{#it{prim}}");
      fPrimaryFrac.GetXaxis()->SetTitle(kAxisTitlePt);
      fPrimaryFrac.GetXaxis()->SetRangeUser(0.5, 5.0);
      fPrimaryFrac.Write();

      system(Form("mkdir %s/primary_plots", kPlotDir));
      TCanvas cPrim("cPrim", "cPrim");
      cPrim.cd();
      fPrimaryFrac.GetXaxis()->SetRangeUser(1., 2.0);
      fPrimaryFrac.GetYaxis()->SetRangeUser(0.0, 1.1);
      fPrimaryFrac.Draw("");
      cPrim.Print(Form("%s/primary_plots/%s.png", kPlotDir, fPrimaryFrac.GetName()));
    }
  }
  outFile->Close();
}