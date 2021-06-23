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

#include "../utils/Utils.h"
#include "../utils/Config.h"

using namespace utils;
using namespace deuteron;

bool use_uniform = false;
bool use_roofit = false;

void Secondary(const char *cutSettings = "", const char *inFileDatName = "AnalysisResults", const char *inFileMCName = "mc", const char *outFileName = "PrimaryDeuteron")
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  // make signal extraction plots directory
  system(Form("mkdir %s/primary_fraction", kPlotDir));

  // open files
  TFile *inFileDat = TFile::Open(Form("%s/%s.root", kDataDir, inFileDatName));
  TFile *inFileMC = TFile::Open(Form("%s/%s.root", kDataDir, inFileMCName));
  TFile *inFileMCsec = TFile::Open(Form("%s/%s_sec.root", kDataDir, inFileMCName));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iMatt = 1; iMatt < 2; ++iMatt)
  {
    // make plot subdirectory
    system(Form("mkdir %s/primary_fraction/%s_%s", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList(s)
    std::string listName = Form("mpuccio_deuterons_%s", cutSettings);
    TTList *listData = (TTList *)inFileDat->Get(listName.data());
    TTList *listMc = (TTList *)inFileMC->Get(listName.data());
    TTList *listMcSec = (TTList *)inFileMCsec->Get(listName.data());

    // get histograms from files
    TH3F *fDCAdat = (TH3F *)listData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAprim = (TH3F *)listMc->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsec = (TH3F *)listMcSec->Get(Form("fMDCASecondaryTOF"));
    TH3F *fDCAsecHighPt = (TH3F *)listMc->Get(Form("fMDCASecondaryTOF"));

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      TH1D fPrimaryFrac(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]), kNPtBins, kPtBins);
      TH1D fSecondaryFrac(Form("f%sSecFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]), kNPtBins, kPtBins);

      int nUsedPtBins = 22;

      for (int iPtBin = 2; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        fPrimaryFrac.SetBinContent(iPtBin, 1.);

        if (iPtBin > 14)
          continue;
        double ptMin = fPrimaryFrac.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fPrimaryFrac.GetXaxis()->GetBinUpEdge(iPtBin);
        int pTbinsIndexMin = fDCAdat->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fDCAdat->GetYaxis()->FindBin(ptMax - 0.005);
        outFile->cd();

        // project TH3 histogram
        TH1D *fDCAdatProj;
        TH1D *fDCAMcProjPrim;
        TH1D *fDCAMcProjSec;
        TString canvTitleTOF;
        TString canvNameTOF;
        TH1D fRatioDCAPrim("", "", kNDCABinsMedium, kDCABinsMedium);
        TH1D fRatioDCASec("", "", kNDCABinsMedium, kDCABinsMedium);

        TString projTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsDeuteron[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsDeuteron[iCent][1]));
        fDCAdatProj = fDCAdat->ProjectionZ(TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsDeuteron[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsDeuteron[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsDeuteron[iCent][0], kCentBinsDeuteron[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAdatProj->SetTitle(projTitle);
        fDCAMcProjPrim = fDCAprim->ProjectionZ(TString::Format("f%sDCAPrimary_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsDeuteron[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsDeuteron[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsDeuteron[iCent][0], kCentBinsDeuteron[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjPrim->SetTitle(projTitle);
        if (iPtBin < 8)
          fDCAMcProjSec = fDCAsec->ProjectionZ(TString::Format("f%sDCASecondary_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsDeuteron[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsDeuteron[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsDeuteron[0][0], kCentBinsDeuteron[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        else
          continue;
        //fDCAMcProjSec = fDCAsecHighPt->ProjectionZ(TString::Format("f%sDCASecondary_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsDeuteron[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsDeuteron[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsDeuteron[0][0], kCentBinsDeuteron[iCent][1], pTbinsIndexMin, pTbinsIndexMax);

        if (iPtBin > 5){
          fDCAdatProj = (TH1D *)fDCAdatProj->Rebin(kNDCABinsLarge, fDCAdatProj->GetName(), kDCABinsLarge);
          fDCAMcProjPrim = (TH1D *)fDCAMcProjPrim->Rebin(kNDCABinsLarge, fDCAMcProjPrim->GetName(), kDCABinsLarge);
          fDCAMcProjSec = (TH1D *)fDCAMcProjSec->Rebin(kNDCABinsLarge, fDCAMcProjSec->GetName(), kDCABinsLarge);
          fRatioDCAPrim = *(TH1D *)fRatioDCAPrim.Rebin(kNDCABinsLarge, fRatioDCAPrim.GetName(), kDCABinsLarge);
          fRatioDCASec = *(TH1D *)fRatioDCASec.Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
        }

        else if (iPtBin > 3 && iPtBin < 7){
          fDCAdatProj = (TH1D *)fDCAdatProj->Rebin(kNDCABinsMedium2, fDCAdatProj->GetName(), kDCABinsMedium2);
          fDCAMcProjPrim = (TH1D *)fDCAMcProjPrim->Rebin(kNDCABinsMedium2, fDCAMcProjPrim->GetName(), kDCABinsMedium2);
          fDCAMcProjSec = (TH1D *)fDCAMcProjSec->Rebin(kNDCABinsMedium2, fDCAMcProjSec->GetName(), kDCABinsMedium2);
          fRatioDCAPrim = *(TH1D *)fRatioDCAPrim.Rebin(kNDCABinsMedium2, fRatioDCAPrim.GetName(), kDCABinsMedium2);
          fRatioDCASec = *(TH1D *)fRatioDCASec.Rebin(kNDCABinsMedium2, fRatioDCASec.GetName(), kDCABinsMedium2);
        }
        
        else{
          fDCAdatProj = (TH1D *)fDCAdatProj->Rebin(kNDCABinsMedium, fDCAdatProj->GetName(), kDCABinsMedium);
          fDCAMcProjPrim = (TH1D *)fDCAMcProjPrim->Rebin(kNDCABinsMedium, fDCAMcProjPrim->GetName(), kDCABinsMedium);
          fDCAMcProjSec = (TH1D *)fDCAMcProjSec->Rebin(kNDCABinsMedium, fDCAMcProjSec->GetName(), kDCABinsMedium);
          fRatioDCAPrim = *(TH1D *)fRatioDCAPrim.Rebin(kNDCABinsMedium, fRatioDCAPrim.GetName(), kDCABinsMedium);
          fRatioDCASec = *(TH1D *)fRatioDCASec.Rebin(kNDCABinsMedium, fRatioDCASec.GetName(), kDCABinsMedium);
        }

        fDCAMcProjSec->SetTitle(projTitle);
        canvTitleTOF = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsDeuteron[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsDeuteron[iCent][1]));
        canvNameTOF = TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsDeuteron[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsDeuteron[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax));
        TCanvas canv(canvNameTOF, canvTitleTOF);

        fDCAdatProj->GetYaxis()->SetTitle("Entries");
        fDCAMcProjPrim->GetYaxis()->SetTitle("Entries");
        fDCAMcProjSec->GetYaxis()->SetTitle("Entries");
        fDCAMcProjSec->Write();

        // fit data distribution with TFractionFitter
        TObjArray *mc = new TObjArray(2); // MC histograms are put in this array
        mc->Add(fDCAMcProjPrim);
        mc->Add(fDCAMcProjSec);

        // check with constant
        if (use_uniform) {
          TH1D *hWeight = (TH1D*)fDCAdatProj->Clone();
          hWeight->Clear();
          for(int i = 1; i < hWeight->GetNbinsX()+1; ++i)
            hWeight->SetBinContent(i, 1.);
          mc->Pop();
          mc->Add(hWeight);
        }

        // check primary fraction with roofit
        double prim_frac_roofit = 1., prim_frac_roofit_err = 0.;
        if (use_roofit) {
          RooRealVar *dca = new RooRealVar("DCA_{xy}", "DCAxy", -1.3, 1.3, "cm");
          RooDataHist *data = new RooDataHist("data", "data", *dca, fDCAdatProj);
          RooDataHist *prim_tmp = new RooDataHist("prim_tmp", "prim_tmp", *dca, fDCAMcProjPrim);
          RooHistPdf *prim = new RooHistPdf("prim", "prim", *dca, *prim_tmp);
          RooDataHist *sec_tmp = new RooDataHist("sec_tmp", "sec_tmp", *dca, fDCAMcProjSec);
          RooHistPdf *sec = new RooHistPdf("sec", "sec", *dca, *sec_tmp);
          RooRealVar *primfrac = new RooRealVar("#it{f_{prim}}","primfrac",0.,1.);
          RooAddPdf *model = new RooAddPdf("model", "model", RooArgList(*prim, *sec), RooArgList(*primfrac));
          model->fitTo(*data);
          RooPlot *xframe = dca->frame(RooFit::Title("dca"), RooFit::Name(fDCAdatProj->GetName()));
          data->plotOn(xframe);
          model->plotOn(xframe, RooFit::Components("prim"), RooFit::LineColor(kBlue));
          model->plotOn(xframe, RooFit::Components("sec"), RooFit::LineColor(kRed));
          model->plotOn(xframe, RooFit::LineColor(kGreen));
          xframe->Write();

          // integrate primaries (-0.12,0.12)
          dca->setRange("intRange",-0.12,0.12);
          Double_t prim_integral = (prim->createIntegral(*dca,RooFit::NormSet(*dca),RooFit::Range("intRange")))->getVal();
          //std::cout << setprecision(8) << prim_integral << std::endl;
          // integrate model
          Double_t tot_integral = (model->createIntegral(*dca,RooFit::NormSet(*dca),RooFit::Range("intRange")))->getVal();
          Double_t ratio = primfrac->getVal()*prim_integral/tot_integral;
          Double_t ratio_err = TMath::Sqrt(ratio * (1. - ratio) / (tot_integral*fDCAdatProj->Integral()));
          
          // save primary fraction and its uncertainty
          prim_frac_roofit = ratio;
          prim_frac_roofit_err = ratio_err;
        }

        TFractionFitter *fit = new TFractionFitter(fDCAdatProj, mc, "Q"); // initialise
        ROOT::Fit::Fitter *fitter = fit->GetFitter();

        // compute wd template fraction
        double dataIntegralDCAcut = fDCAdatProj->Integral(fDCAdatProj->FindBin(-0.1), fDCAdatProj->FindBin(0.09));
        double dataIntegral = fDCAdatProj->Integral();

        fit->SetRangeX(fDCAdatProj->FindBin(-1.3), fDCAdatProj->FindBin(1.29));
        fit->Constrain(0, 0.0, 1.0); // constrain fraction 1 to be between 0 and 1
        fit->Constrain(1, 0.0, 1.0);

        Int_t status = fit->Fit(); // perform the fit

        if (status == 0)
        { // check on fit status
          TH1F *result = (TH1F *)fit->GetPlot();
          TH1F *mc1 = (TH1F *)fit->GetMCPrediction(0);
          TH1F *mc2 = (TH1F *)fit->GetMCPrediction(1);
          mc1->SetName(Form("%s_Prediction", fDCAMcProjPrim->GetName()));
          mc2->SetName(Form("%s_Prediction", fDCAMcProjSec->GetName()));
          double integralData = fDCAdatProj->Integral();
          double fracMc1, fracMc2, errFracMc1, errFracMc2;
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
          mc2->SetLineColor(kRed);

          TLegend leg(0.574499, 0.60552, 0.879628, 0.866667);
          leg.AddEntry(fDCAdatProj, "data");
          leg.AddEntry(result, "fit");
          leg.AddEntry(mc1, "primary deuterons");
          leg.AddEntry(mc2, "secondary deuterons");
          leg.SetTextSize(0.035);
          leg.SetBorderSize(0);

          // draw on canvas
          double chi2 = fit->GetChisquare();
          gStyle->SetOptStat(0);
          fDCAdatProj->Scale(1, "width");
          result->Scale(1, "width");
          // fDCAdatProj->GetYaxis()->SetRangeUser(1.e3, 1.e6);
          fDCAdatProj->Draw("Ep");
          result->Draw("histosame");
          fDCAdatProj->Draw("Epsame");

          mc2->Draw("histosame");
          mc1->Draw("histosame");

          leg.Draw("same");
          TLatex chiSq(0.28, 5.e4, Form("#chi^{2}/NDF=%.2f", chi2 / fit->GetNDF()));
          TLatex prob(-1., 0.1 * result->GetMaximum(), Form("prob=%.7f", fit->GetProb()));
          chiSq.SetTextSize(0.035);
          prob.SetTextSize(0.035);
          chiSq.Draw("same");
          //prob.Draw("same");

          // compute fraction of primaries and material secondaries
          double intPrimDCAcutError = 0.;
          double intPrimDCAcut = mc1->IntegralAndError(result->FindBin(-0.12), result->FindBin(0.115), intPrimDCAcutError);
          double intSecDCAcut = mc2->Integral(result->FindBin(-0.12), result->FindBin(0.115));
          double intResDCAcutError = 0.;
          double intResDCAcut = result->IntegralAndError(result->FindBin(-0.12), result->FindBin(0.115), intResDCAcutError);
          double primaryRatio = intPrimDCAcut / intResDCAcut;
          double primaryRatioError = primaryRatio*TMath::Sqrt(intPrimDCAcutError*intPrimDCAcutError/intPrimDCAcut/intPrimDCAcut+intResDCAcutError*intResDCAcutError/intResDCAcut/intResDCAcut); // TMath::Sqrt(primaryRatio * (1.f - primaryRatio) / intResDCAcut);
          if (primaryRatio < 1.e-7)
            primaryRatioError = 1. / intResDCAcut;
          double secondaryRatio = intSecDCAcut / intResDCAcut;
          double secondaryRatioError = TMath::Sqrt(secondaryRatio * (1.f - secondaryRatio) / intResDCAcut);
          fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatio);
          fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatioError);
          
          // roofit check
          if (use_roofit) {
            fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(ptMin + 0.005f), prim_frac_roofit);
            fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(ptMin + 0.005f), prim_frac_roofit_err);
          }

          fSecondaryFrac.SetBinContent(fSecondaryFrac.FindBin(ptMin + 0.005f), secondaryRatio);
          fSecondaryFrac.SetBinError(fSecondaryFrac.FindBin(ptMin + 0.005f), secondaryRatioError);

          // canv.SetLogy();
          canv.Write();

          // write fitted histograms
          mc1->Write();
          mc2->Write();

          // scale starting histograms
          fDCAMcProjPrim->Scale(fracMc1 * integralData / fDCAMcProjPrim->Integral(), "width");
          fDCAMcProjSec->Scale(fracMc2 * integralData / fDCAMcProjSec->Integral(), "width");

          // compare template with predictions from fit
          fRatioDCAPrim.SetName(Form("f%sRatioDCAPrim_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)));
          fRatioDCAPrim.SetTitle(fRatioDCAPrim.GetName());
          fRatioDCASec.SetName(Form("f%sRatioDCASec_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)));
          fRatioDCASec.SetTitle(fRatioDCASec.GetName());
          for (int iDCA = 1; iDCA < kNDCABins + 1; ++iDCA)
          {
            double primPrediction = mc1->GetBinContent(iDCA);
            double primMc = fDCAMcProjPrim->GetBinContent(iDCA);
            double primPrediction_err = mc1->GetBinError(iDCA);
            double primMc_err = fDCAMcProjPrim->GetBinError(iDCA);
            (primMc > 1.e-1) ? fRatioDCAPrim.SetBinContent(iDCA, primPrediction / primMc) : fRatioDCAPrim.SetBinContent(iDCA, 0);
            (primMc > 1.e-1) ? fRatioDCAPrim.SetBinError(iDCA, TMath::Sqrt(primPrediction_err * primPrediction_err / primPrediction / primPrediction + primMc_err * primMc_err / primMc / primMc)) : fRatioDCAPrim.SetBinError(iDCA, 0.);

            double secPrediction = mc2->GetBinContent(iDCA);
            double secMc = fDCAMcProjSec->GetBinContent(iDCA);
            double secPrediction_err = mc2->GetBinError(iDCA);
            double secMc_err = fDCAMcProjSec->GetBinError(iDCA);
            (secMc > 1.e-1) ? fRatioDCASec.SetBinContent(iDCA, secPrediction / secMc) : fRatioDCASec.SetBinContent(iDCA, 0);
            (secMc > 1.e-1) ? fRatioDCASec.SetBinError(iDCA, TMath::Sqrt(secPrediction_err * secPrediction_err / secPrediction / secPrediction + secMc_err * secMc_err / secMc / secMc)) : fRatioDCASec.SetBinError(iDCA, 0.);
          }
          fRatioDCAPrim.GetXaxis()->SetTitle(kAxisTitleDCA);
          fRatioDCAPrim.GetYaxis()->SetTitle("Prediction / MC");
          fRatioDCASec.GetXaxis()->SetTitle(kAxisTitleDCA);
          fRatioDCASec.GetYaxis()->SetTitle("Prediction / MC");
          fRatioDCAPrim.GetYaxis()->SetRangeUser(0., 8.);
          fRatioDCASec.GetYaxis()->SetRangeUser(0., 8.);
          fRatioDCAPrim.Write();
          fRatioDCASec.Write();

          // write histograms to file
          fDCAdatProj->Write();
          fDCAMcProjPrim->Write();
          fDCAMcProjSec->Write();

          // save canvas plot
          canv.Print(Form("%s/primary_fraction/%s_%s/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutSettings, kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1], ptMin, ptMax));
        }
      } // end of loop on centrality bin

      // primary fraction fit with fFitFunc function
      gStyle->SetStatX(0.85);
      gStyle->SetStatY(0.5);
      gStyle->SetStatFontSize(0.035);
      TF1 fFitFunc(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsDeuteron[iCent][0], kCentBinsLimitsDeuteron[iCent][1]), "[0]+(1-[0])/(1+[1]*exp([2]*x))", 0.8f, 1.6f);
      fFitFunc.SetParName(0, "a");
      fFitFunc.SetParName(1, "b");
      fFitFunc.SetParName(2, "c");
      fFitFunc.SetParLimits(0, 0., 1.);
      fFitFunc.SetParLimits(1, 0., 10000.);
      fFitFunc.SetParLimits(2, -100., 0.);
      fPrimaryFrac.Fit(&fFitFunc, "MRQ");
      fPrimaryFrac.Fit(&fFitFunc, "MRQ");
      fFitFunc.Write();

      fPrimaryFrac.SetMarkerStyle(20);
      fPrimaryFrac.SetMarkerSize(0.8);
      fPrimaryFrac.SetOption("pe");
      fPrimaryFrac.GetYaxis()->SetTitle("#it{f}_{#it{prim}}");
      fPrimaryFrac.GetXaxis()->SetTitle(kAxisTitlePt);
      fPrimaryFrac.GetXaxis()->SetRangeUser(0.7, 5.0);
      fPrimaryFrac.Write();

      system(Form("mkdir %s/primary_plots", kPlotDir));
      TCanvas cPrim("cPrim", "cPrim");
      cPrim.cd();
      fPrimaryFrac.GetXaxis()->SetRangeUser(0.8, 1.6);
      fPrimaryFrac.GetYaxis()->SetRangeUser(0.0, 1.1);
      fPrimaryFrac.Draw("");
      cPrim.Print(Form("%s/primary_plots/%s.png", kPlotDir, fPrimaryFrac.GetName()));
    }
  }
  outFile->Close();
}