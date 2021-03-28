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

#include "Utils.h"
#include "Config.h"

void Secondary(const float cutDCAz = 1.f, const int cutTPCcls = 89, const char *inFileDatName = "TreeOutData", const char *inFileMCName = "TreeOutMC", const char *inFileWDName = "EfficiencyHe3SecWd", const char *outFileName = "PrimaryHe3", const bool useWdInFit = true, const bool rebinLowStatisticsBin = false)
{
  // make signal extraction plots directory
  system(Form("mkdir %s/primary_fraction", kPlotDir));

  // open files
  TFile *inFileDat = TFile::Open(Form("%s/%s.root", kResDir, inFileDatName));
  TFile *inFileMC = TFile::Open(Form("%s/%s.root", kResDir, inFileMCName));
  TFile *inFileWD = TFile::Open(Form("%s/%s.root", kOutDir, inFileWDName));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    // make plot subdirectory
    system(Form("mkdir %s/primary_fraction/%s_%1.1f_%d", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls));

    // get histograms from files
    TH3F *fDCAdat = (TH3F *)inFileDat->Get(Form("%1.1f_%d_/f%sDCAxyTPC", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));
    TH3F *fDCAprim = (TH3F *)inFileMC->Get(Form("%1.1f_%d_/f%sDCAPrimary", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));
    TH3F *fDCAsec = (TH3F *)inFileMC->Get(Form("%1.1f_%d_/fMDCASecondary", cutDCAz, cutTPCcls));
    TH3F *fDCAsecWeak = (TH3F *)inFileMC->Get(Form("%1.1f_%d_/f%sDCASecondaryWeak", cutDCAz, cutTPCcls, kAntimatterMatter[iMatt]));

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      // get wd fraction
      TH1D *fWDfrac = (TH1D *)inFileWD->Get(Form("f%sEff_TPC_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]));
      if (!fWDfrac)
      {
        std::cout << "fWDfrac histogram not found!" << std::endl;
        return;
      }

      // define pt bins
      const int nPtBinsSec = 18;
      double pTbins[nPtBinsSec + 1];
      pTbins[0] = kLowestPt;
      for (int iPtBins = 1; iPtBins < nPtBinsSec + 1; ++iPtBins)
      {
        pTbins[iPtBins] = pTbins[iPtBins - 1] + kDeltaPt;
      }
      TH1D fPrimaryFrac(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), nPtBinsSec, pTbins);
      TH1D fSecondaryFrac(Form("f%sSecFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), nPtBinsSec, pTbins);

      int nUsedPtBins = 12;
      if (iCent == kNCentClasses - 1) // rebin
      {
        int nPtBins = 13;
        double pTbinsNew[] = {1.f, 1.5f, 2.f, 2.5f, 3.f, 3.5f, 4.f, 4.5f, 5.f, 5.5f, 6.f, 7.f, 8.f, 10.f};
        fPrimaryFrac.SetBins(nPtBins, pTbinsNew);
        fSecondaryFrac.SetBins(nPtBins, pTbinsNew);
        nUsedPtBins = 11;
      }

      for (int iPtBin = 3; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        int lowerPtBinIndex = fDCAdat->GetYaxis()->FindBin(fPrimaryFrac.GetBinLowEdge(iPtBin));
        int upperPtBinIndex = fDCAdat->GetYaxis()->FindBin(fPrimaryFrac.GetBinLowEdge(iPtBin + 1) - 0.005f);
        double minPt = fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex);
        double maxPt = fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex);
        outFile->cd();

        if ((iMatt == 1) && (minPt < 2.95f))
        {
          // project TH3 histogram
          TH1D *fDCAdatProj;
          TH1D *fDCAMcProjPrim;
          TH1D *fDCAMcProjSec;
          TH1D *fDCAMcProjSecWeak;
          TString canvTitleTPC;
          TString canvNameTPC;
          TH1D fRatioDCAPrim("", "", kNDCABins, kDCABins);
          TH1D fRatioDCASec("", "", kNDCABins, kDCABins);

          TString projTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsHe3[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsHe3[iCent][1]));
          fDCAdatProj = fDCAdat->ProjectionZ(TString::Format("f%sDCAxyTPC_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsHe3[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsHe3[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex)), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1], lowerPtBinIndex, upperPtBinIndex);
          fDCAdatProj->SetTitle(projTitle);
          fDCAMcProjPrim = fDCAprim->ProjectionZ(TString::Format("f%sDCAPrimary_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsHe3[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsHe3[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex)), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1], lowerPtBinIndex, upperPtBinIndex);
          fDCAMcProjPrim->SetTitle(projTitle);
          fDCAMcProjSec = fDCAsec->ProjectionZ(TString::Format("f%sDCASecondary_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsHe3[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsHe3[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex)), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1], lowerPtBinIndex, upperPtBinIndex);
          fDCAMcProjSec->SetTitle(projTitle);
          fDCAMcProjSecWeak = fDCAsecWeak->ProjectionZ(TString::Format("f%sDCASecondaryWeak_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsHe3[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsHe3[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex)), kCentBinsHe3[iCent][0], kCentBinsHe3[iCent][1], lowerPtBinIndex, upperPtBinIndex);
          fDCAMcProjSecWeak->SetTitle(projTitle);
          canvTitleTPC = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsHe3[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsHe3[iCent][1]));
          canvNameTPC = TString::Format("f%sDCAxyTPC_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsHe3[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsHe3[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex));
          TCanvas canv(canvNameTPC, canvTitleTPC);

          fDCAdatProj->GetYaxis()->SetTitle("Entries");
          fDCAMcProjPrim->GetYaxis()->SetTitle("Entries");
          fDCAMcProjSec->GetYaxis()->SetTitle("Entries");
          fDCAMcProjSecWeak->GetYaxis()->SetTitle("Entries");

          // DCA rebin for low statistics bin
          if ((iCent < kNCentClasses - 1) && rebinLowStatisticsBin)
          {
            fDCAdatProj = (TH1D *)fDCAdatProj->Rebin(kNDCABinsLarge, fDCAdatProj->GetName(), kDCABinsLarge);
            fDCAMcProjPrim = (TH1D *)fDCAMcProjPrim->Rebin(kNDCABinsLarge, fDCAMcProjPrim->GetName(), kDCABinsLarge);
            fDCAMcProjSec = (TH1D *)fDCAMcProjSec->Rebin(kNDCABinsLarge, fDCAMcProjSec->GetName(), kDCABinsLarge);
            fDCAMcProjSecWeak = (TH1D *)fDCAMcProjSecWeak->Rebin(kNDCABinsLarge, fDCAMcProjSecWeak->GetName(), kDCABinsLarge);
            fRatioDCAPrim = *(TH1D *)fRatioDCAPrim.Rebin(kNDCABinsLarge, fRatioDCAPrim.GetName(), kDCABinsLarge);
            fRatioDCASec = *(TH1D *)fRatioDCASec.Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
          }

          // fit data distribution with TFractionFitter
          TObjArray *mc = new TObjArray(2); // MC histograms are put in this array
          mc->Add(fDCAMcProjPrim);
          mc->Add(fDCAMcProjSec);
          if (useWdInFit)
            mc->Add(fDCAMcProjSecWeak);
          TFractionFitter *fit = new TFractionFitter(fDCAdatProj, mc, "Q"); // initialise
          ROOT::Fit::Fitter *fitter = fit->GetFitter();

          // compute wd template fraction
          double dataIntegralDCAcut = fDCAdatProj->Integral(fDCAdatProj->FindBin(-0.1), fDCAdatProj->FindBin(0.09));
          double wdIntegralDCAcut = fDCAMcProjSecWeak->Integral(fDCAMcProjSecWeak->FindBin(-0.1), fDCAMcProjSecWeak->FindBin(0.09));
          double dataIntegral = fDCAdatProj->Integral();
          double wdIntegral = fDCAMcProjSecWeak->Integral();
          double wdTemplateFrac = (fWDfrac->GetBinContent(iPtBin)) * (dataIntegralDCAcut / dataIntegral) / (wdIntegralDCAcut / wdIntegral);

          if (useWdInFit)
          {
            fitter->Config().ParSettings(2).Set("wd", wdTemplateFrac);
            fitter->Config().ParSettings(2).Fix();
          }
          fit->SetRangeX(fDCAdatProj->FindBin(-1.3), fDCAdatProj->FindBin(1.2));
          fit->Constrain(0, 0.0, 1.0); // constrain fraction 1 to be between 0 and 1
          fit->Constrain(1, 0.0, 0.6);

          Int_t status = fit->Fit(); // perform the fit

          if (status == 0)
          { // check on fit status
            TH1F *result = (TH1F *)fit->GetPlot();
            TH1F *mc1 = (TH1F *)fit->GetMCPrediction(0);
            TH1F *mc2 = (TH1F *)fit->GetMCPrediction(1);
            mc1->SetName(Form("%s_Prediction", fDCAMcProjPrim->GetName()));
            mc2->SetName(Form("%s_Prediction", fDCAMcProjSec->GetName()));
            double integralData = fDCAdatProj->Integral();
            double fracMc1, fracMc2, fracMc3, errFracMc1, errFracMc2, errFracMc3;
            fit->GetResult(0, fracMc1, errFracMc1);
            fit->GetResult(1, fracMc2, errFracMc2);
            fit->GetResult(2, fracMc3, errFracMc3);
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

            TLegend leg(0.574499, 0.60552, 0.879628, 0.866667);
            leg.AddEntry(result, "Fit");
            leg.AddEntry(mc1, "Primary ^{3}He");
            leg.AddEntry(mc2, "Secondary ^{3}He");
            // draw on canvas
            double chi2 = fit->GetChisquare();
            gStyle->SetOptStat(0);
            fDCAdatProj->Scale(1, "width");
            result->Scale(1, "width");
            fDCAdatProj->Draw("Ep");
            result->Draw("histosame");
            fDCAdatProj->Draw("Epsame");

            if (useWdInFit)
            {
              TH1F *mc3 = (TH1F *)fit->GetMCPrediction(2);
              mc3->SetName(Form("%s_Prediction", fDCAMcProjSecWeak->GetName()));
              mc3->Scale(fracMc3 * integralData / mc3->Integral(), "width");
              mc3->SetLineColor(kRed);
              mc3->Write();

              leg.AddEntry(mc3, "Secondary ^{3}He, weak decays");
              mc3->Draw("histosame");
            }
            mc2->Draw("histosame");
            mc1->Draw("histosame");

            leg.Draw("same");
            TLatex chiSq(-1., 0.2 * result->GetMaximum(), Form("#chi^{2}/NDF=%.2f", chi2/fit->GetNDF()));
            TLatex prob(-1., 0.1 * result->GetMaximum(), Form("prob=%.7f", fit->GetProb()));
            // TLatex frac1(-1, 0.75 * result->GetMaximum(), Form("%.3f +- %.3f", fracMc1, errFracMc1));
            // TLatex frac2(-1, 0.3 * result->GetMaximum(), Form("%.3f +- %.3f", fracMc2, errFracMc2));
            // TLatex frac3(-1., 0.5 * result->GetMaximum(), Form("%.3f +- %.3f", fracMc3, errFracMc3));
            chiSq.Draw("same");
            prob.Draw("same");
            // frac1.Draw("same");
            // frac2.Draw("same");
            // if (useWdInFit)
            //   frac3.Draw("same");

            // compute fraction of primaries and material secondaries
            double intPrimDCAcut = mc1->Integral(result->FindBin(-0.1), result->FindBin(0.9));
            double intSecDCAcut = mc2->Integral(result->FindBin(-0.1), result->FindBin(0.9));
            double intResDCAcut = result->Integral(result->FindBin(-0.1), result->FindBin(0.9));
            double primaryRatio = intPrimDCAcut / intResDCAcut;
            double primaryRatioError = TMath::Sqrt(primaryRatio * (1.f - primaryRatio) / intResDCAcut);
            double secondaryRatio = intSecDCAcut / intResDCAcut;
            double secondaryRatioError = TMath::Sqrt(secondaryRatio * (1.f - secondaryRatio) / intResDCAcut);
            fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(minPt + 0.005f), primaryRatio);
            fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(minPt + 0.005f), primaryRatioError);
            fSecondaryFrac.SetBinContent(fSecondaryFrac.FindBin(minPt + 0.005f), secondaryRatio);
            fSecondaryFrac.SetBinError(fSecondaryFrac.FindBin(minPt + 0.005f), secondaryRatioError);

            canv.SetLogy();
            canv.Write();

            // write fitted histograms
            mc1->Write();
            mc2->Write();

            // scale starting histograms
            fDCAMcProjPrim->Scale(fracMc1 * integralData / fDCAMcProjPrim->Integral(), "width");
            fDCAMcProjSec->Scale(fracMc2 * integralData / fDCAMcProjSec->Integral(), "width");
            if (useWdInFit)
              fDCAMcProjSecWeak->Scale(fracMc3 * integralData / fDCAMcProjSecWeak->Integral(), "width");

            // compare template with predictions from fit
            fRatioDCAPrim.SetName(Form("f%sRatioDCAPrim_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1], fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex)));
            fRatioDCAPrim.SetTitle(fRatioDCAPrim.GetName());
            fRatioDCASec.SetName(Form("f%sRatioDCASec_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1], fDCAdat->GetYaxis()->GetBinLowEdge(lowerPtBinIndex), fDCAdat->GetYaxis()->GetBinUpEdge(upperPtBinIndex)));
            fRatioDCASec.SetTitle(fRatioDCASec.GetName());
            for (int iDCA = 1; iDCA < kNDCABins + 1; ++iDCA)
            {
              double primPrediction = mc1->GetBinContent(iDCA);
              double primMc = fDCAMcProjPrim->GetBinContent(iDCA);
              (primMc > 1.e-1) ? fRatioDCAPrim.SetBinContent(iDCA, primPrediction / primMc) : fRatioDCAPrim.SetBinContent(iDCA, 0);
              (primMc > 1.e-1) ? fRatioDCAPrim.SetBinError(iDCA, 0.) : fRatioDCAPrim.SetBinError(iDCA, 0.);

              double secPrediction = mc2->GetBinContent(iDCA);
              double secMc = fDCAMcProjSec->GetBinContent(iDCA);
              (secMc > 1.e-1) ? fRatioDCASec.SetBinContent(iDCA, secPrediction / secMc) : fRatioDCASec.SetBinContent(iDCA, 0);
              (secMc > 1.e-1) ? fRatioDCASec.SetBinError(iDCA, 0.) : fRatioDCASec.SetBinError(iDCA, 0.);
            }
            fRatioDCAPrim.GetXaxis()->SetTitle(kAxisTitleDCA);
            fRatioDCAPrim.GetYaxis()->SetTitle("Prediction / MC");
            fRatioDCASec.GetXaxis()->SetTitle(kAxisTitleDCA);
            fRatioDCASec.GetYaxis()->SetTitle("Prediction / MC");
            fRatioDCAPrim.GetYaxis()->SetRangeUser(0., 5.);
            fRatioDCASec.GetYaxis()->SetRangeUser(0., 5.);
            fRatioDCAPrim.Write();
            fRatioDCASec.Write();

            // write histograms to file
            fDCAdatProj->Write();
            fDCAMcProjPrim->Write();
            fDCAMcProjSec->Write();
            if (useWdInFit)
              fDCAMcProjSecWeak->Write();

            // save canvas plot
            canv.Print(Form("%s/primary_fraction/%s_%1.1f_%d/cent_%.0f_%.0f_pt_%.2f_%.2f.png", kPlotDir, kAntimatterMatter[iMatt], cutDCAz, cutTPCcls, kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1], minPt, maxPt));
          }
        }
        else
        {
          fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(minPt + 0.005f), 1.f - fWDfrac->GetBinContent(fWDfrac->FindBin(minPt + 0.005f)));
          fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(minPt + 0.005f), fWDfrac->GetBinError(fWDfrac->FindBin(minPt + 0.005f)));
        }
      } // end of loop on centrality bin

      // primary fraction fit with fSigmoid function
      TF1 fSigmoid(Form("f%sSigmoidFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsHe3[iCent][0], kCentBinsLimitsHe3[iCent][1]), "[0]/(1+exp([1]*x+[2]))", 2.f, 7.f);
      fSigmoid.SetParameter(0, 0.9f);
      if (iMatt == 1)
      {
        fSigmoid.SetParLimits(1, -10., 0.);
        fSigmoid.SetParLimits(2, 0.f, 100.f);
      }
      else
      {
        fSigmoid.FixParameter(1, 0.f);
        fSigmoid.FixParameter(2, 0.f);
      }
      fPrimaryFrac.Fit(&fSigmoid, "RQ");
      fSigmoid.Write();

      fPrimaryFrac.SetMarkerStyle(20);
      fPrimaryFrac.SetMarkerSize(0.8);
      fPrimaryFrac.SetOption("pe");
      fPrimaryFrac.GetYaxis()->SetTitle("#it{f}_{#it{prim}}");
      fPrimaryFrac.GetXaxis()->SetTitle(kAxisTitlePt);
      fPrimaryFrac.Write();
      if (iMatt == 1)
      {
        fSecondaryFrac.SetMarkerStyle(20);
        fSecondaryFrac.SetMarkerSize(0.8);
        fSecondaryFrac.SetOption("pe");
        fSecondaryFrac.GetYaxis()->SetTitle("#it{f}_{#it{sec}}");
        fSecondaryFrac.GetXaxis()->SetTitle(kAxisTitlePt);
        fSecondaryFrac.Write();
      }
    }
  }
  outFile->Close();
}