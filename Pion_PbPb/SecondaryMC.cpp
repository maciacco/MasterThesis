// SecondaryMC.cpp
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

#include "../utils/Utils.h"
#include "../utils/Config.h"

using namespace utils;
using namespace pion;

bool use_uniform = false;
bool use_roofit = false;

void SecondaryMC(const char *cutSettings = "", const double DCAxyCut = 0.12, const char *inFileDatName = "AnalysisResults_largeNsigma", const char *inFileMCName = "mc_20g7_20210929", const char *outFileName = "PrimaryPionMC")
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetTextFont(44);
  // make signal extraction plots directory
  system(Form("mkdir %s/primary_fraction", kPlotDir));

  // open files
  TFile *inFileMCInj = TFile::Open(Form("%s/%s.root", kDataDir, inFileMCName));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    double noSecMaterialThreshold = 1.99f;
    //if (iMatt == 0)noSecMaterialThreshold = 0.f;

    // make plot subdirectory
    system(Form("mkdir %s/primary_fraction/%s_%s", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList(s)
    std::string listName_mcFalse = Form("nuclei_pion_mcFalse_%s", cutSettings);
    std::string listName_mcTrue = Form("nuclei_pion_mcTrue_%s", cutSettings);
    std::string listName = Form("nuclei_pion_%s", cutSettings);
    TTList *listMcInj_likeData = (TTList *)inFileMCInj->Get(listName_mcFalse.data());
    TTList *listMcInj = (TTList *)inFileMCInj->Get(listName_mcTrue.data());

    // get histograms from files
    TH3F *fDCAdat = (TH3F *)listMcInj_likeData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAprim, *fDCAsec, *fDCAsecWD, *fDCAprim_20g7, *fDCAsec_20g7, *fDCAsecWD_20g7;
    
    fDCAprim = (TH3F *)listMcInj->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
    fDCAsec = (TH3F *)listMcInj->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
    fDCAsecWD = (TH3F *)listMcInj->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      TH1D fPrimaryFrac(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), kNPtBins, kPtBins);

      int nUsedPtBins = 27;

      for (int iPtBin = 7; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        fPrimaryFrac.SetBinContent(iPtBin, 1.);

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
        TH1D *fDCAdatProj_20g7;
        TH1D *fDCAMcProjPrim_20g7;
        TH1D *fDCAMcProjSec_20g7;
        TH1D *fDCAMcProjSecWD_20g7;
        TString canvTitleTOF;
        TString canvNameTOF;
        TH1D fRatioDCAPrim("", "", kNDCABinsMedium, kDCABinsMedium);
        TH1D fRatioDCASec("", "", kNDCABinsMedium, kDCABinsMedium);
        TH1D fRatioDCASecWD("", "", kNDCABinsMedium, kDCABinsMedium);

        TString projTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]));
        fDCAdatProj = fDCAdat->ProjectionZ(TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[iCent][0], kCentBinsPion[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAdatProj->SetTitle(projTitle);
        fDCAMcProjPrim = fDCAprim->ProjectionZ(TString::Format("f%sDCAPrimTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[iCent][0], kCentBinsPion[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjPrim->SetTitle(projTitle);
        fDCAMcProjSec = fDCAsec->ProjectionZ(TString::Format("f%sDCASecTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[iCent][0], kCentBinsPion[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSec->SetTitle(projTitle);
        fDCAMcProjSecWD = fDCAsecWD->ProjectionZ(TString::Format("f%sDCASecWDTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsPion[iCent][0], kCentBinsPion[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSecWD->SetTitle(projTitle);

        canvTitleTOF = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]));
        canvNameTOF = TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsPion[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsPion[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax));
        TCanvas canv(canvNameTOF, canvTitleTOF);

        fDCAdatProj->GetYaxis()->SetTitle("Entries");
        fDCAMcProjPrim->GetYaxis()->SetTitle("Entries");

        fDCAdatProj->Write();
        fDCAMcProjPrim->Write();

        // compute wd template fraction
        // compute fraction of primaries and material secondaries
        double intPrimDCAcutError = 0.;
        double intPrimDCAcut = fDCAMcProjPrim->Integral(fDCAMcProjPrim->FindBin(-DCAxyCut), fDCAMcProjPrim->FindBin(DCAxyCut-0.001));
        double intSecDCAcut = fDCAMcProjSec->Integral(fDCAMcProjSec->FindBin(-DCAxyCut), fDCAMcProjSec->FindBin(DCAxyCut-0.001));
        double intSecWDDCAcut = fDCAMcProjSecWD->Integral(fDCAMcProjSecWD->FindBin(-DCAxyCut), fDCAMcProjSecWD->FindBin(DCAxyCut-0.001));
        double dataIntegralDCAcut = intPrimDCAcut+intSecDCAcut+intSecWDDCAcut;
        if (kVerbose) std::cout << "primary integral = " << intPrimDCAcut << std::endl;
        double primaryRatio = (intPrimDCAcut)/(dataIntegralDCAcut);
        double primaryRatioError = TMath::Sqrt(primaryRatio * (1.f - primaryRatio) / dataIntegralDCAcut);
        fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatio);
        fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatioError);
      }

      // primary fraction fit with fFitFunc function
      gStyle->SetStatX(0.85);
      gStyle->SetStatY(0.5);
      gStyle->SetStatFontSize(0.035);
      TF1 fFitFunc(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsPion[iCent][0], kCentBinsLimitsPion[iCent][1]), "pol1", 1.0f, 2.0f);
      fPrimaryFrac.Fit(&fFitFunc, "MRQ");
      fPrimaryFrac.Fit(&fFitFunc, "MRQ");
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
      cPrim.Print(Form("%s/primary_plots/%s.pdf", kPlotDir, fPrimaryFrac.GetName()));
    
    } // end of loop on centrality bin
  }
  outFile->Close();
}