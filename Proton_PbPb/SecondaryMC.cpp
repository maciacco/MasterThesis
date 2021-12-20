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
using namespace proton;

bool use_uniform = false;
bool use_roofit = false;

void SecondaryMC(const char *cutSettings = "", const char *inFileDatName = "AnalysisResults_largeNsigma", const char *inFileMCName = "mc_20g7_20210929", const char *outFileName = "PrimaryProtonMC", const bool useAntiProtonsAsPrimaries = false)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetTextFont(44);
  // make signal extraction plots directory
  system(Form("mkdir %s/primary_fraction", kPlotDir));

  // open files
 //TFile *inFileDat = TFile::Open(Form("%s/%s.root", kDataDir, inFileDatName)); 
  TFile *inFileMC20g7_likeData = TFile::Open(Form("%s/%s.root", kDataDir, "mc_20g7_20210929")); 
  TFile *inFileMC21l5 = TFile::Open(Form("%s/%s.root", kDataDir, inFileMCName));
  TFile *inFileMC20g7 = TFile::Open(Form("%s/%s.root", kDataDir, "mc_20g7_likeData_largeNsigma"));
  //TFile *inFileMC1 = TFile::Open(Form("%s/%s.root", kDataDir, inFileMCName));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    double noSecMaterialThreshold = 1.79f;
    if (iMatt == 0)noSecMaterialThreshold = 0.f;

    // make plot subdirectory
    system(Form("mkdir %s/primary_fraction/%s_%s", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList(s)
    std::string listName_mcFalse = Form("mpuccio_proton_mcFalse_%s", cutSettings);
    std::string listName_mcTrue = Form("mpuccio_proton_mcTrue_%s", cutSettings);
    std::string listName = Form("nuclei_proton_%s", cutSettings);
    TTList *listMc21l5_likeData = (TTList *)inFileMC21l5->Get(listName_mcFalse.data());
    TTList *listMc21l5 = (TTList *)inFileMC21l5->Get(listName_mcTrue.data());
    TTList *listMc20g7_likeData = (TTList *)inFileMC20g7->Get(listName.data());
    TTList *listMc20g7 = (TTList *)inFileMC20g7_likeData->Get(listName.data());

    // get histograms from files
    TH3F *fDCAdat = (TH3F *)listMc21l5_likeData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAdat_20g7 = (TH3F *)listMc20g7_likeData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[iMatt]));
    if (ADD20g7)
      fDCAdat->Add(fDCAdat_20g7);
    TH3F *fDCAprim, *fDCAsec, *fDCAsecWD, *fDCAprim_20g7, *fDCAsec_20g7, *fDCAsecWD_20g7;
    if(useAntiProtonsAsPrimaries){
      fDCAprim = (TH3F *)listMc21l5_likeData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[0]));
    }
    else{
      fDCAprim = (TH3F *)listMc21l5->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
      fDCAsec = (TH3F *)listMc21l5->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
      fDCAsecWD = (TH3F *)listMc21l5->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
      fDCAprim_20g7 = (TH3F *)listMc20g7->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
      fDCAsec_20g7 = (TH3F *)listMc20g7->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
      fDCAsecWD_20g7 = (TH3F *)listMc20g7->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
      //fDCAprim2 = (TH3F *)listMc20e3a_2->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[0]));
      //fDCAprim3 = (TH3F *)listMc3->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[0]));
      //fDCAprim->Add(fDCAprim2);
      //fDCAprim->Add(fDCAprim3);
    }
    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      TH1D fPrimaryFrac(Form("f%sPrimFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);

      int nUsedPtBins = 24;

      for (int iPtBin = 5; iPtBin < nUsedPtBins + 1; ++iPtBin)
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
        TString canvTitleTOF;
        TString canvNameTOF;
        TH1D fRatioDCAPrim("", "", kNDCABinsMedium, kDCABinsMedium);
        TH1D fRatioDCASec("", "", kNDCABinsMedium, kDCABinsMedium);
        TH1D fRatioDCASecWD("", "", kNDCABinsMedium, kDCABinsMedium);

        TString projTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]));
        fDCAdatProj = fDCAdat->ProjectionZ(TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[iCent][0], kCentBinsProton[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAdatProj->SetTitle(projTitle);
        fDCAMcProjPrim = fDCAprim->ProjectionZ(TString::Format("f%sDCAPrimTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[iCent][0], kCentBinsProton[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjPrim->SetTitle(projTitle);
        fDCAMcProjSec = fDCAsec->ProjectionZ(TString::Format("f%sDCASecTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[iCent][0], kCentBinsProton[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSec->SetTitle(projTitle);
        fDCAMcProjSecWD = fDCAsecWD->ProjectionZ(TString::Format("f%sDCASecWDTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[iCent][0], kCentBinsProton[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSecWD->SetTitle(projTitle);

        canvTitleTOF = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]));
        canvNameTOF = TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax));
        TCanvas canv(canvNameTOF, canvTitleTOF);

        fDCAdatProj->GetYaxis()->SetTitle("Entries");
        fDCAMcProjPrim->GetYaxis()->SetTitle("Entries");

        fDCAdatProj->Write();
        fDCAMcProjPrim->Write();

        // compute wd template fraction
        // compute fraction of primaries and material secondaries
        double intPrimDCAcutError = 0.;
        double intPrimDCAcut = fDCAMcProjPrim->Integral(fDCAMcProjPrim->FindBin(-0.12), fDCAMcProjPrim->FindBin(0.115));
        double intSecDCAcut = fDCAMcProjSec->Integral(fDCAMcProjSec->FindBin(-0.12), fDCAMcProjSec->FindBin(0.115));
        double intSecWDDCAcut = fDCAMcProjSecWD->Integral(fDCAMcProjSecWD->FindBin(-0.12), fDCAMcProjSecWD->FindBin(0.115));
        double dataIntegralDCAcut = intPrimDCAcut+intSecDCAcut+intSecWDDCAcut;
        std::cout << "primary integral = " << intPrimDCAcut << std::endl;
        //double intSecDCAcut = mc3->Integral(result->FindBin(-0.12), result->FindBin(0.115));
        double primaryRatio = (intPrimDCAcut)/(dataIntegralDCAcut);
        double primaryRatioError = TMath::Sqrt(primaryRatio * (1.f - primaryRatio) / dataIntegralDCAcut);
        //double secondaryRatio = intSecDCAcut / intResDCAcut;
        //double secondaryRatioError = ;
        fPrimaryFrac.SetBinContent(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatio);
        fPrimaryFrac.SetBinError(fPrimaryFrac.FindBin(ptMin + 0.005f), primaryRatioError);
      }

      // primary fraction fit with fFitFunc function
      gStyle->SetStatX(0.85);
      gStyle->SetStatY(0.5);
      gStyle->SetStatFontSize(0.035);
      TF1 fFitFunc(Form("f%sFunctionFit_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), "pol2", 0.5f, 3.0f);
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