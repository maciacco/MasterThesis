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
using namespace proton;

bool use_uniform = false;

const double fitRange = 1.25;

void ComparePrimaryTemplates(const char *cutSettings = "", const double DCAxyCut=0.12, const char *inFileDatName = "AnalysisResults", const char *inFileMCName = "mc", const char *outFileName = "PrimaryProton", const bool use_roofit = false, const bool useAntiProtonsAsPrimaries = false)
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
  TFile *inFileDat = TFile::Open(Form("%s/%s_largeNsigma_largeBinningDCA.root", kDataDir, inFileDatName));
  TFile *inFileMC21l5 = TFile::Open(Form("%s/%s.root", kDataDir, "AnalysisResults_LHC21l5_full_largeDCA"/* inFileMCName */));
  TFile *inFileMC20e3a_1 = TFile::Open(Form("%s/%s_20e3a_runlist1_20210929.root", kDataDir, "mc"));
  TFile *inFileMC20e3a_2 = TFile::Open(Form("%s/%s_20e3a_runlist2_20210929.root", kDataDir, "mc"));
  //TFile *inFileMC1 = TFile::Open(Form("%s/%s.root", kDataDir, inFileMCName));
  TFile *outFile = TFile::Open(Form("%s/%s.root", kOutDir, outFileName), "recreate");

  for (int iMatt = 0; iMatt < 2; ++iMatt)
  {
    double noSecMaterialThreshold = 1.59f;/* 1.24f; */ //1.59f; // DO NOT USE SECONDARIES FROM MATERIAL
    if (iMatt == 0)
      noSecMaterialThreshold = 0.f;

    // make plot subdirectory
    system(Form("mkdir %s/primary_fraction/%s_%s", kPlotDir, kAntimatterMatter[iMatt], cutSettings));

    // get TTList(s)
    std::string listName_true = Form("nuclei_proton_mcTrue_%s", cutSettings);
    std::string listName = Form("nuclei_proton_%s", cutSettings);
    TTList *listData = (TTList *)inFileDat->Get(listName.data());
    TTList *listMc21l5 = (TTList *)inFileMC21l5->Get(listName_true.data());
    TTList *listMc20e3a_1 = (TTList *)inFileMC20e3a_1->Get(listName.data());
    TTList *listMc20e3a_2 = (TTList *)inFileMC20e3a_2->Get(listName.data());

    // get histograms from files
    TH3F *fDCAdat = (TH3F *)listData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAprim, *fDCAprim1, *fDCAprim2, *fDCAprim3;
    if (useAntiProtonsAsPrimaries)
    {
      fDCAprim = (TH3F *)listData->Get(Form("f%sDCAxyTOF", kAntimatterMatter[0]));
    }
    else
    {
      fDCAprim1 = (TH3F *)listMc20e3a_1->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
      fDCAprim2 = (TH3F *)listMc20e3a_2->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
      //fDCAprim3 = (TH3F *)listMc3->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[0]));
      fDCAprim = (TH3F *)fDCAprim1->Clone(fDCAprim1->GetName());
      fDCAprim->Add(fDCAprim2);
      //fDCAprim->Add(fDCAprim3);
    }
    TH3F *fDCAprim21l5 = (TH3F *)listMc21l5->Get(Form("f%sDCAPrimaryTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsec1 = (TH3F *)listMc20e3a_1->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsec2 = (TH3F *)listMc20e3a_2->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
    //TH3F *fDCAsec3 = (TH3F *)listMc3->Get(Form("f%sDCASecondaryTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsec = (TH3F *)fDCAsec1->Clone(fDCAsec1->GetName());
    fDCAsec->Add(fDCAsec2);
    //fDCAsec->Add(fDCAsec3);
    TH3F *fDCAsecWD1 = (TH3F *)listMc20e3a_1->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsecWD2 = (TH3F *)listMc20e3a_2->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
    //TH3F *fDCAsecWD3 = (TH3F *)listMc3->Get(Form("f%sDCASecondaryWeakTOF", kAntimatterMatter[iMatt]));
    TH3F *fDCAsecWD = (TH3F *)fDCAsecWD1->Clone(fDCAsecWD1->GetName());
    fDCAsecWD->Add(fDCAsecWD2);
    //fDCAsecWD->Add(fDCAsecWD3);

    for (int iCent = 0; iCent < kNCentClasses; ++iCent)
    {
      TH1D fKolmogorov(Form("f%sKolmogorovTest_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);
      TH1D fSecondaryFrac(Form("f%sSecFrac_%.0f_%.0f", kAntimatterMatter[iMatt], kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), Form("%.0f-%.0f%%", kCentBinsLimitsProton[iCent][0], kCentBinsLimitsProton[iCent][1]), kNPtBins, kPtBins);

      int nUsedPtBins = 24;

      for (int iPtBin = 5; iPtBin < nUsedPtBins + 1; ++iPtBin)
      { // loop on pT bins
        //fKolmogorov.SetBinContent(iPtBin, 1.);

        double ptMin = fKolmogorov.GetXaxis()->GetBinLowEdge(iPtBin);
        double ptMax = fKolmogorov.GetXaxis()->GetBinUpEdge(iPtBin);
        int pTbinsIndexMin = fDCAdat->GetYaxis()->FindBin(ptMin);
        int pTbinsIndexMax = fDCAdat->GetYaxis()->FindBin(ptMax - 0.005);
        outFile->cd();

        // project TH3 histogram
        TH1D *fDCAdatProj;
        TH1D *fDCAMcProjPrim;
        TH1D *fDCAMcProjPrim21l5;
        TH1D *fDCAMcProjSec;
        TH1D *fDCAMcProjSecWD;
        TString canvTitleTOF;
        TString canvNameTOF;
        TH1D fRatioDCAPrim("", "", kNDCABinsMediumOld, kDCABinsMediumOld);
        TH1D fRatioDCASec("", "", kNDCABinsMediumOld, kDCABinsMediumOld);

        TString projTitle = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]));
        fDCAdatProj = fDCAdat->ProjectionZ(TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[iCent][0], kCentBinsProton[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAdatProj->SetTitle(projTitle);
        fDCAMcProjPrim = fDCAprim->ProjectionZ(TString::Format("f%sDCAPrimaryTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[iCent][0], kCentBinsProton[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjPrim->SetTitle(projTitle);
        fDCAMcProjPrim21l5 = fDCAprim21l5->ProjectionZ(TString::Format("f%sDCAPrimaryTOF21l5_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[iCent][0], kCentBinsProton[iCent][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjPrim21l5->SetTitle(projTitle);

        fDCAdatProj->GetYaxis()->SetTitle("Entries");
        fDCAdatProj->Write();
        fDCAMcProjPrim->GetYaxis()->SetTitle("Entries");
        fDCAMcProjPrim->Write();
        fDCAMcProjPrim21l5->GetYaxis()->SetTitle("Entries");
        fDCAMcProjPrim21l5->Write();
        // kolmogorov test
        double kolm_test=fDCAMcProjPrim21l5->KolmogorovTest(fDCAMcProjPrim);
        std::cout<<"Cent = "<<iCent<<", Pt = "<<ptMin<<", Kolmogorov test = "<<kolm_test<<std::endl;
        fKolmogorov.SetBinContent(iPtBin,kolm_test);
        fDCAMcProjSec = fDCAsec->ProjectionZ(TString::Format("f%sDCASecondaryTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[3][0], kCentBinsProton[3][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSec->SetTitle(projTitle);
        fDCAMcProjSecWD = fDCAsecWD->ProjectionZ(TString::Format("f%sDCASecondaryWeakTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax)), kCentBinsProton[3][0], kCentBinsProton[3][1], pTbinsIndexMin, pTbinsIndexMax);
        fDCAMcProjSecWD->SetTitle(projTitle);

        // rebin
        /* 
        fDCAdatProj = (TH1D *)fDCAdatProj->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        fDCAMcProjPrim = (TH1D *)fDCAMcProjPrim->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        fDCAMcProjPrim21l5 = (TH1D *)fDCAMcProjPrim21l5->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        fDCAMcProjSec = (TH1D *)fDCAMcProjSec->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld);
        fDCAMcProjSecWD = (TH1D *)fDCAMcProjSecWD->Rebin(kNDCABinsMediumOld, fRatioDCASec.GetName(), kDCABinsMediumOld); */
        /* 
        fDCAdatProj = (TH1D*)fDCAdatProj->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
        fDCAMcProjPrim = (TH1D*)fDCAMcProjPrim->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
        fDCAMcProjSec = (TH1D*)fDCAMcProjSec->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
        fDCAMcProjSecWD = (TH1D*)fDCAMcProjSecWD->Rebin(kNDCABinsLarge, fRatioDCASec.GetName(), kDCABinsLarge);
 */
        canvTitleTOF = TString::Format("%.2f#leq #it{p}_{T}<%.2f GeV/#it{c}, %.0f-%.0f%%", fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax), fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]));
        canvNameTOF = TString::Format("f%sDCAxyTOF_%.0f_%.0f_%.2f_%.2f", kAntimatterMatter[iMatt], fDCAdat->GetXaxis()->GetBinLowEdge(kCentBinsProton[iCent][0]), fDCAdat->GetXaxis()->GetBinUpEdge(kCentBinsProton[iCent][1]), fDCAdat->GetYaxis()->GetBinLowEdge(pTbinsIndexMin), fDCAdat->GetYaxis()->GetBinUpEdge(pTbinsIndexMax));
        TCanvas canv(canvNameTOF, canvTitleTOF);

      }
      fKolmogorov.Write();
    }
  }
  outFile->Close();
}