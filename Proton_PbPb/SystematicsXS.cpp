// SpystematicsXS.cpp
// This macro computes fully corrected spectra

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>

#include "../utils/Utils.h"
#include "../utils/Config.h"

using utils::TTList;
using namespace proton;

double protonCorrectionPt(int iMatt,double pt){
  if (iMatt == 1) 
    return 1;
    //return 0.99876*TMath::Power(pt,0.00036);
  //return 1.03176*TMath::Power(pt,-0.01249);
  return 1;
};

void SystematicsXS(const char *cutSettings = "", const double roi_nsigma = 8., const bool G3G4Prim = true, const bool binCounting = false, const int bkg_shape = 1, const bool sigmoidCorrection = true, const char *histoNameDir = ".", const char *outFileName = "SystematicsXS", const char *outFileOption = "recreate", const char *effFile = "EfficiencyProtonMC_21l5_false_XS")
{
  std::cout << "cutSettings = " << cutSettings << std::endl;
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(44);

  TFile *inFileEff[3];
  inFileEff[0] = TFile::Open(Form("%s/%s.root", kOutDir, effFile));
  inFileEff[1] = TFile::Open(Form("%s/%s-.root", kOutDir, effFile));
  inFileEff[2] = TFile::Open(Form("%s/%s+.root", kOutDir, effFile));
  if (!inFileEff)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

  int iNsigma = 0;
  if (roi_nsigma > 7.8 && roi_nsigma < 8.1) iNsigma = 1;
  else if (roi_nsigma > 8.1) iNsigma = 2;

  const char* XSLowUp[]={"","Low","Up"};

  TFile outFile(Form("%s/%s.root", kOutDir, outFileName), outFileOption);

  TH1D *fRatio[3][2];
  for (int iXS=0; iXS<3; ++iXS)
  {
    for (int iMatt=0; iMatt<2; ++iMatt){
      // std::cout << "read: " << Form("%s_%d_%d/fATOFrawYield_%.0f_%.0f", cutSettings, binCounting, bkg_shape, kCentBinsLimitsProton[iXS][0], kCentBinsLimitsProton[iXS][1]) << std::endl;
      fRatio[iXS][iMatt] = new TH1D(Form("f%sRatio%s",kAntimatterMatter[iMatt],XSLowUp[iXS]),Form("f%sRatio",kAntimatterMatter[iMatt]),kNPtBins,kPtBins);
    }
  }

  int pTbinMax = 24;
  TH1D *eff_A[3];
  TH1D *eff_M[3]; 
  for (int iXS=0;iXS<3;++iXS)
  {
    eff_A[iXS] = (TH1D *)inFileEff[iXS]->Get(Form("%s_/f%sEff_TOF_%.0f_%.0f", cutSettings, kAntimatterMatter[0], kCentBinsLimitsProton[4][0], kCentBinsLimitsProton[4][1]));
    eff_M[iXS] = (TH1D *)inFileEff[iXS]->Get(Form("%s_/f%sEff_TOF_%.0f_%.0f", cutSettings, kAntimatterMatter[1], kCentBinsLimitsProton[4][0], kCentBinsLimitsProton[4][1]));

    if (iXS > 0){
      for (int iPtBin = 5; iPtBin < pTbinMax + 1; ++iPtBin)
      {
        double NumXS = eff_A[iXS]->GetBinContent(iPtBin);
        double XS = eff_A[0]->GetBinContent(iPtBin);
        double NumXSErr = eff_A[iXS]->GetBinError(iPtBin);
        double XSErr = eff_A[0]->GetBinError(iPtBin);
        fRatio[iXS][0]->SetBinContent(iPtBin, NumXS / XS);
        fRatio[iXS][0]->SetBinError(iPtBin, NumXS / XS * TMath::Sqrt(NumXSErr * NumXSErr / NumXS / NumXS + XSErr * XSErr / XS / XS));

        NumXS = eff_M[iXS]->GetBinContent(iPtBin);
        XS = eff_M[0]->GetBinContent(iPtBin);
        NumXSErr = eff_M[iXS]->GetBinError(iPtBin);
        XSErr = eff_M[0]->GetBinError(iPtBin);
        fRatio[iXS][1]->SetBinContent(iPtBin, NumXS / XS);
        fRatio[iXS][1]->SetBinError(iPtBin, NumXS / XS * TMath::Sqrt(NumXSErr * NumXSErr / NumXS / NumXS + XSErr * XSErr / XS / XS));
      }
      for (int iMatt=0;iMatt<2;++iMatt){
        fRatio[iXS][iMatt]->GetXaxis()->SetTitle(kAxisTitlePt);
        fRatio[iXS][iMatt]->Fit("pol0");
        fRatio[iXS][iMatt]->GetXaxis()->SetRangeUser(1.0,2.0);
        fRatio[iXS][iMatt]->Write();
      }
    }
  }

  // antimatter curve
  TCanvas cAntip("cAntip","cAntip");
  double x_antip[]={0.913,1.5};
  double x_err_antip[]={0,0};
  double y_antip[]={fRatio[1][0]->GetFunction("pol0")->GetParameter(0),fRatio[2][0]->GetFunction("pol0")->GetParameter(0)};
  double y_err_antip[]={fRatio[1][0]->GetFunction("pol0")->GetParError(0),fRatio[2][0]->GetFunction("pol0")->GetParError(0)};
  TGraphErrors gAntip(2,x_antip,y_antip,x_err_antip,y_err_antip);
  //TF1 f("f","[0]*x*x+[1]*x+1-[0]-[1]");
  TF1 f("f","[0]*x+1-[0]");
  gAntip.Fit("f");
  double lower_ratio = f.Eval(1.075);
  double upper_ratio = f.Eval(0.925);
  std::cout<<"Error = "<<(upper_ratio-lower_ratio)/2.<<std::endl;
  gAntip.Draw("ape");
  cAntip.Write();
  outFile.Close();
}
