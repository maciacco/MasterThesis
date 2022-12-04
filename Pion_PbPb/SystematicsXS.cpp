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

const char* pi_plus_minus[]={"#pi^{-}","#pi^{+}"};

void SystematicsXS(const char *cutSettings = "", const char *outFileName = "SystematicsXS", const char *outFileOption = "recreate", const char *effFile = "EfficiencyProtonMC_21l5_false_XS")
{
  std::cout << "cutSettings = " << cutSettings << std::endl;
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(44);

  TFile *inFileEff[3];
  inFileEff[0] = TFile::Open(Form("../data/Proton_PbPb/%s.root", effFile));
  inFileEff[1] = TFile::Open(Form("../data/Proton_PbPb/%s-.root", effFile));
  inFileEff[2] = TFile::Open(Form("../data/Proton_PbPb/%s+.root", effFile));
  if (!inFileEff)
  {
    std::cout << "Input files do not exist!" << std::endl;
    return;
  }

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
    eff_A[iXS] = (TH1D *)inFileEff[iXS]->Get(Form("%s_/f%sEff_TOF_-10_0", cutSettings, kAntimatterMatter[0]));
    eff_M[iXS] = (TH1D *)inFileEff[iXS]->Get(Form("%s_/f%sEff_TOF_-10_0", cutSettings, kAntimatterMatter[0]));

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
  double range_XS[][2]={{1.036-0.0246,1.036+0.0246},{1.1348-0.0132,1.1348+0.0132}};
  double x_XS[][3]={{0.8,0.913/0.7,1.5/0.7},{0.8,0.913/0.7,1.5/0.7}};
  for (int iMatt=0;iMatt<2;++iMatt){
    TCanvas cXS(Form("c%s",kAntimatterMatter[iMatt]),Form("c%s",kAntimatterMatter[iMatt]));
    double x_err_antip[]={0,0,0};
    double y_antip[]={0,fRatio[1][iMatt]->GetFunction("pol0")->GetParameter(0)/1.10632,fRatio[2][iMatt]->GetFunction("pol0")->GetParameter(0)/1.10632};
    double y_err_antip[]={0,fRatio[1][iMatt]->GetFunction("pol0")->GetParError(0)/1.10632,fRatio[2][iMatt]->GetFunction("pol0")->GetParError(0)/1.10632};
    TGraphErrors gXS(3,x_XS[iMatt],y_antip,x_err_antip,y_err_antip);
    //TF1 f("f","[0]*x*x+[1]*x+1-[0]-[1]");
    TF1 f("f","1/TMath::Exp([0])*TMath::Exp([0]*x)");
    gXS.Fit("f","QR","",0.86,2.2);
    //std::cout<<"eval = "<<f.Eval(0.7)<<std::endl;
    gXS.SetTitle(pi_plus_minus[iMatt]);
    gXS.SetMinimum(0);
    gXS.SetMaximum(1.1);
    gXS.SetMarkerStyle(20);
    gXS.SetMarkerSize(0.8);
    gXS.GetXaxis()->SetRangeUser(0.86,2.2);
    gXS.GetYaxis()->SetTitle("Efficiency ratio to default");
    gXS.GetXaxis()->SetTitle("Cross section scaling factor");
    gXS.SetMarkerSize(0.8);
    double lower_ratio = f.Eval(range_XS[iMatt][1]);
    double upper_ratio = f.Eval(range_XS[iMatt][0]);
    TLine lsx(range_XS[iMatt][0], 0, range_XS[iMatt][0], upper_ratio);
    lsx.SetLineStyle(kDashed);
    TLine ldx(range_XS[iMatt][1], 0, range_XS[iMatt][1], lower_ratio);
    ldx.SetLineStyle(kDashed);
    std::cout<<"Error = "<<(upper_ratio-lower_ratio)/2.<<std::endl;
    gXS.Draw("ape");
    lsx.Draw("same");
    ldx.Draw("same");
    cXS.Write();
  }
  outFile.Close();
}
