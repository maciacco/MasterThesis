#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLatex.h>

#include "../utils/Config.h"
#include "../utils/Utils.h"

using namespace he3;
using namespace utils;

void fitXsection(){
  gStyle->SetOptFit(00000);
  gStyle->SetOptStat(00000);

  TFile fileIn(Form("%s/He3inelXS_LHC18q_0_10_pass3_20211028.root",kDataDir));
  TFile fileOut("XSOut.root","recreate");

  TCanvas *canv = (TCanvas*)fileIn.Get("c1");
  TH1D *hXSectGEANT4 = (TH1D*)canv->GetPrimitive("He3_Cl_inelXS");
  TGraphErrors *gXSectMeasrd = (TGraphErrors*)canv->GetPrimitive("CrossSec_He3_pos");
  TCanvas c("c","c");

  // convert th1 to tf1
  TH2TF *th2tf=new TH2TF();
  th2tf->SetInputHist(hXSectGEANT4);
  auto inputG4He3=new TF1("inputG4He3",th2tf,&TH2TF::Eval,hXSectGEANT4->GetXaxis()->GetXmin(),hXSectGEANT4->GetXaxis()->GetXmax(),1);

  gXSectMeasrd->Fit("inputG4He3","WR","",0.5,7.);
  gXSectMeasrd->GetFunction("inputG4He3")->SetNpx(10000);
  gXSectMeasrd->GetYaxis()->SetRangeUser(-0.5,2.5);
  c.cd();
  gXSectMeasrd->Draw("ap");
  TLatex chi2(5,2.,Form("#chi^{2}/NDF = %.2f/%d",inputG4He3->GetChisquare(),inputG4He3->GetNDF()));
  TLatex par(5,1.5,Form("C = %.2f #pm %.2f",inputG4He3->GetParameter(0),inputG4He3->GetParError(0)));
  chi2.Draw("same");
  par.Draw("same");
  c.Write("fitToData");
}
