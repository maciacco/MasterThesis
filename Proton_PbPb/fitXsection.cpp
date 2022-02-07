#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>

#include "../utils/Config.h"
#include "../utils/Utils.h"

using namespace proton;
using namespace utils;

void fitXsection(){
  gStyle->SetOptFit(00000);
  gStyle->SetOptStat(00000);

  TFile fileInDataAntiP(Form("%s/HEPData-ins1797442-v1-Table_11.root",kDataDir));
  TFile fileInG4AntiP(Form("%s/CrossSectionsGeant4default.root",kDataDir));
  TFile fileOut("XSOut.root","recreate");

  TH1D *hXSectGEANT4 = (TH1D*)fileInG4AntiP.Get("hCSgeant4_ap_highp");
  TGraphAsymmErrors *gXSectMeasrd = (TGraphAsymmErrors*)fileInDataAntiP.Get("Table 11/Graph1D_y1");
  gXSectMeasrd->Write();
  TCanvas c("c","c");

  // convert th1 to tf1
  TH2TF *th2tf=new TH2TF();
  th2tf->SetInputHist(hXSectGEANT4);
  auto inputG4AntiP=new TF1("inputG4AntiP",th2tf,&TH2TF::Eval,0.95,2.7,1);

  gXSectMeasrd->Fit("inputG4AntiP","WR","",0.95,2.7);
  gXSectMeasrd->GetFunction("inputG4AntiP")->SetNpx(10000);
  gXSectMeasrd->GetYaxis()->SetRangeUser(-0.5,2.5);
  c.cd();
  gXSectMeasrd->Draw("ap");
  TLatex chi2(5,2.,Form("#chi^{2}/NDF = %.2f/%d",inputG4AntiP->GetChisquare(),inputG4AntiP->GetNDF()));
  TLatex par(5,1.5,Form("C = %.2f #pm %.2f",inputG4AntiP->GetParameter(0),inputG4AntiP->GetParError(0)));
  chi2.Draw("same");
  par.Draw("same");
  c.Write("fitToData");
}
