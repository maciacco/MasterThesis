#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLatex.h>

#include "../utils/Config.h"
#include "../utils/Utils.h"

using namespace triton;
using namespace utils;

void fitXsection(){
  gStyle->SetOptFit(00000);
  gStyle->SetOptStat(00000);

  TFile fileIn(Form("%s/FinalInelasticCrossSections_pp_PbPb_antitriton_antiHe3.root",kDataDir));
  TFile fileOut("XSOut.root","recreate");

  TH1D *hXSectGEANT4 = (TH1D*)fileIn.Get("hGEANT4param");
  TGraphErrors *gXSectMeasrdArray[2];
  gXSectMeasrdArray[0] = (TGraphErrors*)fileIn.Get("hInelXS_H3PbPb_SummedUnc");
  gXSectMeasrdArray[1] = (TGraphErrors*)fileIn.Get("hInelXS_H3pp_SummedUnc");
  TGraphErrors *gXSectMeasrd = new TGraphErrors();
  for (int i=0; i<6; ++i){
    gXSectMeasrd->AddPoint(gXSectMeasrdArray[i/3]->GetPointX(i%3),gXSectMeasrdArray[i/3]->GetPointY(i%3));
    gXSectMeasrd->SetPointError(i,0.,gXSectMeasrdArray[i/3]->GetErrorY(i%3));
    std::cout<<"index = "<<i%3<<"; X = "<<gXSectMeasrdArray[i/3]->GetPointX(i%3)<<"; Y = "<<gXSectMeasrdArray[i/3]->GetPointY(i%3)<<"; err = "<<gXSectMeasrdArray[i/3]->GetErrorY(i%3)<<std::endl;
  }
  TCanvas c("c","c");

  // convert th1 to tf1
  TH2TF *th2tf=new TH2TF();
  th2tf->SetInputHist(hXSectGEANT4);
  auto inputG4He3=new TF1("inputG4He3",th2tf,&TH2TF::Eval,hXSectGEANT4->GetXaxis()->GetXmin(),hXSectGEANT4->GetXaxis()->GetXmax(),1);

  gXSectMeasrd->Fit("inputG4He3","WR","",0.5,3.);
  gXSectMeasrd->GetFunction("inputG4He3")->SetNpx(10000);
  gXSectMeasrd->GetYaxis()->SetRangeUser(-0.5,2.5);
  gXSectMeasrd->SetMarkerStyle(20);
  gXSectMeasrd->SetMarkerSize(1);
  c.cd();
  gXSectMeasrd->Draw("ap");
  TLatex chi2(2,2.,Form("#chi^{2}/NDF = %.2f/%d",inputG4He3->GetChisquare(),inputG4He3->GetNDF()));
  TLatex par(2,1.5,Form("C = %.2f #pm %.2f",inputG4He3->GetParameter(0),inputG4He3->GetParError(0)));
  chi2.Draw("same");
  par.Draw("same");
  c.Write("fitToData");
}
