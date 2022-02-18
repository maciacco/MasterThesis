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
  TGraphAsymmErrors gXsectMeasrdShift(*gXSectMeasrd);
  for (int ip=0;ip<gXSectMeasrd->GetN();++ip){
    double tmpX=gXsectMeasrdShift.GetPointX(ip);
    double tmpY=gXsectMeasrdShift.GetPointY(ip);
    double tmpYerr=gXsectMeasrdShift.GetErrorYhigh(ip);
    gXsectMeasrdShift.SetPointY(ip,tmpY+tmpYerr);
  }
  TF1 inputG4AntiPCopy(*inputG4AntiP);
  gXsectMeasrdShift.Fit(&inputG4AntiPCopy,"MR","",1.0,2.7);
  gXSectMeasrd->Fit("inputG4AntiP","MR","",1.0,2.7);
  gXSectMeasrd->GetFunction("inputG4AntiP")->SetNpx(10000);
  gXSectMeasrd->GetYaxis()->SetRangeUser(-0.5,2.5);
  c.cd();
  gXSectMeasrd->Draw("ap");
  inputG4AntiPCopy.Draw("same");
  TLatex chi2(1,2.,Form("#chi^{2}/NDF = %.2f/%d",inputG4AntiP->GetChisquare(),inputG4AntiP->GetNDF()));
  TLatex par(1,1.5,Form("C = %.2f #pm %.2f",inputG4AntiP->GetParameter(0),inputG4AntiP->GetParError(0)));
  chi2.Draw("same");
  par.Draw("same");
  c.Write("fitToData");
  std::cout<<"Error: "<<inputG4AntiPCopy.GetParameter(0)-inputG4AntiP->GetParameter(0)<<std::endl;
}
