#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TLatex.h>

Double_t scalingFactorMine(Double_t *x, Double_t *par){
  Double_t f=1.-(1./0.029)* /* (0.738506/1.058)*(0.738506-1.) */ (1.058-1.)*par[0]*TMath::Power(x[0],-0.19483);
  Double_t f_new=1.-(1./0.029)* (0.738506/1.058)*(0.738506-1.) *par[0]*TMath::Power(x[0],-0.19483);
  Double_t f_note=0.99274*TMath::Power(x[0],0.00143);
  return f_new*(1+(f-f_note)/f_note);
}

Double_t scalingFactorNote(Double_t *x, Double_t *par){
  Double_t note = par[0]*TMath::Power(x[0],-0.01525);
  Double_t f_mine_note = 1.+(0.02088*TMath::Power(x[0],-0.48766))/(0.084)*(1.-0.83);
  return f_mine_note;
}

Double_t scalingFactorNoteM(Double_t *x, Double_t *par){
  return par[0]*TMath::Power(x[0],0.00143);
}

Double_t fpol5(Double_t *x, Double_t *par)
{
  if (x[0]< 2.4)
    return par[32]*(par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4) + par[5]*pow(x[0],5) + par[6]*pow(x[0],6) + par[7]*pow(x[0],7));
  else if (x[0]>2.4 && x[0]<3.15)
    return par[32]*(par[8] + par[9]*x[0] + par[10]*pow(x[0],2) + par[11]*pow(x[0],3) + par[12]*pow(x[0],4) + par[13]*pow(x[0],5) + par[14]*pow(x[0],6) + par[15]*pow(x[0],7));
  else if (x[0]>3.15 && x[0]<4.19)
    return par[32]*(par[16] + par[17]*x[0] + par[18]*pow(x[0],2) + par[19]*pow(x[0],3) + par[20]*pow(x[0],4) + par[21]*pow(x[0],5) + par[22]*pow(x[0],6) + par[23]*pow(x[0],7));
  return par[32]*(par[24] + par[25]*x[0] + par[26]*pow(x[0],2) + par[27]*pow(x[0],3) + par[28]*pow(x[0],4) + par[29]*pow(x[0],5) + par[30]*pow(x[0],6) + par[31]*pow(x[0],7));
}

void fitXsection(){
  gStyle->SetOptFit(00000);
  gStyle->SetOptStat(00000);

  TFile fileIn("He3inelXS_LHC18q_0_10_pass3_20211028.root");
  TFile fileOut("XSOut.root","recreate");

  TCanvas *canv = (TCanvas*)fileIn.Get("c1");
  TH1D *hXSectGEANT4 = (TH1D*)canv->GetPrimitive("He3_Cl_inelXS");
  TGraphErrors *gXSectMeasrd = (TGraphErrors*)canv->GetPrimitive("CrossSec_He3_pos");
  TCanvas c("c","c");

  TF1 f("fitPol","pol7",0.5,10./* ,13 */);
  hXSectGEANT4->Fit("fitPol","WR","",0.8,2.39);
  TF1 f1("fitPol1","pol4",0.5,10./* ,13 */);
  hXSectGEANT4->Fit("fitPol1","WR","",2.41,3.14);
  TF1 f12("fitPol12","pol5",0.5,10./* ,13 */);
  hXSectGEANT4->Fit("fitPol12","WR","",3.16,4.18);
  TF1 f123("fitPol123","pol7",0.5,10./* ,13 */);
  hXSectGEANT4->Fit("fitPol123","WR","",4.20,10.);
  TF1 f2("fitPol2","fpol5",0.5,10.,33);
  f2.FixParameter(0,f.GetParameter(0));
  f2.FixParameter(1,f.GetParameter(1));
  f2.FixParameter(2,f.GetParameter(2));
  f2.FixParameter(3,f.GetParameter(3));
  f2.FixParameter(4,f.GetParameter(4));
  f2.FixParameter(5,f.GetParameter(5));
  f2.FixParameter(6,f.GetParameter(6));
  f2.FixParameter(7,f.GetParameter(7));
  f2.FixParameter(8,f1.GetParameter(0));
  f2.FixParameter(9,f1.GetParameter(1));
  f2.FixParameter(10,f1.GetParameter(2));
  f2.FixParameter(11,f1.GetParameter(3));
  f2.FixParameter(12,f1.GetParameter(4));
  f2.FixParameter(13,0/* f1.GetParameter(5) */);
  f2.FixParameter(14,0/* f1.GetParameter(6) */);
  f2.FixParameter(15,0/* f1.GetParameter(7) */);
  f2.FixParameter(16,f12.GetParameter(0));
  f2.FixParameter(17,f12.GetParameter(1));
  f2.FixParameter(18,f12.GetParameter(2));
  f2.FixParameter(19,f12.GetParameter(3));
  f2.FixParameter(20,f12.GetParameter(4));
  f2.FixParameter(21,f12.GetParameter(5));
  f2.FixParameter(22,0/* f12.GetParameter(6) */);
  f2.FixParameter(23,0/* f12.GetParameter(7) */);
  f2.FixParameter(24,f123.GetParameter(0));
  f2.FixParameter(25,f123.GetParameter(1));
  f2.FixParameter(26,f123.GetParameter(2));
  f2.FixParameter(27,f123.GetParameter(3));
  f2.FixParameter(28,f123.GetParameter(4));
  f2.FixParameter(29,f123.GetParameter(5));
  f2.FixParameter(30,f123.GetParameter(6));
  f2.FixParameter(31,f123.GetParameter(7));
  gXSectMeasrd->Fit("fitPol2","WR","",0.5,10.);
  gXSectMeasrd->GetFunction("fitPol2")->SetNpx(10000);
  gXSectMeasrd->GetYaxis()->SetRangeUser(-0.5,2.5);
  c.cd();
  gXSectMeasrd->Draw("ap");
  TLatex chi2(5,2.,Form("#chi^{2}/NDF = %.2f/%d",f2.GetChisquare(),f2.GetNDF()));
  TLatex par(5,1.5,Form("C = %.2f #pm %.2f",f2.GetParameter(32),f2.GetParError(32)));
  chi2.Draw("same");
  par.Draw("same");
  c.Write("fitToData");

  TCanvas c1("c1","c1");
  hXSectGEANT4->Draw();
  f.Draw("same");
  f1.Draw("same");
  f12.Draw("same");
  f123.Draw("same");
  c1.Write("fitToGeant");
  hXSectGEANT4->Write();

  TF1 fMine("scaleMine",scalingFactorMine,1.,10.,1);
  fMine.SetParameter(0,0.00294);
  TF1 fNote("scaleNote",scalingFactorNote,1.,10.,1);
  fNote.SetParameter(0,1.04948);
  //fNote.SetParameter(0,0.99274);
  fMine.Write();
  fNote.Write();
  TCanvas cc("cc","cc");
  cc.cd();
  fMine.Draw("");
  fNote.Draw("same");
  cc.Write("TryScalingFactor");
}
