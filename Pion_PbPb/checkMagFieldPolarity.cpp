#include "../utils/Config.h"

constexpr const char* f_name = "out/SpectraPion_%s.root";
constexpr const char* f_name_out = "out/checkMagFieldPolarity.root";
constexpr const char* period[] = {"q","r"};
constexpr const char* opt[] = {"pe","samepe"};
constexpr Color_t colors[] = {kBlue,kRed};

using namespace pion;

void checkMagFieldPolarity(){
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TFile *f[2];
  TFile f_out(f_name_out,"recreate");
  for (int i=0; i<2; ++i) f[i]=TFile::Open(Form(f_name,period[i]));
  TH1D h("hSys","hSys",5,0,5);
  for (int iC_{0}; iC_ < kNCentClasses; ++iC_){
    int iC=iC_;
    if (iC==2||iC==3) iC=3-(iC_/3);
    std::cout<<iC<<std::endl;
    TCanvas c(Form("c_%d",iC),Form("c_%d",iC));
    TLatex t;
    t.SetTextFont(44);
    t.SetTextSize(25);
    TH1D* r[2];
    double ratio[2];
    double ratio_err[2];
    TLegend l(0.609023,0.205217,0.87594,0.354783);
    for (int i=0; i<2; ++i) {
      r[i]=(TH1D*)f[i]->Get(Form("fRatio_%.0f_%.0f",kCentBinsLimitsPion[iC][0],kCentBinsLimitsPion[iC][1]));
      //std::cout << Form("fRatio_%.0f_%.0f",kCentBinsLimitsPion[iC][0],kCentBinsLimitsPion[iC][1]) << std::endl;
      r[i]->GetXaxis()->SetRangeUser(0.7,1.6);
      r[i]->GetYaxis()->SetRangeUser(0.95,1.05);
      r[i]->Fit("pol0","Q");
      ratio[i]=r[i]->GetFunction("pol0")->GetParameter(0);
      ratio_err[i]=r[i]->GetFunction("pol0")->GetParError(0);
      r[i]->SetLineColor(colors[i]);
      r[i]->SetMarkerColor(colors[i]);
      r[i]->SetMarkerSize(1.2);
      r[i]->GetFunction("pol0")->SetLineColor(colors[i]);
      r[i]->GetFunction("pol0")->SetLineStyle(7);
      l.AddEntry(r[i],Form("LHC18%s",period[i]));
    }
    c.cd();
    r[0]->Draw(opt[0]);
    r[1]->Draw(opt[1]);
    t.DrawLatex(0.8,1.035,Form("Barlow check: |R_{q}-R_{r}|/#sqrt{#sigma_{q}^{2}+#sigma_{r}^{2}} = %.2f",std::abs(ratio[1]-ratio[0])/sqrt(ratio_err[1]*ratio_err[1]+ratio_err[0]*ratio_err[0])));
    l.Draw("same");
    f_out.cd();
    c.Write();
    h.SetBinContent(iC_+1,std::abs(ratio[1]-ratio[0])*0.5);
  }
  h.Write();
  f_out.Close();
}