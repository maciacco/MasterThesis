#include <Riostream.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TString.h>
#include <TCanvas.h>
#include <ROOT/RDataFrame.hxx>
#include "AliPWGFunc.h"

const char* kInFileMCName = "./AnalysisResults_reweight_BW_30_50.root";
const char* kInFileCentName = "../data/LambdaPrompt_PbPb/StrangenessRatios_summary.root";

constexpr bool reject = false;

int iPart = 1;
double speed_of_light = 2.99792458;
double xi_0_lifetime = 2.90;
double max_ct = 100;

int int_pow(int a, int b){
  if (b < 0) return -999;
  int r = 1;
  for (int i=0;i<b;++i)
    r *= a;
  return r;
}

void ctShaping(const char *inFileMCName=kInFileMCName, const char *inFileCentName=kInFileCentName){
  ROOT::EnableImplicitMT(4);
  TFile inFileCent(inFileCentName);
  ROOT::RDataFrame df("LambdaTree",inFileMCName);
  TH1D *hCent = (TH1D*)inFileCent.Get("Centrality_selected");
  int part_flag = int_pow(2,iPart);
  std::cout << "iPart = " << iPart << std::endl;
  // reweight trees
  std::cout << "Process tree..." << std::endl;
  std::string cut_variable="ctMotherMC";
  TF1 fExp("fExp","[0]*TMath::Exp(-[1]*x)",0,max_ct);
  fExp.SetParameter(1,1/xi_0_lifetime/speed_of_light);
  auto hNorm = df.Filter(Form("(flag & BIT(%d))==%d",iPart,part_flag)).Histo1D({"hNorm","hNorm",20,0,max_ct},cut_variable.data());
  hNorm->GetXaxis()->SetRangeUser(5,max_ct);
  double normMin = hNorm->GetMinimum();
  hNorm->GetXaxis()->SetRangeUser(0.,max_ct);
  TF1 fExpNorm("fExpNorm","[0]*exp(-[1]*x)",0,max_ct);
  TH1D *hTmp = new TH1D(*hNorm);
  fExpNorm.SetParLimits(0,0.,1.e8);
  fExpNorm.SetParLimits(1,0.,.5);
  hTmp->Fit("fExpNorm","MR+","",0,max_ct);
  fExp.SetParameter(0,fExpNorm.Eval(max_ct)/TMath::Exp(-max_ct/xi_0_lifetime/speed_of_light));
  std::cout << "Process tree (2)..." << std::endl;
  std::cout << "cut variable = " << cut_variable.data() << "; iPart = " << iPart << "; part_flag = " << part_flag << std::endl;
  if (reject){
    auto reweight = [fExpNorm, normMin, fExp](float ct){
      double hNormVal = fExpNorm.Eval(ct);
      double exp = fExp.Eval(ct);
      bool cut = (gRandom->Rndm()*hNormVal) < exp;
      return cut;
    };
    df.Filter(Form("(flag & BIT(%d))==%d && %s < %f",iPart,part_flag,cut_variable.data(),max_ct)).Filter(reweight,{cut_variable.data()}).Snapshot("LambdaTree",Form("AnalysisResults_%d.root",iPart));
  }
  else{
    auto reweight = [fExpNorm, normMin, fExp](float ct){
      double hNormVal = fExpNorm.Eval(ct);
      double exp = fExp.Eval(ct);
      double weight = exp/hNormVal;
      return weight;
    };
    df.Filter(Form("(flag & BIT(%d))==%d && %s < %f",iPart,part_flag,cut_variable.data(),max_ct)).Define("weightMC",reweight,{cut_variable.data()}).Snapshot("LambdaTree",Form("AnalysisResults_%d.root",iPart));
  }
  TCanvas canv("c","c");
  canv.SetLogy();
  fExpNorm.Draw();
  fExp.Draw("same");
  canv.Print("c.png");
  gSystem->Exec("mv AnalysisResults_1.root AnalysisResults_Xi0_reweight_BW_30_50.root");
}