#include "../utils/Config.h"

float cent[11] = {0.f,5.f,10.f,20.f,30.f,40.f,50.f,60.f,70.f,80.f,90.f};
float pt[55]={0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5,2.55,2.6,2.65,2.7,2.75,2.8,2.85,2.9,2.95,3};
const char* split[]={"< 0.5", ">0.5"};

void HistFromTree() {
  ROOT::RDataFrame df("LambdaTree","AnalysisResults_PrFromV0.root");
  TFile outFile("ProtonDataITS.root","recreate");
  outFile.mkdir("nuclei_proton_");
  for (int iC = 0; iC < 2; ++iC) {
    const int fITSSigmaNbins = 255;
    const double fITSSigmaLimit = 5.;
    float ITSSigmaBins[fITSSigmaNbins + 1];
    const float deltaITSSigma = 2.f * fITSSigmaLimit / fITSSigmaNbins;
    for (int i = 0; i <= fITSSigmaNbins; ++i)
      ITSSigmaBins[i] = i * deltaITSSigma - fITSSigmaLimit;
    auto df_split = df.Filter(Form("matter %s && cosPA > 0.9999 && dcaPrPV < 1. && mass > 1.105 && mass < 1.13", split[iC]));
    auto h3 = df_split.Histo3D({Form("f%sITSnSigma",kAntimatterMatter[iC]),Form("f%sITSnSigma",kAntimatterMatter[0]),10,cent,54,pt,fITSSigmaNbins,ITSSigmaBins},"centrality","ptPr","itsNsigmaPr");
    outFile.cd("nuclei_proton_");
    h3->Write();
  }
  outFile.Close();
}