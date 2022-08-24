#include <Riostream.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TString.h>
#include <ROOT/RDataFrame.hxx>

constexpr double kMass[] = {1.115683, 1.32171, 1.67245};
const char *partNames[] = {"Lambda", "Xi", "Omega"};

// BW parameters from https://arxiv.org/abs/1910.07678
constexpr int nCent = 10;
constexpr double cent[nCent+1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
constexpr double beta_avg[nCent] = {0.663, 0.660, 0.655, 0.643, 0.622, 0.595, 0.557, 0.506, 0.435, 0.355};
constexpr double T_kin[nCent] = {0.090, 0.091, 0.094, 0.097, 0.101, 0.108, 0.115, 0.129, 0.147, 0.161};
constexpr double n[nCent] = {0.735, 0.736, 0.739, 0.771, 0.828, 0.908, 1.052, 1.262, 1.678, 2.423};
constexpr double avg_mult[nCent] = {1943, 1587, 1180, 786, 512, 318, 183, 96.3, 44.9, 17.5};

constexpr int centSplit[][2]={{0,0},{1,1},{4,5},{0,10}}; // 0-5%, 5-10%, 30-50%, 0-90%
constexpr int nCentSplit[]={1,1,2,10}; // 0-5%, 5-10%, 30-50%, 0-90%
constexpr int centClass = 0; // 0-5% -> 0, 5-10% -> 1, 30-50% -> 2, 0-90% -> 3

constexpr bool reject = true;

TF1* fLastFunc;

Double_t IntegrandBG(const double * x, const double* p){
  // integrand for boltzman-gibbs blast wave
     // x[0] -> r (radius)
     // p[0] -> mass
     // p[1] -> pT (transverse momentum)
     // p[2] -> beta_max (surface velocity)
     // p[3] -> T (freezout temperature)
     // p[4] -> n (velocity profile)


  double x0 = x[0]; 
  
  double mass     = p[0];
  double pT       = p[1];
  double beta_max = p[2];
  double temp     = p[3];
  Double_t n      = p[4];

  // Keep beta within reasonable limits
  Double_t beta = beta_max * TMath::Power(x0, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;

  double mT      = TMath::Sqrt(mass*mass+pT*pT);

  double rho0   = TMath::ATanH(beta);  
  double arg00 = pT*TMath::SinH(rho0)/temp;
  if (arg00 > 700.) arg00 = 700.; // avoid FPE
  double arg01 = mT*TMath::CosH(rho0)/temp;
  double f0 = x0*mT*TMath::BesselI0(arg00)*TMath::BesselK1(arg01);

  //  printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", x0, pT, beta_max, temp, n, mT, beta, rho0, arg00, arg01);

  return f0;
}

Double_t StaticBGdNdPt(const double * x, const double* p) {

  // implementation of BGBW (1/pt dNdpt)

  double pT = x[0];;
  

  double mass    = p[0];
  double beta    = p[1];
  double temp    = p[2];
  double n       = p[3];
  double norm    = p[4];

  static TF1 * fIntBG = 0;
  if(!fIntBG)
    fIntBG = new TF1 ("fIntBG", IntegrandBG, 0, 1, 5);

  fIntBG->SetParameters(mass, pT, beta, temp,n);
  double result = fIntBG->Integral(0,1);
  //  printf ("[%4.4f], Int :%f\n", pT, result);
  return result*norm;//*1e30;;

}

Double_t StaticBGdNdPtTimesPt(const double * x, const double* p) {
  // BGBW dNdpt implementation
  return x[0]*StaticBGdNdPt(x,p);
}

TF1* GetBGBW(Double_t mass, Double_t beta, Double_t T,
			 Double_t n, Double_t norm, const char * name){

  // Boltzmann-Gibbs blast wave
  fLastFunc = new TF1 (name, StaticBGdNdPtTimesPt, 0.0, 10, 5);
  fLastFunc->SetParameters(mass,beta,T,n,norm);    
  fLastFunc->FixParameter(0,mass);
  fLastFunc->SetParNames("mass", "#beta", "temp", "n", "norm");
  return fLastFunc;
  
}

double beta_max(const double beta_avg, const double n){
  return beta_avg*0.5*(2+n);
}

int int_pow(int a, int b){
  if (b < 0) return -999;
  int r = 1;
  for (int i=0;i<b;++i)
    r *= a;
  return r;
}

double full_bw(TF1 *bw_array[], double pT){
  int sz = nCentSplit[centClass];
  double bw = 0;
  for (int iF=0; iF<sz; ++iF){
    bw+=bw_array[iF]->Eval(pT);
  }
  return bw;
}

double full_bw_max(TF1 *bw_array[], double max_steps = 1.e3, double max_pT = 10.){
  TH1D h("h","h",max_steps,0,max_pT);
  for (int i=0;i<max_steps;++i){
    h.SetBinContent(i+1,full_bw(bw_array,i/(max_steps/max_pT)));
  }
  return h.GetMaximum();
}

const char* kInFileMCName = "/data/mciacco/Omega_PbPb/AnalysisResults-mc.root";
const char* kInFileCentName = "/data/mciacco/LambdaPrompt_PbPb/StrangenessRatios_summary.root";
const char* kOutFileName = "pTShapesLXiOm_0_5.root";


void pTShaping(const char *inFileMCName=kInFileMCName, const char *inFileCentName=kInFileCentName, const char *outFileName=kOutFileName){
  ROOT::EnableImplicitMT(4);
  TFile inFileCent(inFileCentName);
  TFile outFile(outFileName, "recreate");
  ROOT::RDataFrame df("XiOmegaTree",inFileMCName);
  TH1D *hCent = (TH1D*)inFileCent.Get("Centrality_selected");
  double nEv = hCent->Integral(cent[centSplit[centClass][0]]+1,cent[centSplit[centClass][1]+1]);
  for (int iPart = 2; iPart < 3; ++iPart){
    int int_pow_part = int_pow(2,iPart);
    std::cout << "iPart = " << iPart << std::endl;
    // build BW
    TString dir(Form("bw_%s",partNames[iPart]));
    outFile.mkdir(dir);
    outFile.cd(dir);
    TF1 *blastWave[nCentSplit[centClass]];
    std::string name[nCentSplit[centClass]];
    for (int iC = centSplit[centClass][0]; iC <= centSplit[centClass][1]; ++iC){
      double nEvInRange = hCent->Integral(cent[iC]+1,cent[iC+1]);
      name[iC-centSplit[centClass][0]] = "bw_";
      name[iC-centSplit[centClass][0]] += std::to_string(iC);
      blastWave[iC-centSplit[centClass][0]] = GetBGBW(kMass[iPart],beta_max(beta_avg[iC],n[iC]),T_kin[iC],n[iC],1,name[iC-centSplit[centClass][0]].data());
      blastWave[iC-centSplit[centClass][0]] = GetBGBW(kMass[iPart],beta_max(beta_avg[iC],n[iC]),T_kin[iC],n[iC],avg_mult[iC]*nEvInRange/nEv/blastWave[iC-centSplit[centClass][0]]->Integral(0,1.e3),name[iC-centSplit[centClass][0]].data());
      std::cout << "weight = " << nEvInRange/nEv << std::endl;
      blastWave[iC-centSplit[centClass][0]]->Write();
    }
    TH1D *hBW=new TH1D("h","h",100000,0,10);
    for (int iB=1;iB<=100000;++iB){
      hBW->SetBinContent(iB,full_bw(blastWave,iB/10000.));
    }
    hBW->Write();
    double bw_max = hBW->GetMaximum();
    // reweight trees
    std::cout << "Process tree..." << std::endl;
    auto dff = df.Define("index","gRandom->Rndm()+mass-mass").Filter("index < 1");
    std::string cut_variable = "ptMC";
    auto hNorm = df.Filter("std::abs(pdg)==3334").Histo1D({"hNorm","hNorm",10000,0,10},cut_variable.data());
    hNorm->Write();
    hNorm->GetXaxis()->SetRangeUser(1.,4.);
    double normMin = hNorm->GetMinimum();
    hNorm->GetXaxis()->SetRangeUser(0.,10.);
    TH1D *hTmp = new TH1D(*hNorm);
    std::cout << "Process tree (2)..." << std::endl;
    std::cout << "cut variable = " << cut_variable.data() << "; iPart = " << iPart << "; int_pow_part = " << int_pow_part << std::endl;
    if (reject){
      auto reweight = [hTmp, normMin, bw_max, hBW](float pT){
        int sz = nCentSplit[centClass];
        double bw = hBW->GetBinContent(hBW->FindBin(pT));
        double hNormVal = hTmp->GetBinContent(hTmp->FindBin(pT));
        bool cut = gRandom->Rndm() < (bw*normMin/hNormVal/bw_max);
        return cut;
      };
      dff.Filter("std::abs(pdg)==3334").Filter(reweight,{cut_variable.data()}).Snapshot("XiOmegaTree","AnalysisResults_0_5.root");
    }
    else {
      auto reweight = [hTmp, normMin, bw_max, hBW](float pT){
        int sz = nCentSplit[centClass];
        double bw = hBW->GetBinContent(hBW->FindBin(pT));
        double hNormVal = hTmp->GetBinContent(hTmp->FindBin(pT));
        double weight=bw*normMin/hNormVal/bw_max;
        return weight;
      };
      dff.Filter("std::abs(pdg)==3334").Define("weightMC",reweight,{cut_variable.data()}).Snapshot("XiOmegaTree","AnalysisResults_0_5.root");
    }
  }

  outFile.Close();
}