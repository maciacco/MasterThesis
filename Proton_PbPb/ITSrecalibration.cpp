// ITS recalibration

#include <TFile.h>
#include <Riostream.h>
#include <TH3F.h>
#include <TH1F.h>
#include <TF1.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooAbsPdf.h>

#include "../utils/Config.h"
#include "../utils/Utils.h"
#include "../utils/RooGausExp.h"

using namespace proton;
using namespace utils;

const char *in_name[] = {"ProtonDataITS","../data/Proton_PbPb/AnalysisResults_low_pt_proton_mc"}; //../data/Proton_PbPb/AnalysisResults_
const char *out_name[] = {"Data","Mc"};
const char *ls_name[] = {"","mcFalse_"};

void ITSrecalibration(bool mc=false, bool save_fits=false){
  TFile *inFile=TFile::Open(Form("%s.root",in_name[(int)mc]));
  if (!inFile) std::cout << "Input file does not exist!" << std::endl;
  TFile outFile(Form("out/ITScalib%s.root",out_name[(int)mc]),"recreate");
  TH3F *fNSigmaITSpTCentA=nullptr, *fNSigmaITSpTCentM=nullptr;
  TString inFileName(in_name[(int)mc]);
  if (mc || !inFileName.CompareTo("../data/Proton_PbPb/AnalysisResults_")) {
    TTList *ls_proton=(TTList*)inFile->Get(Form("nuclei_proton_%s",ls_name[(int)mc]));
    fNSigmaITSpTCentA=(TH3F*)ls_proton->Get(Form("f%sITSnSigma",kAntimatterMatter[0]));
    fNSigmaITSpTCentM=(TH3F*)ls_proton->Get(Form("f%sITSnSigma",kAntimatterMatter[1]));
  }
  else {
    fNSigmaITSpTCentA=(TH3F*)inFile->Get(Form("nuclei_proton_/f%sITSnSigma",kAntimatterMatter[0]));
    fNSigmaITSpTCentM=(TH3F*)inFile->Get(Form("nuclei_proton_/f%sITSnSigma",kAntimatterMatter[1]));
  }
  TH1F fITSCalibMean("fITSCalibMean",";#it{p}_{T} (GeV/#it{c});Mean (a.u.)",18,0.3,1.2);
  TH1F fITSCalibSigma("fITSCalibSigma",";#it{p}_{T} (GeV/#it{c});Sigma (a.u.)",18,0.3,1.2);
  TH1F fITSCalibMeanRoo("fCustomITSpidMu",";#it{p}_{T} (GeV/#it{c});Mean (a.u.)",18,0.3,1.2);
  TH1F fITSCalibSigmaRoo("fCustomITSpidSigma",";#it{p}_{T} (GeV/#it{c});Sigma (a.u.)",18,0.3,1.2);
  TH1F fITSCalibRatio("fITSCalibRatio",";#it{p}_{T} (GeV/#it{c});Ratio (a.u.)",18,0.3,1.2);
  if (save_fits) outFile.mkdir("fits");
  for (int iB=1;iB<54;iB++){
    TH1D* fNSigmaITSpT=fNSigmaITSpTCentA->ProjectionZ(Form("proj_%.2f_%.2f",kPtBins[iB-1],kPtBins[iB]),1,10,iB,iB);
    TH1D* fNSigmaITSpTM=fNSigmaITSpTCentM->ProjectionZ(Form("proj_%.2f_%.2f",kPtBins[iB-1],kPtBins[iB]),1,10,iB,iB);
    fNSigmaITSpT->Add(fNSigmaITSpTM);
    fNSigmaITSpT->Rebin(1);
    fNSigmaITSpT->GetXaxis()->SetRangeUser(-5,-2);
    double max = fNSigmaITSpT->GetBinCenter(fNSigmaITSpT->GetMaximumBin());
    fNSigmaITSpT->GetXaxis()->SetRangeUser(-6,6);
    double upNsigma = 4.;
    double lowNsigma = max+0.8;
    //if (kPtBins[iB-1]>0.39 && kPtBins[iB-1]<0.6) lowNsigma = -3.2;
    /* else if (kPtBins[iB-1]>0.6 && kPtBins[iB-1]<0.89) lowNsigma = -2.7;
    else if (kPtBins[iB-1]>0.89 && kPtBins[iB-1]<0.94) lowNsigma = -2.5;
    else if (kPtBins[iB-1]>0.94) lowNsigma = -2.4;
    if (kPtBins[iB-1]>0.99) lowNsigma = -2.; */
    RooRealVar nSigma("nSigma","nSigma",lowNsigma,upNsigma);
    RooDataHist data("data","data",RooArgList(nSigma),fNSigmaITSpT);
    RooPlot *xframe = nSigma.frame();
    RooRealVar mu("mu","mu",-0.5,-1.,1.);
    RooRealVar sig("sig","sig",0.5,0.,1.);
    RooRealVar tau("tau","tau",0.8,0.6,2.);
    RooRealVar muB("muB","muB",-6.,-1.);
    RooRealVar sigB("sigB","sigB",0.5,0.,1.);
    RooRealVar tauB("tauB","tauB",0.5,0.,2.);
    RooRealVar slope("slope","slope",-10,10);
    RooGausExp signal("signal","signal",nSigma,mu,sig,tau);
    //RooGausExp bkg("bkg","bkg",nSigma,muB,sigB,tauB);
    RooExponential bkg("bkg","bkg",nSigma,slope);
    RooRealVar w("w","w",0.,1.);
    RooAddPdf model("model","model",signal,bkg,w);
    model.fitTo(data);
    data.plotOn(xframe,RooFit::Name("data"));
    model.plotOn(xframe,RooFit::Name("model"));
    model.plotOn(xframe,RooFit::Components("signal"),RooFit::Name("signal"),RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
    model.plotOn(xframe,RooFit::Components("bkg"),RooFit::Name("bkg"),RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
    model.paramOn(xframe,RooFit::Label(TString::Format("#chi^{2}/NDF = %2.2f", xframe->chiSquare("model", "data"))), RooFit::Layout(0.58812,0.911028,0.861955));
    TF1 func("func","gaus(0)+expo(3)");
    func.SetParLimits(0,0,1.e9);
    func.SetParLimits(1,-1.5,0.);
    func.SetParLimits(2,0.5,2.);
    func.SetParameter(1,-0.7);
    func.SetParameter(2,0.7);
    func.SetParLimits(4,-20.,10.);
    if (kPtBins[iB-1]<0.5){for (int i=0;i<2;++i)fNSigmaITSpT->Fit("func","RQL","",-3.5,0.2);}
    else if (kPtBins[iB-1]<0.6){for (int i=0;i<2;++i)fNSigmaITSpT->Fit("func","RQL","",-3.,0.2);}
    else{for (int i=0;i<2;++i)fNSigmaITSpT->Fit("func","RQL","",-2.5,0.2);}
    //if (fNSigmaITSpT->GetFunction("func")->GetChisquare()/fNSigmaITSpT->GetFunction("func")->GetNDF() > 3.) continue;
    fITSCalibMean.SetBinContent(iB,fNSigmaITSpT->GetFunction("func")->GetParameter(1));
    fITSCalibMean.SetBinError(iB,fNSigmaITSpT->GetFunction("func")->GetParError(1));
    fITSCalibSigma.SetBinContent(iB,fNSigmaITSpT->GetFunction("func")->GetParameter(2));
    fITSCalibSigma.SetBinError(iB,fNSigmaITSpT->GetFunction("func")->GetParError(2));
    fITSCalibMeanRoo.SetBinContent(iB,mu.getVal());
    fITSCalibMeanRoo.SetBinError(iB,mu.getError());
    fITSCalibSigmaRoo.SetBinContent(iB,sig.getVal());
    fITSCalibSigmaRoo.SetBinError(iB,sig.getError());
    fITSCalibRatio.SetBinContent(iB,mu.getVal()/sig.getVal());
    fITSCalibRatio.SetBinError(iB,mu.getVal()/sig.getVal()*TMath::Sqrt(mu.getError()*mu.getError()/mu.getVal()/mu.getVal()+sig.getError()*sig.getError()/sig.getVal()/sig.getVal()));
    if (save_fits){
      outFile.cd("fits");
      fNSigmaITSpT->Write();
      xframe->Write();
    }
  }
  outFile.cd("");
  //fITSCalibMean.Write();
  //fITSCalibSigma.Write();
  fITSCalibMeanRoo.Write();
  fITSCalibSigmaRoo.Write();
  //fITSCalibRatio.Write();
}