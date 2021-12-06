// CutTree.cpp
// downscale TTree

#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>

const char *kPathToFile = ".";

void CutTreeLambda(const char *inFileName, const char *outFileName, const char *treeName, const double fraction)
{
  // read tree
  ROOT::RDataFrame df(treeName, Form("%s/%s.root", kPathToFile, inFileName));

  // filter data data frame
  auto df_filtered = df.Filter(Form("(gRandom->Rndm())<%f && pt > 0", fraction));

  // write
  df_filtered.Snapshot(treeName, Form("./%s.root", outFileName));
}

void CutTreeLambda(const char *inFileName, const char *outFileName, const char *treeName){
  // read tree
  ROOT::RDataFrame df(treeName, Form("%s/%s.root", kPathToFile, inFileName));

  // get subset of columns
  df.Snapshot(treeName, Form("./%s.root", outFileName), {"cosPA", "pt", "cosPA", "dcaV0tracks", "dcaPrPV", "dcaPiPV", "dcaV0PV" , "tpcNsigmaPr", "matter", "centrality", "ct", "mass", "isReconstructed", "ctMC", "ptMC"});

}
