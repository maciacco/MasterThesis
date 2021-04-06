// CutTree.cpp
// downscale TTree

#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>

const char *kPathToFile = "../data/Hypertriton_PbPb";

void CutTree(const char *inFileName, const char *outFileName, const char *treeName, const double fraction)
{
  // read tree
  ROOT::RDataFrame df(treeName, Form("%s/%s.root", kPathToFile, inFileName));

  // filter data data frame
  auto df_filtered = df.Filter(Form("( (gRandom->Rndm())>%f )", fraction));

  // write
  df_filtered.Snapshot(treeName, Form("./%s.root", outFileName));
}

void CutTree(const char *inFileName, const char *outFileName, const char *treeName){
  // read tree
  ROOT::RDataFrame df(treeName, Form("%s/%s.root", kPathToFile, inFileName));

  // get subset of columns
  df.Snapshot(treeName, Form("./%s.root", outFileName), {"V0CosPA", "pt", "ProngsDCA", "PiProngPvDCAXY", "He3ProngPvDCAXY", "He3ProngPvDCA", "PiProngPvDCA", "NpidClustersHe3", "TPCnSigmaHe3", "TPCnSigmaPi", "NitsClustersHe3", "Matter", "centrality", "ct", "m"});

}
