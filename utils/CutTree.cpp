// CutTree.cpp
// downscale TTree

#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>

const char *kPathToFile="../data/Hypertriton_PbPb";

void CutTree(const char *inFileName, const char *outFileName, const char *treeName){
  // read tree
  ROOT::RDataFrame df(treeName,Form("%s/%s.root",kPathToFile,inFileName));

  // filter data data frame
  auto df_filtered=df.Filter("( (gRandom->Rndm())>0.78 && centrality<10) || ( (gRandom->Rndm())>0.4 && centrality>30)" );

  // write
  df_filtered.Snapshot(treeName,Form("./%s.root",outFileName));
}
