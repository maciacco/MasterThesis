// Config.h
// configuration file

#ifndef CONFIGFILE_H
#define CONFIGFILE_H

#define PID // enable pid studies

//////////////////////////////////////////////////////////////
// Common
//////////////////////////////////////////////////////////////

// directories
const char *kResDir = "./results";
const char *kOutDir = "./out";
const char *kPlotDir = "./plots";

// antimatter / matter
const char *kAntimatterMatter[2] = {"A", "M"};

// centrality binning
const int kNCentClasses = 3; // 3

const char *kAxisTitleDCA = "DCA_{xy} (cm)";
const char *kAxisTitlePt = "#it{p}_{T} (GeV/#it{c})";

const char *kBoolString[] = {"false", "true"};
const int kNPoints = 1e4;

//////////////////////////////////////////////////////////////
// He3 analysis
//////////////////////////////////////////////////////////////

namespace he3
{
  // directories
  const char *kDataDir = "../data/He3_PbPb";

  // histograms
  const char *kTPC = "TPCcounts";

  // data files
  const int kNDataFiles = 2;
  const char *kDataFileLabel[kNDataFiles] = {"LHC18q", "LHC18r"};

  // pt binning
  const int kNPtBins = 14; // analysis binning
  const double kDeltaPt = 0.5f;
  const double kLowestPt = 1.f;
  const int kNPtBinsFine = 36; // finer pt binning
  const double kDeltaPtFine = 0.25f;
  const int kNPt[2] = {12, 9};

  // centrality binning
  const int kNCentBins = 11; // total number of V0M multiplicity classes
  const double kCentBins[kNCentBins + 1] = {0.f, 5.f, 7.5f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f};
  const int kCentBinsHe3[][2] = {{1, 1}, {2, 3}, {6, 7}};             // centrality classes bin indexes in He3 analysis
  const double kCentBinsLimitsHe3[][2] = {{0, 5}, {5, 10}, {30, 50}}; // centrality classes bin limits in He3 analysis

  // TPC nsigma binning
  const int kNSigmaBins = 240;
  const double kDeltaNSigma = 0.05f;
  const double kLowestNSigma = -6.f;
  const double kNSigmaMin = -4.5f;
  const double kNSigmaMax = 4.5f;
  // DCAxy binning
  const int kNDCABins = 40;
  const int kNDCABinsLarge = 18;
  const int kNDCABinsLarge2 = 14;
  const double kDCABins[kNDCABins + 1] = {-1.30, -1.20, -1.10, -1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.07, -0.05, -0.04, -0.02, 0.00, 0.02, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30};
  const double kDCABinsLarge[kNDCABinsLarge + 1] = {-1.30, -1.10, -0.90, -0.70, -0.50, -0.35, -0.20, -0.10, -0.05, 0.00, 0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90, 1.10, 1.30};
  const double kDCABinsLarge2[kNDCABinsLarge2 + 1] = {-1.30, -1.10, -0.90, -0.70, -0.50, -0.30, -0.10, 0.00, 0.10, 0.30, 0.50, 0.70, 0.90, 1.10, 1.30};
  const int kNDCABinsMedium = 22;
  const double kDCABinsMedium[kNDCABinsMedium + 1] = {-1.30, -1.10, -0.90, -0.70, -0.50, -0.35, -0.20, -0.10, -0.07, -0.04, -0.02, 0.00, 0.02, 0.04, 0.07, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90, 1.10, 1.30};

  // antimatter/matter
  const char *kAntimatterMatterLabel[2] = {"^{3}#bar{He}", "^{3}He"};

  // axis labels
  const char *kAxisTitleCent = "Centrality (%)";
  const char *kAxisTitleNSigma = "^{3}He n#sigma (a.u.)";

  // track selections
  const char *kTrackSelectionsEta = "(std::abs(eta)<0.8)";
  const char *kTrackSelectionsDCAxy = "std::abs(dcaxy)<0.1f";

  // cuts with variations (systematics computation)
  const int kNCutDCAz = 11;
  const double kCutDCAz[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
  const int kNCutTPCClusters = 31;
  const int kCutTPCClusters[] = {74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104};

  const int kNAnalyses = 2976;
}

//////////////////////////////////////////////////////////////
// deuteron analysis
//////////////////////////////////////////////////////////////

namespace deuteron
{
  // directories
  const char *kDataDir = "../data/Deuteron_PbPb";

  // TOF signal binning
  const double kTOFSignalMin = -2.;
  const double kTOFSignalMax = 2.5;

  // pt binning
  const int kNPtBins = 24; // analysis binning
  double kPtBins[25] = {0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.4f, 1.6f, 1.8f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.2f, 4.4f, 5.0f, 6.0f, 8.0f};

  // centrality binning
  const double kCentBinsLimitsDeuteron[][2] = {{0, 5}, {5, 10}, {30, 50}, {0, 90}};
  const int kCentBinsDeuteron[][2] = {{1, 1}, {2, 2}, {5, 6}, {1, 10}}; // centrality classes bin indexes in He3 analysis

  // antimatter / matter
  const char *kAntimatterMatterLabel[2] = {"#bar{d}", "d"};
  const char *kAntimatterMatterLabelExtended[2] = {"antideuterons", "deuterons"};

  // DCAxy binning
  const int kNDCABins = 26;
  const int kNDCABinsLarge = 10;
  const double kDCABins[kNDCABins + 1] = {-1.30f, -1.00f, -0.80f, -0.60f, -0.50f, -0.30f, -0.20f, -0.12f, -0.07f, -0.06f, -0.05f, -0.04f, -0.02f, 0.00f, 0.02f, 0.04f, 0.05f, 0.06f, 0.07f, 0.12f, 0.20f, 0.30f, 0.50f, 0.60f, 0.80f, 1.00f, 1.30f};
  const double kDCABinsLarge[kNDCABinsLarge + 1] = {-1.30f, -1.00f, -0.80f, -0.40f, -0.12f, 0.00f, 0.12f, 0.40f, 0.80f, 1.00f, 1.30f};
  const int kNDCABinsMedium = 16;
  const double kDCABinsMedium[kNDCABinsMedium + 1] = {-1.30f, -1.00f, -0.80f, -0.50f, -0.20f, -0.12f, -0.06f, -0.03f, 0.00f, 0.03f, 0.06f, 0.12f, 0.20f, 0.50f, 0.80f, 1.00f, 1.30f};
  const int kNDCABinsMedium2 = 14;
  const double kDCABinsMedium2[kNDCABinsMedium2 + 1] = {-1.30f, -1.00f, -0.80f, -0.50f, -0.20f, -0.12f, -0.06f, 0.00f, 0.06f, 0.12f, 0.20f, 0.50f, 0.80f, 1.00f, 1.30f};
  /* const int kNDCABinsMedium = 26;
  const double kDCABinsMedium[kNDCABinsMedium + 1] = {-1.30f, -1.10f, -0.90f, -0.70f, -0.50f, -0.40f, -0.30f, -0.20f, -0.15f, -0.10f, -0.07f, -0.04f, -0.02f, 0.00f, 0.02f, 0.04f, 0.07f, 0.10f, 0.15f, 0.20f, 0.30f, 0.40f, 0.50f, 0.70f, 0.90f, 1.10f, 1.30f};
 */

  // systematics variations
  const int kNCutDCAz = 4;
  const double kCutDCAz[] = {0.5, 0.75, 1.5, 2.0};
  const int kNCutTPCClusters = 4;
  const double kCutTPCClusters[] = {59., 64., 74., 79.};
  const int kNTPCPidSigmas = 2;
  const double kTPCPidSigmas[] = {3.25, 3.50};
  const char *cutSettings[] = {"dcaz", "pid", "tpc"};
}

//////////////////////////////////////////////////////////////
// proton analysis
//////////////////////////////////////////////////////////////

namespace proton
{
  const bool ADD20g7 = false;

  // directories
  const char *kDataDir = "../data/Proton_PbPb";

  // TOF signal binning
  const double kTOFnSigmaMin = -20.;
  const double kTOFnSigmaMax = 20.;

  // pt binning
  //const int kNPtBins = 84; // analysis train binning
  //const int kNPtBins = 2;
  const int kNPtBins = 32;
  //const int kNPtBins = 70; // analysis binning
  //const int kNPtBins = 11;
  //double kPtBins[kNPtBins+1] = {1.00f, 1.50f, 2.00f};
  double kPtBins[85]={0.8f, 0.85f, 0.9f, 0.95f, 1.00f, 1.05f, 1.1f, 1.15f, 1.2f, 1.25f, 1.3f, 1.35f, 1.4f, 1.45f, 1.5f, 1.55f, 1.6f, 1.65f, 1.7f, 1.75f, 1.8f, 1.85f, 1.9f, 1.95f, 2.00f, 2.05f, 2.1f, 2.15f, 2.2f, 2.25f, 2.3f, 2.35f, 2.4f, 2.45f, 2.5f, 2.55f, 2.6f, 2.65f, 2.7f, 2.75f, 2.8f, 2.85f, 2.9f, 2.95f, 3.00f, 3.05f, 3.1f, 3.15f, 3.2f, 3.25f, 3.3f, 3.35f, 3.4f, 3.45f, 3.5f, 3.55f, 3.6f, 3.65f, 3.7f, 3.75f, 3.8f, 3.85f, 3.9f, 3.95f, 4.00f, 4.05f, 4.1f, 4.15f, 4.2f, 4.25f, 4.3f, 4.35f, 4.4f, 4.45f, 4.5f, 4.55f, 4.6f, 4.65f, 4.7f, 4.75f, 4.8f, 4.85f, 4.9f, 4.95f, 5.00f}; // binning from trains
  //double kPtBins[45]={0.8f, 0.85f, 0.9f, 0.95f, 1.00f, 1.05f, 1.1f, 1.15f, 1.2f, 1.25f, 1.3f, 1.35f, 1.4f, 1.45f, 1.5f, 1.55f, 1.6f, 1.65f, 1.7f, 1.75f, 1.8f, 1.85f, 1.9f, 1.95f, 2.00f, 2.1f, 2.2f, 2.3f,2.4f, 2.5f, 2.6f, 2.7f, 2.8f, 2.9f, 3.00f, 3.2f, 3.4f, 3.6f, 3.8f, 4.00f, 4.2f, 4.4f, 4.6f, 4.8f, 5.00f};
  //double kPtBins[33]={0.8f, 0.85f, 0.9f, 0.95f, 1.00f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.00f, 2.1f, 2.2f, 2.3f, 2.4f, 2.5f, 2.6f, 2.8f, 3.00f, 3.2f, 3.4f, 3.6f, 3.8f, 4.00f, 4.2f, 4.4f, 4.6f, 4.8f, 5.00f};

  // centrality binning
  const double kCentBinsLimitsProton[][2] = {{0, 5}, {5, 10}, {30, 50}, {0, 90}};
  const int kCentBinsProton[][2] = {{1, 1}, {2, 2}, {5, 6}, {1, 10}, {1,4}}; // centrality classes bin indexes in He3 analysis

  // antimatter / matter
  const char *kAntimatterMatterLabel[2] = {"#bar{p}", "p"};
  const char *kAntimatterMatterLabelExtended[2] = {"antiprotons", "protons"};

  // DCAxy binning
  const int kNDCABinsLarge = 34;
  const double kDCABinsLarge[kNDCABinsLarge + 1] = {-1.30f,-1.20f,-1.10f,-1.00f,-0.90f,-0.80f,-0.70f,-0.60f,-0.50f,-0.40f,
    -0.35f,-0.30f,-0.25f,-0.20f,-0.15f,-0.10f,-0.05f,0.00f,0.05f,0.10f,0.15f, 0.20f,
     0.25f, 0.30f, 0.35f, 0.40f, 0.50f, 0.60f, 0.70f, 0.80f, 0.90f, 1.00f,
     1.10f, 1.20f, 1.30f};
  const int kNDCABinsTask = 52;
  const int kNDCABinsMedium = /* 52; */26;
  const double kDCABinsMedium[kNDCABinsMedium + 1] = {-1.30f, -1.10f, -0.90f, -0.70f, -0.50f, -0.35f, -0.25f, -0.15f, -0.10f, -0.08f, -0.06f, -0.04f, -0.02f, -0.0f, 0.02f, 0.04f, 0.06f, 0.08f, 0.10f, 0.15f, 0.25f, 0.35f, 0.50f, 0.70f, 0.90f, 1.10f, 1.30f};
  const double kDCABinsTask[kNDCABinsTask + 1] =  {-1.30f,-1.20f,-1.10f,-1.00f,-0.90f,-0.80f,-0.70f,-0.60f,-0.50f,-0.40f,
    -0.35f,-0.30f,-0.25f,-0.20f,-0.15f,-0.12f,-0.10f,-0.09f,-0.08f,-0.07f,
    -0.06f,-0.05f,-0.04f,-0.03f,-0.02f,-0.01f, 0.00f, 0.01f, 0.02f, 0.03f,
     0.04f, 0.05f, 0.06f, 0.07f, 0.08f, 0.09f, 0.10f, 0.12f, 0.15f, 0.20f,
     0.25f, 0.30f, 0.35f, 0.40f, 0.50f, 0.60f, 0.70f, 0.80f, 0.90f, 1.00f,
     1.10f, 1.20f, 1.30f};
  
  /* const int kNDCABinsLarge = 26;
  const double kDCABinsLarge[kNDCABinsLarge + 1] = {-1.30f, -1.10f, -0.90f, -0.70f, -0.50f, -0.35f, -0.25f, -0.12f, -0.10f, -0.08f, -0.06f, -0.04f, -0.02f, -0.0f, 0.02f, 0.04f, 0.06f, 0.08f, 0.10f, 0.12f, 0.25f, 0.35f, 0.50f, 0.70f, 0.90f, 1.10f, 1.30f};
 *//* const int kNDCABinsMedium = 16;
  const double kDCABinsMedium[kNDCABinsMedium + 1] = {-1.30f, -1.00f, -0.80f, -0.50f, -0.20f, -0.12f, -0.06f, -0.03f, 0.00f, 0.03f, 0.06f, 0.12f, 0.20f, 0.50f, 0.80f, 1.00f, 1.30f};
   */
  // systematics variations
  const int kNCutDCAz = 5;
  const double kCutDCAz[] = {0.5, 0.75, 1.0, 1.5, 2.0};
  const int kNCutTPCClusters = 5;
  const double kCutTPCClusters[] = {59., 64., 69., 74., 79.};
  const int kNTPCPidSigmas = 3;
  const double kTPCPidSigmas[] = {3.0, 3.25, 3.50};
  const char *cutSettings[] = {"dcaz", "pid", "tpc"};
  const int kNTrackCuts = 11;
  const int trackCutIndexes[kNTrackCuts] = {0,0,1,2,3,0,1,0,1,2,3};
  const char* trackCutSettings[kNTrackCuts] = {"","dcaz","dcaz","dcaz","dcaz","pid","pid","tpc","tpc","tpc","tpc"};
}

#endif // CONFIGFILE_H
