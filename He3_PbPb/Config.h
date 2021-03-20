// Config.h
// configuration file

#ifndef CONFIGFILE_H
#define CONFIGFILE_H

//////////////////////////////////////////////////////////////
// He3 analysis
//////////////////////////////////////////////////////////////

#define PID // enable pid studies

// directories
const char *kDataDir = "../data/He3_PbPb";
const char *kResDir = "./results";
const char *kOutDir = "./out";
const char *kPlotDir = "./plots";

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
const int kNCentClasses = 3;                                        // number of analysed centrality classes in He3 analysis
const int kCentBinsHe3[][2] = {{1, 1}, {2, 3}, {6, 7}};             // centrality classes bin indexes in He3 analysis
const double kCentBinsLimitsHe3[][2] = {{0, 5}, {5, 10}, {30, 50}}; // centrality classes bin limits in He3 analysis

// TPC nsigma binning
const int kNSigmaBins = 240;
const double kDeltaNSigma = 0.05f;
const double kLowestNSigma = -6.f;
const double kNSigmaMin = -6.f;
const double kNSigmaMax = 6.f;

// DCAxy binning
const int kNDCABins = 38;
const int kNDCABinsLarge = 26;
const double kDCABins[kNDCABins + 1] = {-1.30, -1.20, -1.10, -1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.07, -0.04, -0.02, 0.00, 0.02, 0.04, 0.07, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30};
const double kDCABinsLarge[kNDCABinsLarge + 1] = {-1.30, -1.10, -0.90, -0.70, -0.50, -0.40, -0.30, -0.20, -0.15, -0.10, -0.07, -0.04, -0.02, 0.00, 0.02, 0.04, 0.07, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.90, 1.10, 1.30};
// antimatter/matter
const char *kAntimatterMatter[2] = {"A", "M"};
const char *kAntimatterMatterLabel[2] = {"^{3}#bar{He}", "^{3}He"};

// axis labels
const char *kAxisTitleCent = "Centrality (%)";
const char *kAxisTitlePt = "#it{p}_{T} (GeV/#it{c})";
const char *kAxisTitleDCA = "DCA_{xy} (cm)";
const char *kAxisTitleNSigma = "^{3}He n#sigma (a.u.)";

// track selections
const char *kTrackSelectionsEta = "(std::abs(eta)<0.8)";
const char *kTrackSelectionsDCAxy = "std::abs(dcaxy)<0.1f";

// cuts with variations (systematics computation)
const int kNCutDCAz = 12;
const double kCutDCAz[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
const int kNCutTPCClusters = 31;
const int kCutTPCClusters[] = {74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104};
const char *kBoolString[] = {"false", "true"};
const int kNAnalyses = 2976;
const int kNPoints = 1e4;

#endif // CONFIGFILE_H