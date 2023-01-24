#ifndef CONFIG_H
#define CONFIG_H

constexpr bool kDebug = true;
constexpr bool split_polarities = false;
constexpr const char* kReconstructed = "true";
constexpr const char* kOutDir = "results_data_int_merge";
constexpr const char* kOutDirCombine = "results_data_int_merge";
constexpr const char* kChargeFlag[2][2] = {{"Neg", "neg"}, {"Pos", "pos"}};
constexpr int kNCentBins = 5; // 18;
constexpr int kNPtBins = 40;
constexpr double kCentBinLimits[kNCentBins+1] = {0.,5.,10.,30.,50.,90.}; // {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.};
constexpr double kPtBinLimits[kNPtBins+1] = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2};
constexpr double kCutCosPA = 0.9995;
constexpr double kCutDcaV0PV = 0.5;
constexpr double kCutDcaV0tracks = 1.;
constexpr double kCutDCADaughPV = 5.;
constexpr double kCutTPCNsigmaDaugh = 3.;
constexpr double kCutRadius = 3.;
constexpr double kCutMass[] = {0.4985 - 0.0800, 0.4985 + 0.0800};

#endif // CONFIG_H