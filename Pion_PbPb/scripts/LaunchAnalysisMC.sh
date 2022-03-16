# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=1
expFlag=1 # 1->sum of 2 exp, 0 -> sum of exp and pol
sigmoidFlag=0
spectraHistNameId=""
extractRatios=1

# fileData="../AnalysisResults_LHC22b9_1"
# fileMC="../AnalysisResults_LHC22b9_1"
# fileData="AnalysisResults_LHC21l5_full_largeDCA_cutChi2"
# fileMC="AnalysisResults_LHC21l5_full_largeDCA_cutChi2"
fileData="../../data/AnalysisResults_LHC22b9_2"
fileMC="../../data/AnalysisResults_LHC22b9_2"
# fileData="LHC20e3a"
# fileMC="LHC20e3a"
signalName="SignalPionMC_21l5_false__LHC22b9_2"
spectraName="SpectraPionMC_21l5_false_LHC22b9_2"
EfficiencyHe3="EfficiencyPion_LHC22b9_2"
PrimaryHe3="PrimaryPionMC_21l5_false_LHC22b9_2"

# create output directories
DIR_OUT=out
DIR_PLOT=plots

if [ -d "$DIR_OUT" ]; then
    echo "Directory '$DIR_OUT' already exists"
else
    mkdir $DIR_OUT
fi

if [ -d "$DIR_PLOT" ]; then
    echo "Directory '$DIR_PLOT' already exists"
else
    mkdir $DIR_PLOT
fi

# launch analysis

echo $cutSettings
argumentSignal="$binCountingFlag,$expFlag"

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L ../utils/RooGausExp.cxx+
.L ../utils/RooDSCBShape.cxx+
.L ../utils/RooGausDExp.cxx+
.L SignalBinnedMCv2.cpp+
.L SecondaryMC.cpp+
.L EfficiencyNew.cpp+
.L Efficiency.cpp+
//.L Spectra.cpp+
SignalBinnedMCv2("$cutSettings",1.5,11.,8.5,13.5,$argumentSignal,"$fileData","$signalName","recreate")
SecondaryMC("$cutSettings",0.12,"$fileData","$fileMC","$PrimaryHe3")
EfficiencyNew("$cutSettings","$fileMC","$EfficiencyHe3","$signalName","$PrimaryHe3")
//Efficiency("$cutSettings","$fileMC","$EfficiencyHe3")
.q
EOF
fi
