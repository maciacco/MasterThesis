# LaunchAnalysis.sh

# parameters
cutSettings="$1"
binCountingFlag="$2"
expFlag="$3" # 1->sum of 2 exp
roiNsigma="$4"
extractRatios=1

fileData="AnalysisResults"
fileDataEff="AnalysisResults_LHC21l5"
fileMC="mc"
signalNameEff="SignalProtonSysEff"
spectraNameEff="SpectraProtonSysEff"
signalName="SignalProtonSys"
spectraName="SpectraProtonSys"
EfficiencyHe3="EfficiencyProtonSys"
PrimaryHe3="PrimaryProtonSys"

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
if (($cutSettings==999)); then
    cutSettings=""
fi

echo $cutSettings
argumentSignal="$binCountingFlag,$expFlag"

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L ../utils/RooGausExp.cxx+
.L ../utils/RooDSCBShape.cxx+
.L ../utils/RooGausDExp.cxx+
.L SignalBinned.cpp+
.L Efficiency.cpp+
.L Secondary.cpp+
.L Spectra.cpp+
SignalBinned("$cutSettings",$roiNsigma,$argumentSignal,"$fileData","$signalName","update")
Secondary("$cutSettings","$fileData","$fileMC","$PrimaryHe3")
.q
EOF
fi

if [ $extractRatios -eq 1 ]; then
    root -b -l <<EOF
.L ../utils/RooGausExp.cxx+
.L ../utils/RooDSCBShape.cxx+
.L ../utils/RooGausDExp.cxx+
.L SignalBinnedMC.cpp+
.L EfficiencyNew.cpp+
.L SecondaryMC.cpp+
SignalBinnedMC("$cutSettings",$argumentSignal,"$fileDataEff","$signalNameEff","recreate")
SecondaryMC("$cutSettings","$fileDataEff","$fileMC","$PrimaryHe3Eff")
EfficiencyNew("$cutSettings","$fileMC","$EfficiencyHe3","$signalNameEff","$PrimaryHe3Eff")
.q
EOF
fi