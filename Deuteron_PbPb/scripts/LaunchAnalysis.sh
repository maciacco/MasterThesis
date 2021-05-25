# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=0
expFlag=1 # 1->sum of 2 exp
sigmoidFlag=1
spectraHistNameId="_1_1_1"
extractRatios=1

fileData="AnalysisResults"
signalName="SignalDeuteron"
spectraName="SpectraDeuteron"
EfficiencyHe3="EfficiencyDeuteron"
EfficiencyHe3SecWD="EfficiencyDeuteronSecWD"
PrimaryHe3="PrimaryDeuterons"

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
.L SignalBinned.cpp+
//.L Efficiency.cpp+
//.L EfficiencySec.cpp+
//.L Secondary.cpp+
//.L Spectra.cpp+
SignalBinned("$cutSettings",$argumentSignal,"$fileData","$signalName","recreate")
// Efficiency($argumentCuts,"$treeMC","$EfficiencyHe3")
// EfficiencySec($argumentCuts,"$treeMC","$EfficiencyHe3SecWD",$HYPER_TO_HE3_RATIO)
// Secondary($argumentCuts,"$treeData","$treeMC","$EfficiencyHe3SecWD","$PrimaryHe3")
// Spectra($argumentCuts,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","recreate","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
.q
EOF
fi
