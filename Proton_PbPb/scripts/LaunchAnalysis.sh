# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=1
expFlag=1 # 1->sum of 2 exp, 0 -> sum of exp and pol
sigmoidFlag=1
spectraHistNameId=""
extractRatios=1

fileData="AnalysisResults_15o"
fileMC="mc_15Injected"
signalName="SignalProtonGausDExpSignal1_15"
spectraName="SpectraProtonGausDExpSignal1_Corrected_noMB_15"
EfficiencyHe3="EfficiencyProton_15"
PrimaryHe3="PrimaryProton_15"

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
.L Efficiency.cpp+
.L Secondary.cpp+
.L Spectra.cpp+
//SignalBinned("$cutSettings",$argumentSignal,"$fileData","$signalName","recreate")
//Efficiency("$cutSettings","$fileMC","$EfficiencyHe3")
Secondary("$cutSettings","$fileData","$fileMC","$PrimaryHe3")
Spectra("$cutSettings",$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","recreate","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
.q
EOF
fi
