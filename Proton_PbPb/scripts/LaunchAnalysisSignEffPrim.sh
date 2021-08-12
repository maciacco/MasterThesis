# LaunchAnalysis.sh

# parameters
cutSettings="$1"
binCountingFlag="$2"
expFlag="$3" # 1->sum of 2 exp
extractRatios=1

fileData="AnalysisResults"
fileMC="mc"
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
SignalBinned("$cutSettings",$argumentSignal,"$fileData","$signalName","update")
//Efficiency("$cutSettings","$fileMC","$EfficiencyHe3")
//Secondary("$cutSettings","$fileData","$fileMC","$PrimaryHe3")
.q
EOF
fi
