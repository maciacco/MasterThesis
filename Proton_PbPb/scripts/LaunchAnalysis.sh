# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=1
expFlag=1 # 1->sum of 2 exp, 0 -> sum of exp and pol
sigmoidFlag=1
spectraHistNameId=""
extractRatios=1

fileData="AnalysisResults"
fileMC="mc"
signalName="SignalProtonGausDExpSignal1_LongMCTracks"
spectraName="SpectraProtonGausDExpSignal1_LongMCTracks_newPrimary"
EfficiencyHe3="EfficiencyProton_LongMCTracks_new"
PrimaryHe3="PrimaryProtong"

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
.L Secondary.cpp+
.L Spectra.cpp+
.L AbsorptionError.cpp+
//SignalBinned("$cutSettings",8,$argumentSignal,"$fileData","$signalName","recreate")
Secondary("$cutSettings","$fileData","$fileMC","$PrimaryHe3")//,true)
Spectra("$cutSettings",8,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","recreate","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
//AbsorptionError("AbsErrorMCorrection","recreate","$spectraName")
.q
EOF
fi
