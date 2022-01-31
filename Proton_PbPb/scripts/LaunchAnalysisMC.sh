# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=1
expFlag=1 # 1->sum of 2 exp, 0 -> sum of exp and pol
sigmoidFlag=1
spectraHistNameId=""
extractRatios=1

fileData="AnalysisResults_LHC21l5_full_largeDCA"
fileMC="AnalysisResults_LHC21l5_full_largeDCA"
signalName="SignalProtonMC_21l5_false"
spectraName="SpectraProtonMCMC_21l5_false"
EfficiencyHe3="EfficiencyProtonMC_21l5_false_"
PrimaryHe3="PrimaryProtonMC_21l5_false"

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
.L SignalBinnedMC.cpp+
.L SecondaryMC.cpp+
.L EfficiencyNew.cpp+
//.L Spectra.cpp+
//SignalBinnedMC("$cutSettings",8.,$argumentSignal,"$fileData","$signalName","recreate")
//SecondaryMC("$cutSettings",0.12,"$fileData","$fileMC","$PrimaryHe3")
EfficiencyNew("$cutSettings","$fileMC","$EfficiencyHe3","$signalName","$PrimaryHe3")
.q
EOF
fi
