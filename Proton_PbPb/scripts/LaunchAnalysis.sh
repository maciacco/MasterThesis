# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=1
expFlag=1 # 1->sum of 2 exp, 0 -> sum of exp and pol
sigmoidFlag=0
spectraHistNameId=""
extractRatios=1

fileData="AnalysisResults"
fileMC="AnalysisResults_LHC21l5_full_largeDCA_cutChi2"
signalName="SignalProtonGausDExpSignal1_LongMCTracks_1"
spectraName="SpectraProton_MC21l5_raw_primary"
spectraNameTPC="SpectraProtonTPC_MC21l5_raw_primary"
EfficiencyHe3="EfficiencyProtonMC_21l5_false__"
PrimaryHe3="PrimaryProton_large"
PrimaryHe3TPC="PrimaryProtonTPC_large"

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

root -b -l <<EOF
.L SecondaryTPC.cpp+
SecondaryTPC("$cutSettings",0.12,"$fileData","$fileMC","$PrimaryHe3TPC",false)
EOF

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
SignalBinned("$cutSettings",8,$argumentSignal,"$fileData","$signalName","recreate")
Secondary("$cutSettings",0.12,"$fileData","$fileMC","$PrimaryHe3",false)
Spectra("$cutSettings",8,0,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraName","recreate","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3")
//AbsorptionError("AbsErrorMCorrection","recreate","$spectraName")
.q
EOF
fi

root -b -l <<EOF
.L SpectraTPC.cpp+
SpectraTPC("$cutSettings",8,0,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","$spectraNameTPC","recreate","AnalysisResults","$signalName","$EfficiencyHe3","$PrimaryHe3TPC")
EOF
