# LaunchAnalysis.sh

# parameters
cutSettings=""
binCountingFlag=1
expFlag=1 # 1->sum of 2 exp, 0 -> sum of exp and pol
sigmoidFlag=1
spectraHistNameId=""
extractRatios=1

fileData="mc_20g7_likeData"
fileMC="mc"
signalName="SignalProtonGausDExpSignal1_LongMCTracks_MClikeData"
spectraName="SpectraProtonGausDExpSignal1_LongMCTracks_newPrimary_MClikeData"
EfficiencyHe3="EfficiencyProton_LongMCTracks_new"
PrimaryHe3="PrimaryProtonMC"

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
.L Spectra.cpp+
//SignalBinnedMC("$cutSettings",$argumentSignal,"$fileData","$signalName","recreate")
//SecondaryMC("$cutSettings","$fileData","$fileMC","$PrimaryHe3")
EfficiencyNew("$cutSettings","$fileMC","$EfficiencyHe3","$signalName","$PrimaryHe3")
.q
EOF
fi
