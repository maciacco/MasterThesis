# AnalysisSys.sh
# This scripts launchs the whole analysis for systematics computation

cutDCAz="$1"
cutTPCcls="$2"
cutDCAxy="$3"
cutChi2TPC="$4"
binCountingFlag="$5"
expFlag="$6"
sigmoidFlag="$7"
spectraHistNameId="$8"

# create output directories
DIR_OUT=out
DIR_RES=results

if [ -d "$DIR_OUT" ]; then
    echo "Directory '$DIR_OUT' already exists"
else
    mkdir $DIR_OUT
fi

if [ -d "$DIR_RES" ]; then
    echo "Directory '$DIR_RES' already exists"
else
    mkdir $DIR_RES
fi

# launch analysis

argumentCuts="$cutDCAz,$cutTPCcls,$cutDCAxy,$cutChi2TPC"
argumentSignal="$binCountingFlag,$expFlag"

root -b -l <<EOF
.L Spectra.cpp+
Spectra($argumentCuts,$argumentSignal,$sigmoidFlag,"$spectraHistNameId","SpectraHe3Syst_extend2","update","AnalysisResults","SignalHe3Sys_extend2")
.q
EOF
