# AnalysisSys.sh
# This scripts launchs the whole analysis for systematics computation

cutDCAz="$1"
cutTPCcls="$2"
cutDCAxy="$3"
binCountingFlag="$4"
expFlag="$5"
sigmoidFlag="$6"
spectraHistNameId="$7"

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

argumentCuts="$cutDCAz,$cutTPCcls,$cutDCAxy"
argumentSignal="$binCountingFlag,$expFlag"

root -b -l <<EOF
.L SignalUnbinned.cpp+
.L Efficiency.cpp+
.L EfficiencySec.cpp+
.L Secondary.cpp+
SignalUnbinned($argumentCuts,$argumentSignal,"TreeOutDataSys","SignalHe3Sys","update")
Efficiency($argumentCuts,"TreeOutMCSys")
EfficiencySec($argumentCuts,"TreeOutMCSys")
Secondary($argumentCuts,"TreeOutDataSys","TreeOutMCSys")
.q
EOF
