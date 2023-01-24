DATA_PATH="/data/mciacco/MatchingStudies/data"
MC_PATH="/data/mciacco/MatchingStudies/mc"
TREE="K0sTree"
SPLIT_PERIODS=0
MAKE_HIST_FROM_TREE=0
GET_EFF=0
COMBINE=1

mkdir results_mc_int_merge
mkdir results_data_int_merge

# data analysis for split magnetic field polarity periods

if [ $SPLIT_PERIODS -eq 1 ]; then

for PERIOD in q r
do

if [ $MAKE_HIST_FROM_TREE -eq 1 ]; then
root -b -l <<EOF
.x MakeHist.cpp("$MC_PATH/AnalysisResults_${PERIOD}_K0s.root", "out_${PERIOD}.root", "$TREE")
EOF
fi

if [ $GET_EFF -eq 1 ]; then
root -b -l <<EOF
.x Efficiency.cpp("out_${PERIOD}.root", "MatchingTPCTOF_${PERIOD}.root")
EOF
fi
done

if [ $COMBINE -eq 1 ]; then
root -b -l <<EOF
.x CombineResults.cpp
EOF
fi

fi

# data analysis for merged magnetic field polarities

if [ $SPLIT_PERIODS -eq 0 ]; then

if [ $MAKE_HIST_FROM_TREE -eq 1 ]; then
root -b -l <<EOF
.x MakeHist.cpp("$DATA_PATH/AnalysisResults_q_K0s.root", "out_qr.root", "$TREE", "$DATA_PATH/AnalysisResults_r_K0s.root")
EOF
fi

if [ $GET_EFF -eq 1 ]; then
root -b -l <<EOF
.x Efficiency.cpp("out_qr.root", "MatchingTPCTOF_qr.root")
EOF
fi

if [ $COMBINE -eq 1 ]; then
root -b -l <<EOF
.x CombineResults.cpp("qr")
EOF
fi

fi