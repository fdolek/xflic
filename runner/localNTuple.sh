#!/bin/sh

mkdir run_$1
cd run_$1

# cp /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/NTupler/NTupler .
# wait
# 
# ./NTupler $1
# wait

cp /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/NTupler/PMTNTupler .
wait

./PMTNTupler $1
wait


cd ..
# rm -r run_$1
