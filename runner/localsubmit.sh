#!/bin/sh

date
mkdir run_$2
cd run_$2

# cp /afs/cern.ch/work/f/flic2019/flicanalyzer/analyze.C .
# wait
# cp /afs/cern.ch/work/f/flic2019/flicanalyzer/Makefile .
# wait
# make
# wait

cp /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/analyze .
wait
cp /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/RunList.txt .
wait
cp /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/BaselineRuns.txt .
wait


./analyze $1 $2 $3
wait

cd ..
#rm -r run_$2

date