#!/bin/sh

mkdir run_$1
cd run_$1

if [ $2 -eq 1 ]
then
	cp /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/NTupler/NTupler .
	wait

	./NTupler $1 >> /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/logs/AllLogs.txt
	wait
elif [ $2 -eq 2 ]
then
	cp /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/NTupler/PMTNTupler .
	wait

	./PMTNTupler $1 >> /afs/cern.ch/work/b/bilki/public/Flic/flicanalyzer/logs/AllLogs.txt
	wait
fi

cd ..
rm -r run_$1
