#!/bin/sh

for i in $@
do
	echo $i
	if [ $i -lt 7000 ]
	then
		scp icadaq@pcatlasst01.cern.ch:/home/icadaq/CAENV1751Daq/PMTDAQCode_updated/mcrbee10/waveforms_ch*_run${i}.txt /eos/project/f/flic2019/Data/PMT/WaveForms/
		wait
		echo "executable = batchNTuple.sh" > submit.sh
		echo "arguments = $i 2" >> submit.sh
		echo "log = log.\$(Cluster).txt" >> submit.sh
		echo "output = output.\$(Cluster).txt" >> submit.sh
		echo "error = error.\$(Cluster).txt" >> submit.sh
		echo "+JobFlavour = \"workday\"" >> submit.sh
		echo "getenv = True" >> submit.sh
		echo "queue" >> submit.sh
		
		condor_submit submit.sh
		wait
	else
		scp icadaq@pcatlasst01.cern.ch:/Data1/completed/Run00${i}.dat /eos/project/f/flic2019/Data/TPC/Runs/
		wait
		scp icadaq@pcatlasst01.cern.ch:/home/icadaq/CAENV1751Daq/PMTDAQCode_updated/mcrbee10/waveforms_ch*_run${i}.txt /eos/project/f/flic2019/Data/PMT/WaveForms/
		wait
		echo "executable = batchNTuple.sh" > submit.sh
		echo "arguments = $i 1" >> submit.sh
		echo "log = log.\$(Cluster).txt" >> submit.sh
		echo "output = output.\$(Cluster).txt" >> submit.sh
		echo "error = error.\$(Cluster).txt" >> submit.sh
		echo "+JobFlavour = \"workday\"" >> submit.sh
		echo "getenv = True" >> submit.sh
		echo "queue" >> submit.sh
		
		condor_submit submit.sh
		wait
	fi
done

