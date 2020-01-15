#!/bin/sh

A=( 208 210 213 217 223 228 231 236 241 245 248 252 256 260 264 273 280 286 295 297 302 306 310 315 321 329 338 )

for(( i = 0 ; i < ${#A[@]} ; i++ ))
do
	echo $1 ${A[$i]}
	
	echo "executable = PMTcalibsubmit.sh" > submit.sh
	echo "arguments = $1 ${A[$i]}" >> submit.sh
	echo "log = log.\$(Cluster).txt" >> submit.sh
	echo "output = output.\$(Cluster).txt" >> submit.sh
	echo "error = error.\$(Cluster).txt" >> submit.sh
	echo "+JobFlavour = \"workday\"" >> submit.sh
	echo "getenv = True" >> submit.sh
	echo "queue" >> submit.sh
	
	condor_submit submit.sh
done





