#!/bin/sh

filename='../RunList.txt'
filenameB='../BaselineRuns.txt'

declare -a A
declare -a B
declare -a C

while IFS=" " read xx yy;do
	A+=($xx)
	B+=($yy)
done < $filename

while IFS=" " read xx;do
	C+=($xx)
done < $filenameB

# # list the runs
# for(( i = 0 ; i < ${#A[@]} ; i++ ))
# do
# 	echo $i ${A[$i]} ${B[$i]}
# done

if [ $# -eq 1 ]
then
	if [ $1 -eq 9 ]
	then
		for(( i = 49 ; i < 140 ; i++ ))
		do
			echo $1 $i
				
			echo "executable = batchsubmit.sh" > submit.sh
			echo "arguments = $1 $i" >> submit.sh
			echo "log = log.\$(Cluster).txt" >> submit.sh
			echo "output = output.\$(Cluster).txt" >> submit.sh
			echo "error = error.\$(Cluster).txt" >> submit.sh
			echo "+JobFlavour = \"workday\"" >> submit.sh
			echo "getenv = True" >> submit.sh
			echo "queue" >> submit.sh
			
			condor_submit submit.sh
		done
	else
		for(( i = 0 ; i < ${#A[@]} ; i++ ))
	# 	for(( i = 0 ; i < ${#C[@]} ; i++ ))
	# 	for(( i = 0 ; i < 106 ; i++ ))
		do
			echo "$1 ${A[$i]} ${B[$i]}"
			
			echo "executable = batchsubmit.sh" > submit.sh
			echo "arguments = $1 ${A[$i]} ${B[$i]}" >> submit.sh
			echo "log = log.\$(Cluster).txt" >> submit.sh
			echo "output = output.\$(Cluster).txt" >> submit.sh
			echo "error = error.\$(Cluster).txt" >> submit.sh
			echo "+JobFlavour = \"workday\"" >> submit.sh
			echo "getenv = True" >> submit.sh
			echo "queue" >> submit.sh
			
			condor_submit submit.sh
		
		
# # # 		baselines
# 		echo "$1 ${C[$i]}"
# 		
# 		echo "executable = batchsubmit.sh" > submit.sh
# 		echo "arguments = $1 ${C[$i]}" >> submit.sh
# 		echo "log = log.\$(Cluster).txt" >> submit.sh
# 		echo "output = output.\$(Cluster).txt" >> submit.sh
# 		echo "error = error.\$(Cluster).txt" >> submit.sh
# 		echo "+JobFlavour = \"workday\"" >> submit.sh
# 		echo "getenv = True" >> submit.sh
# 		echo "queue" >> submit.sh
# 		
# 		condor_submit submit.sh
		
		done
	fi
else
	if [ $1 -eq 0 ]
	then
		echo $1 $2
		
		echo "executable = batchsubmit.sh" > submit.sh
		echo "arguments = $1 $2" >> submit.sh
		echo "log = log.\$(Cluster).txt" >> submit.sh
		echo "output = output.\$(Cluster).txt" >> submit.sh
		echo "error = error.\$(Cluster).txt" >> submit.sh
		echo "+JobFlavour = \"workday\"" >> submit.sh
		echo "getenv = True" >> submit.sh
		echo "queue" >> submit.sh
		
		condor_submit submit.sh
	else
		for(( i = 0 ; i < ${#A[@]} ; i++ ))
		do
			if [ ${A[$i]} -eq $2 ]
			then
				echo $1 ${A[$i]} ${B[$i]}
				
				echo "executable = batchsubmit.sh" > submit.sh
				echo "arguments = $1 ${A[$i]} ${B[$i]}" >> submit.sh
				echo "log = log.\$(Cluster).txt" >> submit.sh
				echo "output = output.\$(Cluster).txt" >> submit.sh
				echo "error = error.\$(Cluster).txt" >> submit.sh
				echo "+JobFlavour = \"workday\"" >> submit.sh
				echo "getenv = True" >> submit.sh
				echo "queue" >> submit.sh
				
				condor_submit submit.sh
			fi
		done
	fi
fi






# if [ $1 -gt -10 ]
# then
# 	echo "executable = batchsubmit.sh" > submit.sh
#         echo "arguments = $1 $2 $3" >> submit.sh
#         echo "log = log.\$(Cluster).txt" >> submit.sh
#         echo "output = output.\$(Cluster).txt" >> submit.sh
#         echo "error = error.\$(Cluster).txt" >> submit.sh
# # 	echo "+JobFlavour = \"longlunch\"" >> submit.sh
#        echo "+JobFlavour = \"workday\"" >> submit.sh
#         echo "getenv = True" >> submit.sh
#         echo "queue" >> submit.sh
# 
#         condor_submit submit.sh
# fi

# f=(7274 7276 7278 7280 7282 7284 7286)

# for(( i = 0 ; i < ${#a[@]} ; i++ ))
# do
# 	echo "$i ${a[$i]} ${b[$i]}"
# done

# for(( i = 0 ; i < ${#b[@]} ; i++ ))
# do
# 	echo "$i ${b[$i]}"
# done

# # 	echo "executable = batchsubmit.sh" > submit.sh
# # 	echo "arguments = $1 $2 $3" >> submit.sh
# # 	echo "log = log.\$(Cluster).txt" >> submit.sh
# #         echo "output = output.\$(Cluster).txt" >> submit.sh
# #         echo "error = error.\$(Cluster).txt" >> submit.sh
# # #	echo "+JobFlavour = \"longlunch\"" >> submit.sh
# #  	echo "+JobFlavour = \"workday\"" >> submit.sh
# # 	echo "getenv = True" >> submit.sh
# # 	echo "queue" >> submit.sh
# # 	
# #         condor_submit submit.sh


# for(( i = 36 ; i < ${#a[@]} ; i++ ))
# for i in 49
# do
# 	echo ${a[$i]}
# 	echo "executable = batchsubmit.sh" > submit.sh
#  	echo "arguments = 1 ${a[$i]} ${b[$i]}" >> submit.sh
# # 	echo "arguments = 0 ${b[$i]}" >> submit.sh
# 	echo "log = log.\$(Cluster).txt" >> submit.sh
#         echo "output = output.\$(Cluster).txt" >> submit.sh
#         echo "error = error.\$(Cluster).txt" >> submit.sh
# # 	echo "+JobFlavour = \"longlunch\"" >> submit.sh
# 	echo "+JobFlavour = \"workday\"" >> submit.sh
# 	echo "getenv = True" >> submit.sh
# 	echo "queue" >> submit.sh
# 
#         condor_submit submit.sh
# done
# 
#baselines
# for(( i = 0 ; i < ${#c[@]} ; i++ ))
# do
# 	echo "executable = batchsubmit.sh" > submit.sh
# 	echo "arguments = 0 ${c[$i]}" >> submit.sh
# 	echo "log = log.\$(Cluster).txt" >> submit.sh
#         echo "output = output.\$(Cluster).txt" >> submit.sh
#         echo "error = error.\$(Cluster).txt" >> submit.sh
# 	echo "+JobFlavour = \"longlunch\"" >> submit.sh
# # 	echo "+JobFlavour = \"workday\"" >> submit.sh
# 	echo "getenv = True" >> submit.sh
# 	echo "queue" >> submit.sh
# 	
#         condor_submit submit.sh
# done


# # #TPC/PMT baselines
# # for(( i = 0 ; i < 17 ; i++ ))
# for(( i = 0 ; i < 7 ; i++ ))
# do
# # 	bsub -q 8nh -J run_0_$i "batchsubmit.sh 0 ${k[$i]}"	#GetTPCBaseline
# #	bsub -q 8nh -J run_1_$i "batchsubmit.sh 1 ${k[$i]} ${k[$i]}"
# 	bsub -q 8nh -J run_2_$i "batchsubmit.sh 2 ${k[$i]} ${k[$i]}"
# # 	bsub -q 8nh -J run_6_$i "batchsubmit.sh 6 ${c[$i]}"	#GetPMTBaseline
# done

# # #PMT calibration
# for(( i = 0 ; i <= 22 ; i++ ))
# do
# 	bsub -q 8nh -J run_9_$i "batchsubmit.sh 9 ${c[$i]}"	#GetPMTCalibration2
# done

#k=(5176 5177 5178 5179 5180 5181 5182 5197 5198 5200 5201 5202 5207)

#l=(5208 5209 5210 5211 5212 5219 5220 5221 5222 5223 5224 5225 5226 5227 5228 5229 5230 5231 5232 5233)

#others
# for(( i = 0 ; i <= 35 ; i++ ))
# do
# 	echo  ${a[$i]} ${b[$i]}
# # 	bsub -q 8nh -J run_1_$i "batchsubmit.sh 1 ${a[$i]} ${b[$i]}"	#WriteTPCWF
# # 	bsub -q 8nh -J run_2_$i "batchsubmit.sh 2 ${a[$i]} ${b[$i]}"	#FindHits
# # 	bsub -q 1nh -J run_3_$i "batchsubmit.sh 3 ${a[$i]}"		#PlotHitsSummary
# # 	bsub -q 1nh -J run_4_$i "batchsubmit.sh 4 ${a[$i]}"		#PlotHits
# # 	bsub -q 1nh -J run_5_$i "batchsubmit.sh 5 ${a[$i]}"		#FitTracks
# 	
# # 	bsub -q 8nh -J run_8_$i "batchsubmit.sh 8 ${a[$i]} ${b[$i]}"	#WritePMTWF
# # 	bsub -q 8nh -J run_10_$i "batchsubmit.sh 10 ${a[$i]} ${b[$i]}"	#GetPMTIntegral
# 	
#  	bsub -q 1nh -J run_104_$i "batchsubmit.sh 104 ${a[$i]} ${b[$i]}"	#TrackAnalysis
# #	bsub -q 8nh -J run_105_$i "batchsubmit.sh 105 ${a[$i]}"	#MultiViewPlot
# 	
# 	
# 	
# # #	./analyze 1 ${a[$i]} ${b[$i]}
# # 	./analyze 2 ${a[$i]} ${b[$i]}
# # #	./analyze 4 ${a[$i]}
# # #	./analyze 101 ${a[$i]}
# # 	./analyze 5 ${a[$i]}
# done





