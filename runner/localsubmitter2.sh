#!/bin/sh

# a=(7074 7076 7078 7081 7082 7084 7087 7089 7092 7094 7096 7099 7104 7106 7107 7110 7121 7123 7125 7126 7128 7129 7131 7132 7133 7134 7135 7136 7137 7138 7139 7140 7141 7142 7143 7144)
# b=(7075 7075 7077 7080 7080 7083 7086 7088 7090 7093 7095 7098 7103 7105 7105 7105 7120 7122 7127 7127 7127 7127 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130)

# a=(7074 7076 7078 7081 7082 7084 7087 7089 7092 7096 7099 7104 7106 7107 7110 7121 7123 7125 7126 7128 7129 7131 7132 7133 7135 7136 7137 7140 7141 7142 7143 7144)
# b=(7075 7075 7077 7080 7080 7083 7086 7088 7090 7095 7098 7103 7105 7105 7105 7120 7122 7127 7127 7127 7127 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130)

# a=(7106 7107 7121 7123 7125 7126 7128 7129 7131 7132 7133 7135 7136 7137 7140 7141 7142 7143 7144)
# b=(7105 7105 7120 7122 7127 7127 7127 7127 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130 7130)

a=(7074 7076 7078 7081 7082 7084 7087 7089 7092 7096 7099 7104)
b=(7075 7075 7077 7080 7080 7083 7086 7088 7090 7095 7098 7103)

c=(7075 7077 7080 7083 7086 7088 7090 7093 7095 7098 7101 7103 7105 7120 7122 7127 7130)

d=(50 51 52 57 59 60 61 68 70 71 74 75 76 77 78 79 80 81 82 83 84 85 86)

# #TPC/PMT baselines
# for(( i = 0 ; i < 17 ; i++ ))
# do
# # 	bsub -q 8nh -J run_0_$i "batchsubmit.sh 0 ${c[$i]}"	#GetTPCBaseline
# 	bsub -q 8nh -J run_6_$i "batchsubmit.sh 6 ${c[$i]}"	#GetPMTBaseline
# done

# # #PMT calibration
# for(( i = 0 ; i <= 22 ; i++ ))
# do
# 	bsub -q 8nh -J run_9_$i "batchsubmit.sh 9 ${c[$i]}"	#GetPMTCalibration2
# done


#others
for(( i = 0 ; i <= 18 ; i++ ))
do
	echo  ${a[$i]} ${b[$i]}
# 	bsub -q 8nh -J run_1_$i "batchsubmit.sh 1 ${a[$i]} ${b[$i]}"	#WriteTPCWF
# 	bsub -q 8nh -J run_2_$i "batchsubmit.sh 2 ${a[$i]} ${b[$i]}"	#FindHits
# 	bsub -q 1nh -J run_3_$i "batchsubmit.sh 3 ${a[$i]}"		#PlotHitsSummary
# 	bsub -q 1nh -J run_4_$i "batchsubmit.sh 4 ${a[$i]}"		#PlotHits
# 	bsub -q 1nh -J run_5_$i "batchsubmit.sh 5 ${a[$i]}"		#FitTracks
	
# 	bsub -q 8nh -J run_8_$i "batchsubmit.sh 8 ${a[$i]} ${b[$i]}"	#WritePMTWF
# 	bsub -q 8nh -J run_10_$i "batchsubmit.sh 10 ${a[$i]} ${b[$i]}"	#GetPMTIntegral
	
# 	bsub -q 8nh -J run_104_$i "batchsubmit.sh 104 ${a[$i]} ${b[$i]}"	#TrackAnalysis
# 	bsub -q 8nh -J run_105_$i "batchsubmit.sh 105 ${a[$i]}"	#MultiViewPlot
	
	
# 	./localsubmit.sh 5 ${a[$i]}
	./localsubmit.sh 104 ${a[$i]}
	wait
	
	
	
	
	
# #	./analyze 1 ${a[$i]} ${b[$i]}
# 	./analyze 2 ${a[$i]} ${b[$i]}
# #	./analyze 4 ${a[$i]}
# #	./analyze 101 ${a[$i]}
# 	./analyze 5 ${a[$i]}
done





