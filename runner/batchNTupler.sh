#!/bin/sh

#second argument: 1- full NTupler 1 - PMTNTupler

k=$1
l=$2

a=(7074 7075 7076 7077 7078 7079 7080 7081 7082 7083 7084 7085 7086 7087 7088 7089 7090 7091 7092 7093 7094 7095 7096 7098 7099 7100 7101 7102 7103 7104 7105 7106 7107 7108 7109 7110 7121 7122 7123 7124 7125 7126 7127 7128 7129 7130 7131 7132 7133 7134 7135 7136 7137 7138 7139 7140 7141 7142 7143 7144 7240 7241 7242 7243 7244 7245 7246 7247 7248 7249 7250 7251 7252 7253 7254 7255 7256 7257 7258 7259 7260 7261 7262 7263 7264 7265 7266 7267 7268 7269 7270)

b=(50 51 52 57 59 60 61 68 70 71 74 75 76 77 78 79 80 81 82 83 84 85 86 88 91 92 93 94 95 96 97 98 99 100 101 102)

# for(( i = 0 ; i < ${#a[@]} ; i++ ))
# # for(( i = 0 ; i < 1 ; i++ ))
# do
# 	echo "executable = batchNTuple.sh" > submit.sh
#         echo "arguments = ${a[$i]} 1" >> submit.sh
#         echo "log = NTuple.log.\$(Cluster).txt" >> submit.sh
#         echo "output = NTuple.output.\$(Cluster).txt" >> submit.sh
#         echo "error = NTuple.error.\$(Cluster).txt" >> submit.sh
# # 	echo "+JobFlavour = \"longlunch\"" >> submit.sh
#        echo "+JobFlavour = \"workday\"" >> submit.sh
#         echo "getenv = True" >> submit.sh
#         echo "queue" >> submit.sh
# 
#         condor_submit submit.sh
# done

# #PMTNtupler
# for(( i = 0 ; i < ${#b[@]} ; i++ ))
# do
# #	bsub -q 8nh -J run_2_$i "batchNTuple.sh ${b[$i]} 2"
# 
# 	echo "executable = batchNTuple.sh" > submit.sh
#         echo "arguments = ${b[$i]} 2" >> submit.sh
#         echo "log = PMTlog.txt" >> submit.sh
#         echo "output = PMToutput.txt" >> submit.sh
#         echo "error = PMTerror.txt" >> submit.sh
#         echo "+JobFlavour = \"longlunch\"" >> submit.sh
# #         echo "+JobFlavour = \"workday\"" >> submit.sh
#         echo "getenv = True" >> submit.sh
#         echo "queue" >> submit.sh
# 
#         condor_submit submit.sh
# 
# done



if [ $k -gt 0 ]
then
	echo "executable = batchNTuple.sh" > submit.sh
        echo "arguments = $k $l" >> submit.sh
        echo "log = log.\$(Cluster).txt" >> submit.sh
#         echo "output = output.txt" >> submit.sh
        echo "output = output.\$(Cluster).txt" >> submit.sh
        echo "error = error.\$(Cluster).txt" >> submit.sh
# 	echo "+JobFlavour = \"longlunch\"" >> submit.sh
       echo "+JobFlavour = \"workday\"" >> submit.sh
        echo "getenv = True" >> submit.sh
        echo "queue" >> submit.sh

        condor_submit submit.sh
fi














