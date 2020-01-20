import numpy
from matplotlib.pyplot import figure, show

# read the saved ASCII data
chan1=numpy.loadtxt("../waveform/waveforms_ch0_run454.txt")
chan2=numpy.loadtxt("../waveform/waveforms_ch1_run454.txt")

tsz=len(chan1) 
print (tsz) # should be "1000001"

t=(1/25e6)*numpy.arange(0,tsz,1)
print(t)