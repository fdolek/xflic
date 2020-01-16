import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib inline
#with open("waveforms_ch0_run421.txt") as m:
#with open("waveforms_ch0_run421.txt") as f:
#    print(f.readline())
    
    
df = pd.read_csv("waveforms_ch2_run443.txt",
                 delim_whitespace=True,
                 usecols=[1,],
                 names=["t",])
df.plot()