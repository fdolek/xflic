import matplotlib.pyplot as plt
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks
import numpy as np
import pandas as pd

df = pd.read_csv("../waveform/waveforms_ch0_run454.txt",
                 delim_whitespace=True,
                 usecols=[1,])
df.to_numpy()

