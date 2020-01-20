import numpy
from matplotlib.pyplot import figure, show
import pandas as pd

# read the saved ASCII data
df = pd.read_csv("Cosmic.dat",
                 delim_whitespace=True,)
                 #usecols=[1,])
                 #names=["y",])
df