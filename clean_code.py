from matplotlib import pyplot as plt
from pandas import *
import numpy as np
from rockit import *
from casadi import *

data = read_csv("leuven_october2022_16-22.csv")
# Time dependent parameters
time = np.asarray(data['time'].tolist())
time = np.reshape(time,[len(time),1])

temp = np.asarray(data['temp'].tolist()) # deg C
for i in range(len(temp)):
    temp[i] = temp[i] + 273.15  # deg K
temp = np.reshape(temp,[len(temp),1])

Qsun = np.asarray(data['solrad'].tolist())  # W/m2
Qsun = np.reshape(Qsun,[len(Qsun),1])

historical = np.append(time, temp, axis=1)
historical = np.append(historical, Qsun, axis=1)



