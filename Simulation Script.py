# This script performs the whole simulation in one go: neuron-->fluorescence signal


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit

# bring the necessary functions
from s_functions import simulate_neuron, simulate_nm_conc, simulate_fluorescence_signal

# simulate a neuron
timesteps = 70000
rate = 13
firing_neuron = simulate_neuron(n_timesteps=timesteps,firing_rate=rate)

# check exactly how many spikes were produced: to see if it works
n_spikes = np.size(np.nonzero(firing_neuron))

# print simulated neuron summary:
print('Simulated neuron with {} spikes in {} timesteps ({} Hz).'.format(n_spikes, timesteps, rate))

