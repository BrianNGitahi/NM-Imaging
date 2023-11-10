# This script performs the whole simulation in one go: neuron-->fluorescence signal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

# bring the necessary functions
from s_functions import simulate_neuron, simulate_nm_conc, plot_nm_conc,simulate_fluorescence_signal, plot_f_signal

# FUTURE UPDATE: Take as arguments the timesteps and rate -- from the command line 

# STEP1: simulate a neuron
timesteps = 70000
rate = 13
firing_neuron = simulate_neuron(n_timesteps=timesteps,firing_rate=rate)

# check exactly how many spikes were produced: to see if it works
n_spikes = np.size(np.nonzero(firing_neuron))

# print simulated neuron summary:
print('Simulated neuron with {} spikes in {} timesteps ({} Hz).'.format(n_spikes, timesteps, rate))

# STEP2: simulate nm dynamics
nm_conc_in, nm_b_conc_in, nm_r_conc_in = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)
print('Simulated nm dynamics from neuron')
plot_nm_conc(nm_conc_in, nm_b_conc_in, nm_r_conc_in)

# FUTURE UPDATE: take in as arguments from the command line the bleach time constant for the dye and nm
bleach_time = 10e5

# STEP3: simulate fluorescence signal from these dynamics and plot it
progression_in, progression_sub_in = simulate_fluorescence_signal(tau_d=bleach_time, tau_nm=bleach_time, tau_tissue=10e9, nm_conc=nm_conc_in)
plot_f_signal(progression_in, progression_sub_in, nm_conc_in)
print('Simulated fluorescence signal from nm dynamics')

