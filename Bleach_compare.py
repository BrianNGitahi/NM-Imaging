# This script is used to compare the bleaching effects established in the 
# Bleaching Analysis and bleaching_bline scripts


# import necessary packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

# import functions from the simulation and bleaching libraries
from s_functions import simulate_neuron, simulate_nm_conc, simulate_fluorescence_signal
from b_factory import bleach_nm, bleach_dnm, bleach_t, bleach_all, bleach_dnm_heat


# ANALYSIS 1: 

# choose the timescale for the bleach factor
chosen_tau = 10e4

# generate an input nm conc array from firing neuron
firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# plot bleached signal for the bleach factor acting on different components of f

# create timesteps array for the plot
t = np.linspace(0,nm_conc.size-1,nm_conc.size)

b1 = bleach_nm(K_D = 1000, tau=chosen_tau, F_max = 45, F_min = 10, nm_conc=nm_conc, bline_len=5000)
b2 = bleach_dnm(K_D = 1000, tau=chosen_tau, F_max = 45, F_min = 10, nm_conc=nm_conc, bline_len=5000)
b3 = bleach_t(K_D = 1000, tau=chosen_tau, F_max = 45, F_min = 10, nm_conc=nm_conc, bline_len=5000)
b4 = bleach_all(K_D = 1000, tau=chosen_tau, F_max = 45, F_min = 10, nm_conc=nm_conc, bline_len=5000)



#ANALYSIS 2:

# check the effect of changing the variance of the gaussian noise on ftissue on the snr 
var_values = np.array([0.0001,0.001,0.01,0.1,1,3])
var_v = np.array([1,3,10])

# bleach time constants for heatmap
specific_taus = np.logspace(5,7,20)

# generate a firing neuron
neuron = simulate_neuron(n_timesteps=70000,firing_rate=13)

# generate nm_conc 
nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# for different variances, get the heatmap
plt.figure()
for i in range(len(var_v)):
    plt.subplot(3,2,i+1)
    bleach_dnm_heat(specific_taus,nm_conc_input=nm_conc, var = var_v[i])

plt.suptitle('SNR vs bleach strength at different variance for ftissue', size = 16)
plt.tight_layout()
plt.show()






# ANALYSIS 3:
# check this bleaching effect 1 for different values of tau -- different bleach factor
#different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])




# PREVIOUS ANALYSIS: comparing the different contributions at one timescale for the bleach factor
# plt.plot(t,b1,label='bleach 1')
# plt.plot(t,b2,label='bleach 2')

# plt.plot(t,b3,label='bleach 3')
# #plt.plot(t,b4,label='bleach 4')
# plt.plot(b4-fit4,label='subtracted 4')

# plt.xlabel('time(ms)')
# plt.ylabel('Delta F/F0')
# plt.title('Flourescence intensity signal over time')
# plt.legend()


