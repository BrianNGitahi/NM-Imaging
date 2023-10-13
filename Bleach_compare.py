# This script is used to compare the bleaching effects established in the 
# Bleaching Analysis and bleaching_bline scripts


# import necessary packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

# import functions from the simulation and bleaching libraries
from s_functions import simulate_neuron, simulate_nm_conc, simulate_flourescence_signal
from b_factory import bleach_nm, bleach_dnm, bleach_t, bleach_all


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


# create a fit for the initial values of the bleach 4 result and subtract it -- revise this later
# it should be a fit of the f value

poly4 = np.polyfit(t,b4,2)
fit4 = np.polyval(poly4,t)

plt.show()



# ANALYSIS 2: comparing how the different bleaches affect the signal to noise ratio

# define the signal to noise ratio

# get the signal t



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


