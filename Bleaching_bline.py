# This is an identical script to the bleachign analysis
# The only difference is how we calculate the df/f values
# Here we use a moving baseline: F0 is the average of x previous f values 

# import necessary pacakges
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

# import the necessary functions
from Functions import simulate_neuron, simulate_nm_conc


# define the background fluorescence due to the tissue
bg_tissue = 1.5

# ANALYSIS 1: bleaching factor acting on the [NM] contribution to F
def bleach_1(K_D, tau, F_max, F_min, nm_conc, bline_len=5000):

    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # bleaching factor -- starts off as 1 then exponentially decreases 
    # we set tau to be a very large constant so this is a slow decrease
    bleach = np.exp(-t/tau)

    # calculate bleached f values: derived in part from eq 2 in Neher/Augustine
    f = bg_tissue + (K_D*F_min + bleach*nm_conc*F_max)/(K_D + bleach*nm_conc)

    # define the signal array
    delta_ft_f0 = np.zeros(f.size)

    # calculate f0 values and populate the signal array
    for i in range(f.size):

        # calculate f0 by averaging the previous x number of f values
        # if x is bigger than the current index then use all the prev f values
        # where x is the length of the moving baseline

        if i < bline_len:
            f0 = np.average(f[:i])
        else: 
            f0 = np.average(f[i-bline_len:i])

        # calculate normalized signal using the calculated f0
        delta_ft_f0[i] = (f[i]-f0)/(f0)


    # plot the normalized signal delta f/f0 at the different t
    plt.plot(t,delta_ft_f0, label = tau)
    plt.xlabel('time(ms)')
    plt.ylabel('Delta F/F0')
    plt.title('Flourescence intensity signal over time (bleach 1)')
    plt.legend()

    return delta_ft_f0


# check this bleaching effect 1 for different values of tau -- different bleach factor
different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])

# generate an input nm conc array from firing neuron
firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# plot bleached signal with different time constants for the bleach factor
for i in range(len(different_taus)):
   bleach_1(K_D = 1000, tau=different_taus[i], F_max = 45, F_min = 10, nm_conc=nm_conc, bline_len=5000)
plt.show()


# ANALYSIS 2: bleaching factor acting on the contributions of the dye, F0 and [NM]
def bleach_2(K_D, tau, F_max, F_min, nm_conc, bline_len =5000):

    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # bleaching factor -- starts off as 1 then exponentially decreases 
    # we set tau to be a very large constant so this is a slow decrease
    bleach = np.exp(-t/tau)

    # calculate bleached f values: derived in part from eq 2 in Neher/Augustine
    f= bg_tissue+ bleach*(K_D*F_min + nm_conc*F_max)/(K_D + nm_conc)
    
    # define the signal array
    delta_ft_f0 = np.zeros(f.size)

    # calculate f0 values and populate the signal array
    for i in range(f.size):

        # calculate f0 by averaging the previous x number of f values
        # if x is bigger than the current index then use all the prev f values
        # where x is the length of the moving baseline

        if i < bline_len:
            f0 = np.average(f[:i])
        else: 
            f0 = np.average(f[i-bline_len:i])

        # calculate normalized signal using the calculated f0
        delta_ft_f0[i] = (f[i]-f0)/(f0)

    # plot the normalized signal delta f/f0 at the different t
    plt.plot(t,delta_ft_f0, label = tau)
    plt.xlabel('time(ms)')
    plt.ylabel('Delta F/F0')
    plt.title('Flourescence intensity signal over time (bleach 2)')
    plt.legend()

    return delta_ft_f0


# check this bleaching effect 1 for different values of tau -- different bleach factor
different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])

# generate an input nm conc array from firing neuron
firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# plot bleached signal with different time constants for the bleach factor
for i in range(len(different_taus)):
   bleach_2(K_D = 1000, tau=different_taus[i], F_max = 45, F_min = 10, nm_conc=nm_conc)
plt.show()


# ANALYIS 3: bleaching factor acting on the background fluores
# cence 
def bleach_3(K_D, tau, F_max, F_min, nm_conc, bline_len=5000):

    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # bleaching factor -- starts off as 1 then exponentially decreases 
    # we set tau to be a very large constant so this is a slow decrease
    bleach = np.exp(-t/tau)

    # calculate bleached f values: derived in part from eq 2 in Neher/Augustine
    f= bleach*bg_tissue + (K_D*F_min + nm_conc*F_max)/(K_D + nm_conc)
    
   # define the signal array
    delta_ft_f0 = np.zeros(f.size)

    # calculate f0 values and populate the signal array
    for i in range(f.size):

        # calculate f0 by averaging the previous x number of f values
        # if x is bigger than the current index then use all the prev f values
        # where x is the length of the moving baseline

        if i < bline_len:
            f0 = np.average(f[:i])
        else: 
            f0 = np.average(f[i-bline_len:i])

        # calculate normalized signal using the calculated f0
        delta_ft_f0[i] = (f[i]-f0)/(f0)

    # plot the normalized signal delta f/f0 at the different t
    plt.plot(t,delta_ft_f0, label = tau)
    plt.xlabel('time(ms)')
    plt.ylabel('Delta F/F0')
    plt.title('Flourescence intensity signal over time (bleach 3)')
    plt.legend()

    return delta_ft_f0


# check this bleaching effect 1 for different values of tau -- different bleach factor
different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])

# generate an input nm conc array from firing neuron
firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# plot bleached signal with different time constants for the bleach factor
for i in range(len(different_taus)):
   bleach_3(K_D = 1000, tau=different_taus[i], F_max = 45, F_min = 10, nm_conc=nm_conc)
plt.show()