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
from s_functions import simulate_neuron, simulate_nm_conc, simulate_fluorescence_signal


# define the background fluorescence due to the tissue
bg_tissue = 1.5


# functino to generate heatmap
def bleach_dnm_heat(tau_values, nm_conc_input, var):

    # create an array to store the snr values
    snr = np.zeros((tau_values.size, tau_values.size))

    # BONUS: find way to do it that's more effecient -- this is O(n^2)
    for i in range(len(tau_values)):
        for j in range(len(tau_values)):
            progression, progression_sub = simulate_fluorescence_signal(nm_conc=nm_conc_input,variance=var, tau_d=tau_values[i],tau_nm=tau_values[i], tau_tissue=tau_values[j])
            snr[i,j] = np.abs(np.mean(progression_sub[3])/np.std(progression_sub[3]))
            #print(np.std(signal[3])) -- for testing
    
    # use log tau values to make a simple scale
    log_tau_values = np.log10(tau_values)

    

    # Generate the heatmap
    start, stop = log_tau_values[0], log_tau_values[-1]
    plt.imshow(snr, cmap='magma', extent=[start, stop, stop, start])
    plt.colorbar(label='snr')
    # plt.xlabel('bleach time constant in tissue fluorescence')
    # plt.ylabel('bleach time constant in dye + nm flourescence')
    plt.title('variance = {}'.format(var))
    #plt.show()

    return snr

# ANALYSIS 1: bleaching factor acting on the [NM] contribution to F
def bleach_nm(K_D, tau, F_max, F_min, nm_conc, bline_len=5000):

    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # bleaching factor -- starts off as 1 then exponentially decreases 
    # we set tau to be a very large constant so this is a slow decrease
    bleach = np.exp(-t/tau)

    # calculate bleached f values: derived in part from eq 2 in Neher/Augustine
    f = bg_tissue + (K_D*F_min + bleach*nm_conc*F_max)/(K_D + nm_conc)

    # define the signal array
    delta_ft_f0 = np.zeros(f.size)

    # calculate f0 values and populate the signal array
    for i in range(f.size):
        
        # calculate f0 by getting the median value of the bottom 70% of previous f values
        
        if i == 0:
            f0 = f[0]
        elif i < bline_len:
            percentile_mark = np.percentile(f[:i],70)
            f0 = np.median(f[f<percentile_mark])
        else: 
            prev_values = f[i-bline_len:i]
            percentile_mark = np.percentile(prev_values,70)
            f0 = np.median(prev_values[prev_values<percentile_mark])


        # calculate normalized signal using the calculated f0
        delta_ft_f0[i] = (f[i]-f0)/(f0)


    # plot the normalized signal delta f/f0 at the different t
    # plt.plot(t,delta_ft_f0, label = tau + 1)
    # plt.xlabel('time(ms)')
    # plt.ylabel('Delta F/F0')
    # plt.title('Flourescence intensity signal over time (bleach 1)')
    # plt.legend()

    return delta_ft_f0

# TESTING IF THE FUNCTION WORKS
# # check this bleaching effect 1 for different values of tau -- different bleach factor
# different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])

# # generate an input nm conc array from firing neuron
# firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
# nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# # plot bleached signal with different time constants for the bleach factor
# for i in range(len(different_taus)):
#    bleach_1(K_D = 1000, tau=different_taus[i], F_max = 45, F_min = 10, nm_conc=nm_conc, bline_len=5000)
# plt.show()



# ANALYSIS 2: bleaching factor acting on the contributions of the dye, F0 and [NM]
def bleach_dnm(K_D, tau, F_max, F_min, nm_conc, bline_len =5000):

    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # bleaching factor -- starts off as 1 then exponentially decreases 
    # we set tau to be a very large constant so this is a slow decrease
    bleach = np.exp(-t/tau)

    # calculate bleached f values: derived in part from eq 2 in Neher/Augustine
    f = bg_tissue+ bleach*(K_D*F_min + nm_conc*F_max)/(K_D + nm_conc)
    
    # fit a polynmial to f and subtract it from f
    poly = np.polyfit(t,f,2)
    fit = np.polyval(poly,t)
    f_sub = f-fit


    # define the df and the df/f signal array
    delta_f = np.zeros(f.size)
    delta_ft_f0 = np.zeros(f.size)

    # calculate f0 values and populate the signal array
    for i in range(f.size):

        # calculate f0 by getting the median value of the bottom 70% of previous f values
        
        if i == 0:
            f0 = f[0]
        elif i < bline_len:
            percentile_mark = np.percentile(f[:i],70)
            f0 = np.median(f[f<percentile_mark])
        else: 
            prev_values = f[i-bline_len:i]
            percentile_mark = np.percentile(prev_values,70)
            f0 = np.median(prev_values[prev_values<percentile_mark])

        # calculate df, and df/f, the normalized signal, using the calculated f0
        delta_f[i] = f_sub[i]-f0
        delta_ft_f0[i] = delta_f[i]/(f0)

    

    # plot f, f - fit, df then df/f
    # plt.figure()
    # plt.subplot(2,2,1)
    # plt.plot(t,f, label='f')
    # plt.xlabel('time (ms)')
    # plt.ylabel('f')
    # plt.legend()

    # plt.subplot(2,2,2)
    # plt.plot(f_sub, label='f subtracted')
    # plt.plot(fit, label='fit')
    # plt.xlabel('time (ms)')
    # plt.ylabel('f-exp')

    # plt.legend()

    # plt.subplot(2,2,3)
    # plt.plot(t,delta_f, label='df')
    # plt.xlabel('time (ms)')
    # plt.ylabel('delta f')
    # plt.legend()


    # plt.subplot(2,2,4)
    # plt.plot(t,delta_ft_f0, label = 'df/f')
    # plt.xlabel('time(ms)')
    # plt.ylabel('Delta F/F0')
    # plt.tight_layout()
    # plt.legend()


    return delta_ft_f0


# TESTING IF THE FUNCTION WORKS
# # check this bleaching effect 1 for different values of tau -- different bleach factor
# different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])

# # generate an input nm conc array from firing neuron
# firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
# nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# # plot bleached signal with different time constants for the bleach factor
# for i in range(len(different_taus)):
#    bleach_2(K_D = 1000, tau=different_taus[i], F_max = 45, F_min = 10, nm_conc=nm_conc)
# plt.show()



# ANALYSIS 3: bleaching factor acting on the background fluorescence 
def bleach_t(K_D, tau, F_max, F_min, nm_conc, bline_len=5000):

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
        # where x is the length of the moving baseline (first element if f[0])

        if i == 0:
            f0 = f[0]
        elif i < bline_len:
            f0 = np.average(f[:i])
        else: 
            f0 = np.average(f[i-bline_len:i])

        # calculate normalized signal using the calculated f0
        delta_ft_f0[i] = (f[i]-f0)/(f0)

    # plot the normalized signal delta f/f0 at the different t
    # plt.plot(t,delta_ft_f0, label = tau + 3)
    # plt.xlabel('time(ms)')
    # plt.ylabel('Delta F/F0')
    # plt.title('Flourescence intensity signal over time (bleach 3)')
    # plt.legend()

    return delta_ft_f0

# TESTING IF THE FUNCTION WORKS
# # check this bleaching effect 1 for different values of tau -- different bleach factor
# different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])

# # generate an input nm conc array from firing neuron
# firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
# nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# # plot bleached signal with different time constants for the bleach factor
# for i in range(len(different_taus)):
#    bleach_3(K_D = 1000, tau=different_taus[i], F_max = 45, F_min = 10, nm_conc=nm_conc)
# plt.show()




# ANALYSIS 4: bleaching on all contributions
def bleach_all(K_D, tau, F_max, F_min, nm_conc, bline_len=5000):

    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # bleaching factor -- starts off as 1 then exponentially decreases 
    # we set tau to be a very large constant so this is a slow decrease
    bleach = np.exp(-t/tau)

    # calculate bleached f values: derived in part from eq 2 in Neher/Augustine
    f= bleach*(bg_tissue + (K_D*F_min + nm_conc*F_max)/(K_D + nm_conc))
    
   # define the signal array
    delta_ft_f0 = np.zeros(f.size)

    # calculate f0 values and populate the signal array
    for i in range(f.size):

        # calculate f0 by averaging the previous x number of f values
        # if x is bigger than the current index then use all the prev f values
        # where x is the length of the moving baseline (first element if f[0])

        if i == 0:
            f0 = f[0]
        elif i < bline_len:
            f0 = np.average(f[:i])
        else: 
            f0 = np.average(f[i-bline_len:i])

        # calculate normalized signal using the calculated f0
        delta_ft_f0[i] = (f[i]-f0)/(f0)

    # plot the normalized signal delta f/f0 at the different t
    # plt.plot(t,delta_ft_f0, label = tau + 4)
    # plt.xlabel('time(ms)')
    # plt.ylabel('Delta F/F0')
    # plt.title('Flourescence intensity signal over time (bleach 4)')
    # plt.legend()

    return delta_ft_f0

# TESTING IF THE FUNCTION WORKS
# # check this bleaching effect 1 for different values of tau -- different bleach factor
# different_taus = np.array([10e6,10e5,10e4,10e3,10e2,10e1,10e0])

# # generate an input nm conc array from firing neuron
# firing_neuron = simulate_neuron(n_timesteps=70000,firing_rate=1)
# nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# # plot bleached signal with different time constants for the bleach factor
# for i in range(len(different_taus)):
#    bleach_4(K_D = 1000, tau=different_taus[i], F_max = 45, F_min = 10, nm_conc=nm_conc)
# plt.show()

