# This script is going to be used to compare the df/f vs firing rate plot for different bleaches

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

from s_functions import simulate_neuron, simulate_nm_conc, simulate_fluorescence_signal
from b_factory import bleach_dnm_heat


# define the base function

def dff_v_activity(neuron_activity, progression_array, window=100):

    number_blocks = int(progression_array[3].size/window)
    signal_points = []


    # array of neural activity
    firing_rate = []

    for block in range(number_blocks):

        # define the interval
        start = (block*window)
        stop = ((block+1)*window+1)

        # get average signal for 100 ms sliding window
        signal_average = np.average(progression_array[3][start:stop])

        # store it
        signal_points.append(signal_average)

        # Check the activity of the firing neuron for current window
        spikes = np.size(np.nonzero(neuron_activity[start:stop]))

        # convert it to Hz and add it to the firing rate array
        firing_rate.append(spikes*1000/100)


    # get the unique firing rates
    unique_frates = np.unique(firing_rate)

    # array to store the average dff values
    average_signal_points = []

    # convert the list to an array
    signal_points = np.array(signal_points)


    # average the df/f values at each of the unique values
    for i in range(len(unique_frates)):

        # get the indices where the firing rate was a given value, x
        indices = np.where(firing_rate==unique_frates[i])[0]
        
        # average all the corresponding df/f values
        ave_dff = np.average(signal_points[indices])
        average_signal_points.append(ave_dff)


    # fit the average signal points 
    poly = np.polyfit(unique_frates,average_signal_points,1)
    fit = np.polyval(poly,unique_frates)

    return unique_frates, average_signal_points, fit


def plot_dff_v_activity(firing_rates, dff, fit):


    plt.plot(firing_rates,dff,'o')
    plt.plot(firing_rates,fit, label='fit')
    plt.xlabel('Firing rates(Hz)')
    plt.ylabel('Average df/f signal')
    plt.title('Signal vs activity plot (average)')
    plt.legend()
    plt.show()


# Test 1

# At different bleaches simulate a neuron firing and get the signal vs firing rate plot
# bleaches = np.logspace(5,7,5)

# # generate a firing neuron
# neuron = simulate_neuron(n_timesteps=70000,firing_rate=13)

# # generate nm_conc 
# nm_conc_input, nm_b_conc, nm_r_conc = simulate_nm_conc(neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# # plot the df/f vs firing rate at diff bleaches
# plt.figure()

# for i in range(len(bleaches)):

#     # genrate the progression array for fluorescence
#     progression, progression_sub = simulate_fluorescence_signal(nm_conc=nm_conc_input,tau_d=bleaches[i],tau_nm=bleaches[i], tau_tissue=10e9)
   
#     # plot the signal vs firing rate
#     x, y, fit = dff_v_activity(neuron,progression_sub)
#     #plt.scatter(x,y,label=np.round(bleaches[i]))
#     plt.plot(x,fit,label=np.round(bleaches[i]))


# plt.xlabel('Firing rates(Hz)')
# plt.ylabel('Average df/f signal')
# plt.title('Signal vs activity plot')
# plt.legend(title='bleach time constants')
# plt.show()




