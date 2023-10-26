# This script is used for analysis of the effect of 
# changing the firing rate of the neuron on the measured signal

# import necessary pacakges
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

from s_functions import simulate_neuron, simulate_nm_conc, simulate_fluorescence_signal


# ANALYSIS 1:

# Define function that makes generates several firing neurons with different firing rates
# and then compares the average df/f signal across them

def signal_vs_activity(firing_rates, bleach_time):

    # create array to store the average signals generated by the neurons
    average_signals = []

    # for each firing rate, simulate a neuron, the nm dynamics from it and the signal
    # then get the average signal
    for i in range(firing_rates.size):
        
        # simulate neuron with different firing rate
        guinea_neuron = simulate_neuron(70000,firing_rates[i])

        # generate the nm_conc
        guinea_nm_conc, guinea_b_conc, guinea_c_conc, guinea_nm_tot = simulate_nm_conc(guinea_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

        # then generate the signal
        guinea_signal = simulate_fluorescence_signal(tau_d=bleach_time, tau_nm=bleach_time, tau_tissue=10e7, nm_conc=guinea_nm_conc)

        # get the average of this signal and add it the original array
        average_signals.append(np.average(guinea_signal))



    return average_signals

def plot_different_bleach(firing_rates,bleach_times):

    " this function takes in a list of bleach time constants and runs the previous function to "
    " obtain average signal plots at the different bleach strengths "

    # create array of the average signal plots
    average_signal_plots = []

    # generate the different plots for the bleach time constants
    for i in range(len(bleach_times)):

        # generate the average signal plot at the specific bleach time constant
        average_signal_plot = signal_vs_activity(firing_rates,bleach_times[i])

        # store the values in the array
        average_signal_plots.append(average_signal_plot)


    # make the average signal plots at different bleach time constants
    plt.figure()
    for i in range(len(bleach_times)):
        plt.plot(firing_rates, average_signal_plots[i], label=np.round(bleach_times[i],0))
   
    plt.xlabel('Firing rates(Hz)')
    plt.ylabel('Average df/f signal')
    plt.title('Signal vs activity plot')
    plt.legend()
    plt.show()


    # plot at log scale
    plt.figure()
    for i in range(len(bleach_times)):
        plt.loglog(firing_rates, average_signal_plots[i], label=np.round(bleach_times[i],0))

    plt.xlabel('Firing rates(Hz)')
    plt.ylabel('Average df/f signal')
    plt.title('Signal vs activity plot')
    plt.legend()
    plt.show()


different_firing_rates = np.linspace(1,50,15)
# Check 1: -- It works
# signal_vs_activity(different_firing_rates,10e7)


# Check 2:
list_of_bleaches = np.logspace(-4,13,5)
plot_different_bleach(different_firing_rates,bleach_times=list_of_bleaches)






# TO ADD:
# then also get an r2 estimate

