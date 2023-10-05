import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

# This script contains the functions in the Simulation notebook
# Its purpose is to run all the functions and produce a df/f plot

# TO DO: make separate functions to plot the results of each of these functions
# print output line after function 2 to summarize the result: average conc for each of them maybe 
# make the function 1 able to make multiple neurons at once!
# add the firing rate to the plots of flourescence

# check the inputs to each of the functions and raise errors if something's wrong

# Bonus: find way to show the figures all at once


# Function 1: firing neuron
def simulate_neuron(n_timesteps, firing_rate, number=1):
    
    # generate a base array with random numbers between 0-1
    x = np.random.rand(n_timesteps)

    # Then populate the bins with signals - firing & not firing -- with a specific probability that you choose
    firing_neuron = x < 0.001*firing_rate
    firing_neuron = firing_neuron.astype(int)


    # # then make a plot of it!
    # plt.plot(firing_neuron)
    # plt.xlabel('timesteps')
    # plt.ylabel('Neuron activity')
    # plt.title('Neuron Activity over {} timesteps'.format(n_timesteps))
    # plt.show()


    # check exactly how many spikes were produced: to see if it works
    n_spikes = np.size(np.nonzero(firing_neuron))

    # print simulated neuron summary:
    print('Simulated neuron with {} spikes in {} timesteps ({} Hz).'.format(n_spikes, n_timesteps, firing_rate))
  

    return firing_neuron

# output 1
#firing_neuron = simulate_neuron(70000,1)


# Function 2: takes in an array of neuron activity and gives corresponding [NM]
def simulate_nm_conc(neuron_activity,nm_conc0, k_b,k_r,gamma):


    # create array of [NM] with same size as neuron activity
    nm_conc = np.zeros(neuron_activity.size)

    # define delta_nm, the increase in [NM] due to a spike 
    # this will be a constant value  -- calculate the amplitude of the exponential function, A
    delta_nm = 1
    
    # first define tau the time constant
    tau = (1+k_r+k_b)/gamma

    # create a for-loop where we update the value of [NM] for current timestep
    for t in range(neuron_activity.size):
        
        # first timebin condition
        if t == 0 : 
            nm_conc[t] = nm_conc0 
        else: 
            nm_conc[t] = nm_conc[t-1]

        # update [NM] value

        # define d_nm/dt, the inifinitesimal decay in [NM] in one timestep 
        d_nm_dt =  (nm_conc[t]-nm_conc0)/tau

        # if there's a spike add Delta_nm else subtract d_nm/dt
        if neuron_activity[t]==1: nm_conc[t] = nm_conc[t] + delta_nm
        else: nm_conc[t] = nm_conc[t] - d_nm_dt 

    # plot the [NM] at all timesteps
    n_timesteps = neuron_activity.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # Calculate the concentrations of the bound forms of the NM
    # start with [NM B] the NM bound to the sensor
    nm_b_conc = k_b*nm_conc

    # then [NM R], the NM bound to the receptor
    nm_r_conc = k_r*nm_conc

    # # plot [NM], [NM B] and [NM R] simulataneously
    # plt.plot(t,nm_conc, color = 'b', label='[NM]')
    # plt.plot(t,nm_b_conc, color = 'g', label='[NM B]')
    # plt.plot(t,nm_r_conc, color = 'r', label='[NM R]')

    # # label the axes and make legend
    # plt.xlabel('time (ms)')
    # plt.ylabel('Concentration')
    # plt.title('NM concentration across {} ms'.format(n_timesteps))
    # plt.legend()
    # plt.show() 
   

    # return the array of the [NM], [NM B], and [NM R]
    return nm_conc, nm_b_conc, nm_r_conc

# output 2
#nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

# Function 2.1: produces zoomed in plot
def plot_nm_conc(nm,start,stop,colour='b', plotlabel = ''):

    # define the timesteps to plot the [NM]
    timesteps = stop - start + 1
    t = np.linspace(start,stop,timesteps)

    # get that section of the [NM] array
    nm_section = nm[start:stop+1]

    # plot the [NM] 
    plt.plot(t,nm_section, color=colour, label=plotlabel)
    plt.xlabel('time (ms)')
    plt.ylabel('NM concentration')
    plt.title('NM {} concentration from {} to {} ms'.format(plotlabel, start,stop))
    plt.show()

# output 2.1
# plot_nm_conc(nm_conc, start = 22000,stop = 26000)
# plot_nm_conc(nm_b_conc, start = 22000,stop = 26000, colour='g', plotlabel='B')
# plot_nm_conc(nm_r_conc, start = 22000,stop = 26000, colour='r', plotlabel='R')


# function 3: get that signal!
def simulate_flourescence_signal(K_D, F_max, F_min, nm_conc):
    
    # define K_D prime as
    K_Dp = K_D*(F_min/F_max)

    # the initial/steady state concentration, [NM]i,0, of the neuromdultor
    # CONFIRM VALUE FROM KENTA
    nm_conc_0 = 0 

    # define the numerator and denominator
    numerator = (K_Dp + nm_conc)/(K_D + nm_conc)
    denominator = (K_Dp + nm_conc_0)/(K_D + nm_conc_0)

    # derive delta f/f0 by plugging in
    delta_ft_f0 = (numerator/denominator) - 1

    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # plot the normalized signal delta f/f0 at the different t
    # plt.plot(t,delta_ft_f0)
    # plt.xlabel('time(ms)')
    # plt.ylabel('Delta F/F0')
    # plt.title('Flourescence intensity signal over time')
    # plt.show()

    print('completed f_signal simulation')

    return delta_ft_f0


# output 3
#f_signal = simulate_flourescence_signal(K_D = 1000, F_max = 45, F_min = 10, nm_conc=nm_conc)


# ANALYSIS 1:

# Define function that makes generates several firing neurons with different firing rates
# and then compares the average df/f signal across them

def signal_vs_activity(number_neurons,firing_rates):

    # create array to store the average signals generated by the neurons
    average_signals = []

    # for each firing rate, simulate a neuron, the nm dynamics from it and the signal
    # then get the average signal
    for i in range(firing_rates.size):
        
        # simulate neuron with different firing rate
        guinea_neuron = simulate_neuron(70000,firing_rates[i])

        # generate the nm_conc
        guinea_nm_conc, guinea_b_conc, guinea_c_conc = simulate_nm_conc(guinea_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)

        # then generate the signal
        guinea_signal = simulate_flourescence_signal(K_D = 1000, F_max = 45, F_min = 10, nm_conc=guinea_nm_conc)

        # get the average of this signal and add it the original array
        average_signals.append(np.average(guinea_signal))

        print('simulated {} neuron(s)'.format(i+1))


    # plot the different firing rates vs their signals
    plt.scatter(firing_rates,average_signals)
    plt.xlabel('Firing rates(Hz)')
    plt.ylabel('Average df/f signal')
    plt.title('Signal vs activity plot')
    plt.show()

    # plot at log scale
    plt.loglog(firing_rates,average_signals, 'o')
    plt.xlabel('Firing rates(Hz)')
    plt.ylabel('Average df/f signal')
    plt.title('Signal vs activity plot -Log')
    plt.show()
    



# Check 1:
signal_vs_activity(2,np.array([1,2,4,10,50,70,100,200,300,400,500,600,700,800,900,1000,1500,2000,2500,3000,5000]))

# results: looks plausible -- tho at a big scale it has a linear/exponential trend then it levels off 




