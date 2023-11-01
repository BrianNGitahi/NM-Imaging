import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit

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
    # print('Simulated neuron with {} spikes in {} timesteps ({} Hz).'.format(n_spikes, n_timesteps, firing_rate))
  

    return firing_neuron



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
        
        # if there's a spike add Delta_nm else subtract d_nm/dt
        if neuron_activity[t]==1: 
            nm_conc[t] = nm_conc[t] + delta_nm
        else: 
            d_nm_dt =  (nm_conc[t]-nm_conc0)/tau
            nm_conc[t] = nm_conc[t] - d_nm_dt 

    # plot the [NM] at all timesteps
    n_timesteps = neuron_activity.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    # Calculate the concentrations of the bound forms of the NM
    # start with [NM B] the NM bound to the sensor
    nm_b_conc = k_b*nm_conc

    # then [NM R], the NM bound to the receptor
    nm_r_conc = k_r*nm_conc

    # then get the total nm concentration - both bound and unbound
    nm_tot = nm_conc + nm_b_conc + nm_r_conc

  
    # return the array of the [NM], [NM B], and [NM R]
    return nm_conc, nm_b_conc, nm_r_conc, nm_tot

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




def simulate_fluorescence_signal(tau_d, tau_nm, tau_tissue, nm_conc, K_D = 1000, F_max = 45, F_min = 10, bline_len=5000):

    # autofluorescence
    f_tissue = 0.02

    # create timesteps 
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps) 

    # define bleach factors for the autofluorescence and fluorescence from dye + nm
    bleach_d = np.exp(-t/tau_d)
    bleach_nm = np.exp(-t/tau_nm)
    bleach_tissue = np.exp(-t/tau_tissue)
    
    # calculate F: derived from eq 2 in Neher/Augsutine
    f = bleach_tissue*f_tissue + (bleach_d*K_D*F_min + bleach_nm*nm_conc*F_max)/(K_D + nm_conc)


    # fit an exponential to remove the bleachign trend 

    # define an exponential function that we'll use as the basis for the fit
    def exp_decay(t,a,b):
        return a*np.exp(-t/b)

    # perform the fit
    params, covariance = curve_fit(exp_decay,t,f)

    # get the parameters
    a_fit, b_fit = params

    # define the fitted function
    fit = exp_decay(t,a_fit,b_fit)
    
    # subtracted f
    f_subtracted = f - fit

    # to correct for negative values in the fluorescence that result in a -ve df/f
    f_alt = f_subtracted + np.max(np.abs(f_subtracted))

    
    # calculate f0 by getting the median value of the bottom 70% of previous f values
    percentile_mark = np.percentile(f,70)
    f0 = np.median(f[f<percentile_mark])
    
    # df calc -- median value method
    df = f-f0
    df_f_med = df/f0


    # df/f with the subtracted formula
    # subtracted signal
    percentile_mark_prime = np.percentile(f_alt,70)
    f0_sub = np.median(f_alt[f_alt<percentile_mark_prime])
    df_sub = f_alt - f0_sub
    
    df_f_med_sub = df_sub/f0_sub
    


    # define the delta f and the df/f signal arrays
    df_f_ave = np.zeros(f.size)
    f0_averages = np.zeros(f.size)
    

    # calculate f0 values and populate the signal array
    for i in range(f.size):

        # calculate f0 using the average method
        if i==0:
            f0_ave=f[0]
        
        elif i < bline_len:
            f0_ave = np.average(f[:i]) 
        else: 
            f0_ave = np.average(f[i-bline_len:i]) 

        # calculate normalized signal using the calculated f0
    
        
        # average value
        df_f_ave[i] = (f[i] - f0_ave)/f0_ave
        f0_averages[i]=f0_ave

    # define progression arrays for f, df, df/f
    progression = []
    progression_sub = []

    # without the subtracted exponent
    progression.append(f)
    progression.append(df)
    progression.append(df_f_med)
    progression.append(df_f_ave)

    # with the subtracted exponent
    progression_sub.append(f)
    progression_sub.append(f_alt)
    progression_sub.append(df_sub)
    progression_sub.append(df_f_med_sub)


    return progression, progression_sub
# output 3
#f_signal = simulate_flourescence_signal(K_D = 1000, F_max = 45, F_min = 10, nm_conc=nm_conc)



    





