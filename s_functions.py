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


    # check exactly how many spikes were produced: to see if it works
    n_spikes = np.size(np.nonzero(firing_neuron))

    # print simulated neuron summary:
    #print('Simulated neuron with {} spikes in {} timesteps ({} Hz).'.format(n_spikes, n_timesteps, firing_rate))
  
    return firing_neuron

# test 1
# plt.figure(1)
# plt.subplot(2,1,1)
# firing_neuron = simulate_neuron(70000,1)
# plt.subplot(2,1,2)
# firing_neuron2 = simulate_neuron(70000,10)
# plt.tight_layout()
# plt.show()

# function 1.1: regularly firing neuron
def reg_neuron(n_timesteps, firing_rate, number=1):

    # get the numner of seconds: n_timesteps is in ms
    seconds = n_timesteps/1000

    # get the number of spikes: firing rate is in Hz
    n_spikes = seconds*firing_rate

    # generate indices where spiking occurs
    firing_indices = np.linspace(0,n_timesteps-1,n_spikes)

    # make an array to reflect this activity
    firing_neuron = np.zeros(n_timesteps)
    firing_neuron[firing_indices]=1

    print('number of spikes {}.'.format(np.size(np.nonzero(firing_neuron))))


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

    # plot [NM], [NM B] and [NM R] simulataneously
    # plt.plot(t,nm_conc, color = 'b', label='[NM]')
    # plt.plot(t,nm_b_conc, color = 'g', label='[NM B]')
    # plt.plot(t,nm_r_conc, color = 'r', label='[NM R]')


    # # label the axes and make legend
    # plt.xlabel('time (ms)')

    # # # to zoom in on a plot
    # # plt.xlim(5000,15000)

    # plt.ylabel('(Change in) Concentration -- arbitrary units')
    # plt.title('NM concentration across {} ms'.format(n_timesteps))
    # plt.legend()
    
  
    # return the array of the [NM], [NM B], and [NM R]
    return nm_conc, nm_b_conc, nm_r_conc

# test 2
# plt.figure(2)
# plt.subplot(2,1,1)
# nm_conc, nm_b_conc, nm_r_conc = simulate_nm_conc(firing_neuron,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)
# plt.subplot(2,1,2)
# nm_conc2, nm_b_conc2, nm_r_conc2 = simulate_nm_conc(firing_neuron2,nm_conc0=0,k_b=0.6, k_r=0.4,gamma=0.004)
# plt.tight_layout()
# plt.show()

# Function 2.1: produces zoomed in plot
def plot_nm_conc(nm,nmb, nmr, colour='b'):

    # define the timesteps to plot the [NM]
    t = np.linspace(0,nm.size-1,nm.size)

    # plot the [NM] 
    plt.plot(t,nm, color=colour, label='[NM]')
    plt.plot(t,nmb, color = 'g', label='[NM B]')
    plt.plot(t,nmr, color = 'r', label='[NM R]')
    plt.xlabel('time (ms)')
    plt.ylabel('NM concentration')
    plt.title('NM dynamics from firing neuron')
    plt.legend()
    plt.show()

    



def simulate_fluorescence_signal(tau_d, tau_nm, tau_tissue, nm_conc, variance=0.0005, K_D = 1000, F_max = 45, F_min = 10, bline_len=5000):

     # noisy autofluorescence
    f_tissue = np.random.normal(0.002,variance,nm_conc.size)

    # create timesteps 
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps) 

    # define bleach factors for the autofluorescence and fluorescence from dye + nm
    bleach_d = np.exp(-t/tau_d)
    bleach_nm = np.exp(-t/tau_nm)
    bleach_tissue = np.exp(-t/tau_tissue)
    
    # calculate F: derived from eq 2 in Neher/Augustine
    f = bleach_tissue*f_tissue + (bleach_d*K_D*F_min + bleach_nm*nm_conc*F_max)/(K_D + nm_conc)


    # fit an exponential to remove the bleaching trend 

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
    #progression_multi = []

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



# plot the progression from f -- df/f
def plot_f_signal(progression, progression_sub, nm_conc):


    # create timesteps array for the plot
    n_timesteps = nm_conc.size
    t = np.linspace(0,n_timesteps-1,n_timesteps)

    plt.figure(1)
    plt.subplot(2,2,1)
    plt.plot(t,progression[0], label='f')
    plt.xlabel('time (ms)')
    plt.ylabel('f')
    plt.title('f  vs time')
    plt.legend()
    
    plt.subplot(2,2,2)
    plt.plot(t,progression[1], label='df')
    plt.xlabel('time (ms)')
    plt.ylabel('df')
    plt.title('df vs time (f0:median)')
    plt.legend()

    plt.subplot(2,2,3)
    plt.plot(t,progression[2], label = 'df/f')
    plt.xlabel('time(ms)')
    plt.ylabel(' df/f')
    plt.title('df/f vs time (f0:median)')
    plt.legend()
   
    plt.subplot(2,2,4)
    plt.plot(t,progression[3], label = 'df/f')
    plt.xlabel('time(ms)')
    plt.ylabel(' df/f')
    plt.title('df/f vs time (f0:average)')
    plt.legend()
    plt.suptitle('Progression from f to df/f', size = 16)
    plt.tight_layout()

    plt.figure(2)
    plt.subplot(2,2,1)
    plt.plot(t,progression_sub[0], label='f')
    plt.xlabel('time (ms)')
    plt.ylabel('f')
    plt.title('f  vs time')
    plt.legend()
    
    plt.subplot(2,2,2)
    plt.plot(t,progression_sub[2], label = 'df')
    plt.xlabel('time(ms)')
    plt.ylabel(' df')
    plt.title('df vs time (f0:median)')
    plt.legend()
    
    plt.subplot(2,2,3)
    plt.plot(t,progression_sub[3], label = 'df/f')
    plt.xlabel('time(ms)')
    plt.ylabel(' df/f')
    plt.title('df/f vs time (f0:median)')
    plt.legend()

    plt.subplot(2,2,4)
    plt.plot(t,progression_sub[1], label='f - exponential (+ constant)')
    plt.xlabel('time (ms)')
    plt.ylabel('f_sub')
    plt.title('f_sub vs time')
    plt.legend()
    
    plt.suptitle('Progression from f to df/f (using subtracted f)', size = 16)
    plt.tight_layout()
    plt.show()



# # output 3
# progression, progression_sub  = simulate_fluorescence_signal(tau_d=10e9,tau_nm=10e9, tau_tissue=10e9, nm_conc = nm_conc)
# progression2, progression_sub2 = simulate_fluorescence_signal(tau_d=10e9,tau_nm=10e9, tau_tissue=10e9, nm_conc = nm_conc2)

# plot_f_signal(progression, progression_sub)
# print('average signal is {}'.format(np.average(progression[2])))
# plt.show()
# plot_f_signal(progression2, progression_sub2)
# print('average signal is {}'.format(np.average(progression2[2])))
# plt.show()





    





