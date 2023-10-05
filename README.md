# NM-Imaging

This repository hosts the notebooks and scripts that we'll be using to simulate the NM  imaging pipeline: starting from neuron activity to the generated df/f signal.
The main goal is to assess how the calculation of f can affect the df/f signal produced. 

- Simulation.ipynb is the notebook with the documentation of this NM imaging model
- Simulation.py will the script that calls the functions in Functions.py

This repository is still new, so there's probably still errors in the code -- if you find some please raise an issue and we can resolve it.

**Specifics on the simulation**

We're modelling (simulating) the dynamics of neuromodulator concentration [x] and neuromodulator bound to sensor [x s] and to receptor [x r]. 
This is equivalent to intracellular calcium measurements in the presence of endogenous ca buffers and exogenous buffers (the sensors). 
Here, we assume that we are dealing with one compartment in equilibrium to start and everything is linear, including uptake. 
The Source of x is a train of spikes, each producing a fixed delta_x. -- the inspiration for the model comes from Neher's work with calcium that Karel suggested on email.
