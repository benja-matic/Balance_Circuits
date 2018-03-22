# Balanced State Networks
Here are some codes used for simulating networks of integrate and fire neurons. I used these codes to study the balanced state, winner-take-all, and rivalry.

## Local Circuits
Local circuits include networks with random or structured connections, delta pulse synapses or synapses with a time constant, and networks that respect, or don't respect, Dale's law.

## Non-Local Circuits
Non-local circuits include two E/I pools that are coupled.
Simplest*.jl contain networks with explicit long-range connections from one excitatory pool to the opposing inhibitory pool.

Euler_W.jl contains all the machinery to actually simulate a network. 
Run_Euler_W.jl is for the user. Use this to play with parameters, etc. 
Analyze.jl contains all the machinery to analyze simulations.

Zero out the parameter g_a to simulate with no spike-frequency-adaptation. If you're in winner-take-all, add g_a back in to get rivalry.

The python code was used to model the distribution of inputs coming into a neuron.

## Various Functions
Junkyard of functions I may have used at one point.
