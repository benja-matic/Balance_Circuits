###get the autocovariance of a vector across time lags tau
import numpy as np
import matplotlib.pyplot as plt

#simulate a single neuron with OU noise
def single_neuron(runtime, drive, mu, sigma, tau):

    vth = 20. #threshold
    tau_m = 20. #membrane time constant
    h = 0.1 #dt
    ht = h/tau #
    htm = h/tau_m
    ntotal = int(runtime/h) #time points

    inputs = np.random.normal(mu*h, sigma*h, ntotal)
    inputs *= np.sqrt(h)
    spikes = []
    v0 = 15.
    v = v0
    d = drive*h
    for i in range(ntotal):
        I = inputs[i] + d -v*htm
        v += I
        if v >= vth:
            spikes.append(i)
            v -= vth
    return spikes, inputs


def FI(runtime, mu):

    vth = 20. #threshold
    tau_m = 20. #membrane time constant
    h = 0.1 #dt
    htm = h/tau_m
    ntotal = int(runtime/h) #time points
    sigma = np.sqrt(mu)
    inputs = np.random.normal(mu, sigma, ntotal)
    inputs *= np.sqrt(h)
    spikes = []
    v0 = 15.
    v = v0
    for i in range(ntotal):
        I = inputs[i]
        v += I - v*htm
        if v >= vth:
            spikes.append(i)
            v -= vth
    return spikes, inputs
#read in a time series from a file
def parse_inputs(f):
    a = open(f, 'r')
    b = a.read()
    a.close()
    c = b.split(',')
    c.pop()
    d = np.zeros(len(c))
    for i in range(len(c)):
        d[i] = float(c[i])
    return d

#Convolves an array against itself and gives the inner product at each increment
def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]


#Convolves an array against itself and gives the inner product at each increment normalized by the length of the non-0 part
def autocovariance(x):
    n = len(x)
    autocov = np.zeros(n-1)
    for i in range(n-1):
        autocov[i] = np.inner(x[i:], x[:n-i])/(n-i-1)
    return autocov

#My implementation of numpy's autocor but convolve from 0-infinity not -infinity to infinity
def autocovariance_n(x):
    n = len(x)
    autocov = np.zeros(n-1)
    for i in range(n-1):
        autocov[i] = np.inner(x[i:], x[:n-i])
    return autocov


spikes, inputs = FI(10000, .2)

plt.plot(spikes, np.ones(len(spikes)), "g.")
plt.show()


# los = parse_inputs("e_bot_l1.txt")
# win = parse_inputs("e_top_w1.txt")
# nrm = parse_inputs("e_top_norm1.txt")
#
# los_m
#
# la = autocorr(los)
# wa = autocorr(win)
# na = autocorr(nrm)
#
# los_m = los - np.mean(los)
# win_m = win - np.mean(win)
# nrm_m = nrm - np.mean(nrm)
#
# la_m = autocovariance(los_m)
# wa_m = autocovariance(win_m)
# na_m = autocovariance(nrm_m)
# time = np.linspace(0, 14000, len(la_m))
# plt.plot(time, la_m, "b")
# plt.plot(time, wa_m, "r")
# plt.plot(time, na_m, "g")
# plt.title("Autocovariance: r = win, b = lose, g = normalization")
# plt.xlabel("Time-Lag (ms)")
# plt.ylabel("Inner Product")
# plt.show()

# spikes, inputs = single_neuron(14000, 3.08, m_norm, s_norm, 1.)
