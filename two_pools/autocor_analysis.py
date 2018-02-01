###get the autocovariance matrix of a vector across all time lags tau
import numpy as np
import matplotlib.pyplot as plt

def single_neuron(runtime, drive, mu, sigma, tau):
    vth = 20.
    tau_m = 20.
    h = 0.1
    ht = h/tau
    htm = h/tau_m
    ntotal = int(runtime/h)

    #get weiner process first
    r = np.random.normal(0, 1, ntotal)
    sh = np.sqrt(h)/tau
    r *= sh
    rn = np.zeros(ntotal)
    x0 = np.random.normal()
    x = x0
    for i in range(ntotal):
        dxh = -x*ht + r[i]
        rn[i] = x + dxh
        x += dxh

    #rescale process
    scale = sigma*tau
    rn *= scale
    v0 = np.random.uniform(0, 1) * vth

    inputs = np.zeros(ntotal)
    spikes = []
    v = v0
    d = drive*h
    for i in range(ntotal):
        inputs[i] = mu*h + rn[i]
        I = inputs[i] + d -v*htm
        v += I
        if v >= vth:
            spikes.append(i)
            v -= vth
    return spikes, inputs

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

def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]

def autocovariance(x):
    n = len(x)
    autocov = np.zeros(n-1)
    for i in range(n-1):
        autocov[i] = np.inner(x[i:], x[:n-i])/(n-i-1)
    return autocov





los = parse_inputs("e_bot_l1.txt")
win = parse_inputs("e_top_w1.txt")
nrm = parse_inputs("e_top_norm1.txt")

los_m

la = autocorr(los)
wa = autocorr(win)
na = autocorr(nrm)

los_m = los - np.mean(los)
win_m = win - np.mean(win)
nrm_m = nrm - np.mean(nrm)

la_m = autocovariance(los_m)
wa_m = autocovariance(win_m)
na_m = autocovariance(nrm_m)
time = np.linspace(0, 14000, len(la_m))
plt.plot(time, la_m, "b")
plt.plot(time, wa_m, "r")
plt.plot(time, na_m, "g")
plt.title("Autocovariance: r = win, b = lose, g = normalization")
plt.xlabel("Time-Lag (ms)")
plt.ylabel("Inner Product")
plt.show()

# spikes, inputs = single_neuron(14000, 3.08, m_norm, s_norm, 1.)
