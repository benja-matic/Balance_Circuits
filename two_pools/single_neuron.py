#single_neuron.py


def single_neuron(runtime, mu, sigma, tau):
    vth = 20.
    tau_m = 20.
    h = 0.1
    ht = h/tau
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
    for i in range(ntotal):
        inputs[i] = -v*h/tau_m + mu*h + rn[i]
        v += inputs[i]
        if v >= vth:
            spikes.append(i)
            v -= vth
    return spikes, inputs
