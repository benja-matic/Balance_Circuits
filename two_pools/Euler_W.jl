function von_mises_dist(x, k, mu, N)
  a = exp(k*cos(x-mu))/(N*besseli(0, k))
end

#

function Weights(Ne,Ni,kee,kei,kie_L,kii,Aee,Aei,Aie_L,Aii,p)
# Construct weight functions
N = Ne+Ni

we = 2*collect(0:Ne-1)*pi/Ne
wi = 2*collect(0:Ni-1)*pi/Ni

vee = Aee*von_mises_dist(we, kee, 0, Ne)
vei = -Aei*von_mises_dist(wi, kei, 0, Ni)
vie = Aie_L*von_mises_dist(we, kie, 0, Ne)
vii = -Aii*von_mises_dist(wi, kii, 0, Ni)

W = zeros(N,N)
for i = 1:Ne
    W[i,1:Ne]=circshift(vee, i-1)
    W[i,Ne+1:end]=circshift(vei,div(Ni*i,Ne)-1) #
end
for i = 1:Ni
    W[Ne+i,1:Ne]=circshift(vie,div(Ne*i,Ni)-1)
    W[Ne+i,Ne+1:end]=circshift(vii,i-1)
end

Wx = rand(N,N) .<= p
W .*= Wx

return W
end

function sparse_rep(W,N)
    flat = []
    for i = 1:N
        mi = find(W[:,i])
        push!(flat, mi)
    end
    return flat
end

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  return t
end

function euler_lif_CSR(h, total, Ne, W, CSR, s1, s2, vth, tau_m, tau_s, tau_a, g_a, Angle)

  ntotal = round(Int64, total/h)
  syn = zeros(N) #synapse
  # syn_e = zeros(N) #for measuring just excitatory synapses
  # syn_i = zeros(N) #for measuring just inhibitory synapses
  A = zeros(N) #adaptation
  V = rand(N)*vth #voltage
  V_buff = V
  # Input = zeros(N, ntotal)
  # adapt = zeros(N, ntotal)
  # E_input = zeros(N, ntotal) #parsed excitatory inputs coming into a neuron
  # I_input = zeros(N, ntotal) #parsed inhibitory inputs coming into a neuron
  # V_store = zeros(N, ntotal) #store the voltage

  drive = zeros(N) #feedforward input
  FPe = div(Ne,5)
  P1s = round(Int,Ne/4-FPe) + Angle
  P1e = round(Int,Ne/4+FPe) + Angle
  P2s = round(Int,3*Ne/4-FPe) - Angle
  P2e = round(Int,3*Ne/4+FPe) - Angle

  drive[P1s:P1e] = s1*h
  drive[P2s:P2e] = s2*h

  time = Float64[0]
  raster = Float64[0]

  m_leak = h/tau_m
  s_leak = h/tau_s
  a_leak = h/tau_a

  kill_flag = false
  for iter = 1:ntotal

      incoming = drive .+ (h .* syn) .- (A .* g_a)
      # Input[:, iter] = incoming
      V .+= incoming .- (V .* m_leak)

      # E_input[:,iter] = drive .+ (h .* syn_e)
      # I_input[:,iter] = (h .* syn_i) .- (A .* g_a)

      syn .-= (syn .* s_leak)
      # syn_e .-= (syn_e .* s_leak)
      # syn_i .-= (syn_i .* s_leak)
      A .-= A*a_leak

      # adapt[:, iter] = A
      # V_store[:,iter] = V

      vs = (V .> vth)
      vsm = sum(vs)

      if vsm > 0
        spe = find(vs)
        for j = 1:vsm
          js = spe[j]
          delta_h = interpolate_spike(V[js], V_buff[js], vth)
          lx = exp(delta_h/tau_m)
          syn[CSR[js]] += W[CSR[js], js] .* lx
          push!(raster, js)
          push!(time, iter-delta_h)
          if spe[j] <= Ne
            A[spe[j]] += 1.
            # syn_e[CSR[js]] .+= W[CSR[js], js] .* lx
          # else
            # syn_i[CSR[js]] .+= W[CSR[js], js] .* lx
          end
        end
      end

      V -= vth*vs
      V_buff = V

    end

    # return time, raster, Input, adapt, E_input, I_input, V_store
    return time, raster
  end
