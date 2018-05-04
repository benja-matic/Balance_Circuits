function homogenous_4x4_weights(N, IFRAC, k, Aee, Aei, Aie, Aie_NL, Aii)

  N_local = Int64(round(N/2)) #divide the network into two EI circuits
  Ni_local = Int64(round(N_local/IFRAC))
  Ne_local = Int64(N_local-Ni_local)
  Ne2 = Ne_local*2
  ks = sqrt(k)

  Jee = Aee/ks
  Jei = -Aei/ks
  Jie = Aie/ks
  Jie_NL = Aie_NL/ks
  Jii = -Aii/ks

  W = zeros(N,N)

  ###EE AND EI FOR BOTH CIRCUITS
  for i = 1:Ne_local

    ee1_inds = rand(1:Ne_local, k) #EE local
    ei1_inds = rand(Ne2+1:Ne2+Ni_local, k) #IE local
    ee2_inds = rand(Ne_local+1:Ne2, k) #other EE local
    ei2_inds = rand(N_local + Ne_local+1:N, k) #other IE local

    for j in eachindex(ee1_inds) #pool1 EE
      W[i, ee1_inds[j]] += Jee
      W[i+Ne_local, ee2_inds[j]] += Jee
      W[i, ei1_inds[j]] += Jei
      W[i+Ne_local, ei2_inds[j]] += Jei
    end

  end

  ###IE AND II FOR BOTH CIRCUITS
  for i = Ne2+1:Ne2+Ni_local

    ie1_inds = rand(1:Ne_local, k)
    ii1_inds = rand(Ne2+1: Ne2+Ni_local, k)
    ie2_inds = rand(Ne_local+1:Ne2, k)
    ii2_inds = rand(Ne2+Ni_local+1:N, k)

    ieNL1_inds = rand(Ne_local+1:Ne2, k)
    ieNL2_inds = rand(1:Ne_local, k)

    for j in eachindex(ie1_inds)
      W[i, ie1_inds[j]] += Jie
      W[i, ii1_inds[j]] += Jii
      W[i+Ni_local, ie2_inds[j]] += Jie
      W[i+Ni_local, ii2_inds[j]] += Jii
      W[i, ieNL1_inds[j]] += Jie_NL
      W[i+Ni_local, ieNL2_inds[j]] += Jie_NL
    end

  end

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

function zero_long_range(W)
  W[N2+1:N2+NiL, NeL+1:N2] = 0.
  return W
end


function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  return t
end

function molda_euler_lif_CSR(h, total, N, IFRAC, W, CSR, fe1, fi1, fe2, fi2, vth, tau_m, tau_s, tau_a, g_a)

  #Divide each chunk of the network
  N2 = Int64(round(N/2))
  NiL = Int64(round(N2/IFRAC))
  NeL = Int64(N2-NiL)
  Ne2 = NeL*2
  Ni2 = NiL*2

  #Feedforward input
  drive = zeros(N)
  drive[1:NeL] .+= fe1 * h
  drive[NeL+1:Ne2] .+= fe2 * h
  drive[Ne2+1:Ne2+NiL] .+= fi1 * h
  drive[Ne2+NiL+1:N] .+= fi2 * h

  #time
  ntotal = round(Int64, total/h)
  time = Float64[0]
  raster = Float64[0]

  #State Variables
  V = rand(N)*vth
  V_buff = V
  A = zeros(N)
  syn = zeros(N)

  #Compute leak terms once in advance (state_variable = state_variable - state_variable*h/tau)
  m_leak = h/tau_m
  s_leak = h/tau_s
  a_leak = h/tau_a
  g_h = g_a * h

  kill_flag = false

  for iter = 1:ntotal

      V .+= drive .+ (h*syn) .- (A*g_h) .- (V*m_leak)
      syn .-= (syn .* s_leak)
      A .-= (A .* a_leak)

      vs = (V .> vth)
      vsm = sum(vs)

      if vsm > 0
        sp = find(vs)
        for j = 1:vsm
          js = sp[j]
          delta_h = interpolate_spike(V[js], V_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn[CSR[js]] .+= W[CSR[js], js] .* lx
          if js <= Ne2
            A[js] .+= 1
          end
          push!(raster, js)
          push!(time, iter-delta_h)
        end
      end

      V .-= vth*vs
      V_buff = V

    end
    return time, raster
  end
