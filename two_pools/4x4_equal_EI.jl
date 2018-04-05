function top_hat_dist(N_in, width, A, p, circ)
  incoming = zeros(N_in)
  n_con = Int64(round(width*N_in))
  nc2 = Int64(round(n_con/2.))
  x = rand(n_con) .<= p
  f = zeros(n_con) .+ A/(n_con*p)
  f .*= x
  incoming[1:n_con] = f
  return circshift(incoming, circ-nc2)
end

function homogenous_4x4_weights(N, p, Aee, Aei, Aie, Aie_NL, Aii)

  N_local = Int64(round(N/2)) #divide the network into two EI circuits
  half = round(Int64, N_local/2)
  k_local = round(Int64, half*p)
  kie = round(Int64, half*(p/2)) #half from your local circuit, half from the other local circuit

  Jee = Aee/k_local
  Jei = -Aei/k_local
  Jie = Aie/kie
  Jie_NL = Aie_NL/kie
  Jii = -Aii/k_local

  W = zeros(N,N)

  ###EE AND EI FOR BOTH CIRCUITS
  for i = 1:half

    ee1_inds = rand(1:half, k_local) #EE local
    ei1_inds = rand(N_local+1:N_local+half, k_local) #IE local

    ee2_inds = rand(half+1:N_local, k_local) #other EE local
    ei2_inds = rand(N_local+half+1:N, k_local) #other IE local

    ie1_inds = rand(1:half, kie)
    ie2_inds = rand(half+1:N_local, kie)

    ieL1_inds = rand(half+1:N_local, kie)
    ieL2_inds = rand(1:half, kie)

    ii1_inds = rand(N_local+1:N_local+half, k_local)
    ii2_inds = rand(N_local+half+1:N, k_local)

    for j in eachindex(ee1_inds) #pool1 EE
      W[i, ee1_inds[j]] += Jee
    end

    for j in eachindex(ee2_inds) #pool2 EE
      W[i+half, ee2_inds[j]] += Jee
    end

    for j in eachindex(ei1_inds) #pool1 EI
      W[i, ei1_inds[j]] += Jei
    end

    for j in eachindex(ei2_inds) #pool2 EI
      W[i+half, ei2_inds[j]] += Jei
    end

    for j in eachindex(ie1_inds) #pool1 IE
      W[i+N_local, ie1_inds[j]] += Jie
    end

    for j in eachindex(ie2_inds) #pool2 IE
      W[i+N_local+half, ie2_inds[j]] += Jie
    end

    for j in eachindex(ieL1_inds) #pool1 IE Non_Local
      W[i+N_local, ieL1_inds[j]] += Jie_NL
    end

    for j in eachindex(ieL2_inds) #pool2 IE Non_Local
      W[i+N_local+half, ieL2_inds[j]] += Jie_NL
    end

    for j in eachindex(ii1_inds) #pool1 II
      W[i+N_local, ii1_inds[j]] += Jii
    end

    for j in eachindex(ii2_inds) #pool2 II
      W[i+N_local+half, ii2_inds[j]] += Jii
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

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  return t
end

function euler_lif_CSR_4x4_s(h, total, N, W, CSR, s, vth, tau_m, tau_s, tau_a, g_a)

  N2 = div(N, 2)
  NeL = round(Int64, N2/2)
  Ne2 = NeL*2
  NiL = round(Int64, N2/2)
  Ni2 = NiL*2
  drive = zeros(Ne2) .+ s*h

  ntotal = round(Int64, total/h)
  #Voltage
  Ve = rand(Ne2)*vth
  Vi = rand(Ni2)*vth
  Ve_buff = Ve
  Vi_buff = Vi

  #Store synaptic inputs to infer theoretical parameters
  SEE = zeros(100, ntotal)
  SEI = zeros(100, ntotal)
  SIE = zeros(100, ntotal)
  SIEL = zeros(100, ntotal)
  SII = zeros(100, ntotal)

  #incoming synapse variables
  syn_ee = zeros(Ne2)
  syn_ei = zeros(Ne2)
  syn_ie = zeros(Ni2)
  syn_ieL = zeros(Ni2)
  syn_ii = zeros(Ni2)
  #Adaptation
  Ae = zeros(Ne2)

  #store synaptic inputs

  time_e = Float64[0]
  raster_e = Float64[0]
  time_i = Float64[0]
  raster_i = Float64[0]
  # LX = []

  m_leak = h/tau_m
  s_leak = h/tau_s
  a_leak = h/tau_a

  kill_flag = false

  for iter = 1:ntotal

      Ve .+= drive .+ h*syn_ee .+ h*syn_ei .- (Ae*g_a) .- (Ve*m_leak)
      Vi .+= h*syn_ie .+ h*syn_ieL .+ h*syn_ii .- (Vi*m_leak)
      #Store synaptic inputs
      SEE[1:50,iter] = syn_ee[1:50]#*h
      SEE[51:100, iter] = syn_ee[NeL+1:NeL+50]#*h

      SEI[1:50,iter] = syn_ei[1:50]#*h
      SEI[51:100, iter] = syn_ei[NeL+1:NeL+50]#*h

      SIE[1:50,iter] = syn_ie[1:50]#*h
      SIE[51:100, iter] = syn_ie[NiL+1:NiL+50]#*h

      SIEL[1:50,iter] = syn_ieL[1:50]#*h
      SIEL[51:100, iter] = syn_ieL[NiL+1:NiL+50]#*h

      SII[1:50,iter] = syn_ii[1:50]#*h
      SII[51:100, iter] = syn_ii[NiL+1:NiL+50]#*h

      syn_ee .-= (syn_ee .* s_leak)
      syn_ei .-= (syn_ei .* s_leak)
      syn_ie .-= (syn_ie .* s_leak)
      syn_ieL .-= (syn_ieL .* s_leak)
      syn_ii .-= (syn_ii .* s_leak)
      Ae .-= (Ae*a_leak)

      ves = (Ve .> vth)
      vesm = sum(ves)

      if vesm > 0
        spe = find(ves)
        for j = 1:vesm
          js = spe[j]
          delta_h = interpolate_spike(Ve[js], Ve_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn_ee .+= W[1:Ne2, js] .* lx
          if js <= NeL
            syn_ie[1:NiL] .+= W[Ne2+1:Ne2+NiL, js] .* lx
            syn_ieL[NiL+1:Ni2] .+= W[Ne2+NiL+1:N, js] .* lx
          else
            syn_ie[NiL+1:Ni2] .+= W[Ne2+NiL+1:N, js] .* lx
            syn_ieL[1:NiL] .+= W[Ne2+1:Ne2+NiL, js] .* lx
          end
          Ae .+= 1
          push!(raster_e, js)
          push!(time_e, iter-delta_h)
        end
      end

      vis = (Vi .> vth)
      vism = sum(vis)

      if vism > 0
        spi = find(vis)
        for j = 1:vism
          js = spi[j]
          delta_h = interpolate_spike(Vi[js], Vi_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn_ii .+= W[Ne2+1:N, js+Ne2] .* lx
          syn_ei .+= W[1:Ne2, js+Ne2] .* lx
          push!(raster_i, js)
          push!(time_i, iter-delta_h)
        end
      end

      Ve .-= vth*ves
      Ve_buff = Ve
      Vi .-= vth*vis
      Vi_buff = Vi

    end
    return time_e, raster_e, time_i, raster_i, SEE, SEI, SIE, SIEL, SII
  end
