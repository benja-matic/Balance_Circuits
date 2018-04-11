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

#k is a fixed number of inputs for a local circuit
#
function homogenous_4x4_weights(N, IFRAC, k, Aee, Aei, Aie, Aie_NL, Aii)

  N_local = Int64(round(N/2)) #divide the network into two EI circuits
  Ni_local = Int64(round(N_local/IFRAC))
  Ne_local = Int64(N_local-Ni_local)
  Ne2 = Ne_local*2
  ks = sqrt(k)
  ke = k
  ki = k
  ke_local = round(Int64, k)
  ki_local = round(Int64, k)
  kie = round(Int64, k/2) #half from your local circuit, half from the other local circuit

  Jee = Aee/ks
  Jei = -Aei/ks
  Jie = Aie/ks
  Jie_NL = Aie_NL/ks
  Jii = -Aii/ks

  # Jee = Aee/ke_local
  # Jei = -Aei/ki_local
  # Jie = Aie/kie
  # Jie_NL = Aie_NL/kie
  # Jii = -Aii/ki_local

  W = zeros(N,N)

  ###EE AND EI FOR BOTH CIRCUITS
  for i = 1:Ne_local

    ee1_inds = rand(1:Ne_local, ke_local) #EE local
    ei1_inds = rand(Ne2+1:Ne2+Ni_local, ki_local) #IE local
    ee2_inds = rand(Ne_local+1:Ne2, ke_local) #other EE local
    ei2_inds = rand(N_local + Ne_local+1:N, ki_local) #other IE local

    for j in eachindex(ee1_inds) #pool1 EE
      W[i, ee1_inds[j]] += Jee
    end

    for j in eachindex(ee2_inds) #pool2 EE
      W[i+Ne_local, ee2_inds[j]] += Jee
    end

    for j in eachindex(ei1_inds) #pool 1 EI
      W[i, ei1_inds[j]] += Jei
    end

    for j in eachindex(ei2_inds) #pool 2 EI
      W[i+Ne_local, ei2_inds[j]] += Jei
    end

  end

  ###IE AND II FOR BOTH CIRCUITS
  for i = Ne2+1:Ne2+Ni_local

    ie1_inds = rand(1:Ne_local, kie)
    ii1_inds = rand(Ne2+1: Ne2+Ni_local, ki_local)
    ie2_inds = rand(Ne_local+1:Ne2, kie)
    ii2_inds = rand(Ne2+Ni_local+1:N, ki_local)

    ieNL1_inds = rand(Ne_local+1:Ne2, kie)
    ieNL2_inds = rand(1:Ne_local, kie)

    for j in eachindex(ie1_inds)
      W[i, ie1_inds[j]] += Jie
    end

    for j in eachindex(ii1_inds)
      W[i, ii1_inds[j]] += Jii
    end

    for j in eachindex(ie2_inds)
      W[i+Ni_local, ie2_inds[j]] += Jie
    end

    for j in eachindex(ii2_inds)
      W[i+Ni_local, ii2_inds[j]] += Jii
    end

    for j in eachindex(ieNL1_inds)
      W[i, ieNL1_inds[j]] += Jie_NL
    end

    for j in eachindex(ieNL2_inds)
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

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  return t
end

function euler_lif_CSR_4x4_s(h, total, N, IFRAC, W, CSR, fe, fi, vth, tau_m, tau_s, tau_a, g_a)

  N2 = Int64(round(N/2)) #divide the network into two EI circuits
  NiL = Int64(round(N2/IFRAC))
  NeL = Int64(N2-NiL)
  Ne2 = NeL*2
  Ni2 = NiL*2
  drive = zeros(Ne2) .+ fe*h
  drivi = zeros(Ni2) .+ fi*h
  # drive[1:NeL] .+= s*h

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

  m_leak = h/tau_m
  s_leak = h/tau_s
  a_leak = h/tau_a

  kill_flag = false

  for iter = 1:ntotal

      Ve .+= drive .+ h*syn_ee .+ h*syn_ei .- (Ae*g_a) .- (Ve*m_leak)
      Vi .+= drivi .+ h*syn_ie .+ h*syn_ieL .+ h*syn_ii .- (Vi*m_leak)
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
