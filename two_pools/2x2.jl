function local_random_2x2(N, IFRAC, p, Aee, Aei, Aie, Aii)

  Ni = Int64(round(N/IFRAC))
  Ne = N - Ni
  ke = round(Int64, Ne*p)
  ki = round(Int64, Ni*p)

  Jee = Aee/ke
  Jei = -Aei/ki
  Jie = Aie/ke
  Jii = -Aii/ki

  W1 = zeros(N,N)
  for i = 1:Ne
    e_inds = rand(1:Ne, ke)
    i_inds = rand(Ne+1:N, ki)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jee
    end
    for j in eachindex(i_inds)
      W1[i, i_inds[j]] += Jei
    end
  end

  for i = Ne+1:N
    e_inds = rand(1:Ne, ke)
    i_inds = rand(Ne+1:N, ki)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jie
    end
    for j in eachindex(i_inds)
      W1[i, i_inds[j]] += Jii
    end
  end

  return W1
end

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  return t
end

function euler_lif_2x2_s(h, total, N, IFRAC, W, s, vth, tau_m, tau_s)

  Ni = Int64(round(N/IFRAC))
  Ne = N - Ni
  drive = zeros(Ne) .+ s*h

  ntotal = round(Int64, total/h)
  #Voltage
  Ve = rand(Ne)*vth
  Vi = rand(Ni)*vth
  Ve_buff = Ve
  Vi_buff = Vi

  #Store synaptic inputs to infer theoretical parameters
  SEE = zeros(100, ntotal)
  SEI = zeros(100, ntotal)
  SIE = zeros(100, ntotal)
  SII = zeros(100, ntotal)

  #incoming synapse variables
  syn_ee = zeros(Ne)
  syn_ei = zeros(Ne)
  syn_ie = zeros(Ni)
  syn_ii = zeros(Ni)

  #store synaptic inputs

  time_e = Float64[0]
  raster_e = Float64[0]
  time_i = Float64[0]
  raster_i = Float64[0]

  m_leak = h/tau_m
  s_leak = h/tau_s

  kill_flag = false

  for iter = 1:ntotal

      Ve .+= drive .+ h*syn_ee .+ h*syn_ei .- (Ve*m_leak)
      Vi .+= h*syn_ie .+ h*syn_ii .- (Vi*m_leak)
      #Store synaptic inputs
      SEE[:,iter] = syn_ee[1:100]#*h
      SEI[:,iter] = syn_ei[1:100]#*h
      SIE[:,iter] = syn_ie[1:100]#*h
      SII[:,iter] = syn_ii[1:100]#*h
      syn_ee .-= (syn_ee .* s_leak)
      syn_ei .-= (syn_ei .* s_leak)
      syn_ie .-= (syn_ie .* s_leak)
      syn_ii .-= (syn_ii .* s_leak)

      ves = (Ve .> vth)
      vesm = sum(ves)

      if vesm > 0
        spe = find(ves)
        for j = 1:vesm
          js = spe[j]
          delta_h = interpolate_spike(Ve[js], Ve_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn_ee .+= W[1:Ne, js] .* lx
          syn_ie .+= W[Ne+1:N, js] .* lx
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
          syn_ii .+= W[Ne+1:N, js+Ne] .* lx
          syn_ei .+= W[1:Ne, js+Ne] .* lx
          push!(raster_i, js)
          push!(time_i, iter-delta_h)
        end
      end

      Ve .-= vth*ves
      Ve_buff = Ve
      Vi .-= vth*vis
      Vi_buff = Vi

    end
    return time_e, raster_e, time_i, raster_i, SEE, SEI, SIE, SII
  end
