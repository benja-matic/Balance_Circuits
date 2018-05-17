function local_random_2x2_symmetric(N, IFRAC, k, Aee, Aei, Aie, Aii)

  Ni = Int64(round(N/IFRAC))
  Ne = N - Ni
  ks = sqrt(k)
  k2 = round(Int64, k/2)
  Ne2 = round(Int64, Ne/2)
  Ni2 = round(Int64, Ni/2)

  Jee = Aee/ks
  Jei = -Aei/ks
  Jie = Aie/ks
  Jii = -Aii/ks

  W1 = zeros(N,N)
  for i = 1:Ne
    e_inds = rand(1:Ne2, k2)
    i_inds = rand(Ne + 1:Ne + Ni2, k2)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jee
      W1[i, e_inds[j] + Ne2] += Jee
      W1[i, i_inds[j]] += Jei
      W1[i, i_inds[j] + Ni2] += Jei
    end
  end

  for i = Ne+1:N
    e_inds = rand(1:Ne2, k2)
    i_inds = rand(Ne + 1:Ne + Ni2, k2)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jie
      W1[i, e_inds[j] + Ne2] += Jie
      W1[i, i_inds[j]] += Jii
      W1[i, i_inds[j] + Ni2] += Jii
    end
  end

  return W1
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

# function find_spikes_update_synapses(v, v_buff, W, compartments, s1, s2, s3, s4)
#   #find spikes
#   vs = (v .> vth)
#   vsm = sum(vs)
#   #update synapses
#   if vesm1 > 0
#     sp = find(vsm)
#     for j = 1:vesm1
#       js = spe[j]
#       delta_h = interpolate_spike(v[js], v_buff[js], vth)
#       lx = exp(-delta_h*h/tau_s)
#       s1 .+= W[1:Ne, js] .* lx
#       s2 .+= W[1:Ne, js] .* lx
#       s3 .+= W[Ne+1:N, js] .* lx
#       s4 .+= W[Ne+1:N, js] .* lx
#       push!(raster_e, js)
#       push!(time_e, iter-delta_h)
#     end
#   end


function euler_lif_2x2_CSR(h, total, N, IFRAC, W, CSR, fe1, fi1, fe2, fi2, vth, tau_m, tau_s)

  #set up all state variables, storage

  #how many neurons/group
  Ni = Int64(round(N/IFRAC))
  Ne = N - Ni
  Ne2 = round(Int64, Ne/2)
  Ni2 = round(Int64, Ni/2)
  #set up feedforward inputs
  drive_e1 = zeros(Ne2) .+ (fe1 * h)
  drive_i1 = zeros(Ni2) .+ (fe2 * h)
  drive_e2 = zeros(Ne2) .+ (fi1 * h)
  drive_i2 = zeros(Ni2) .+ (fi2 * h)
  #simulation time in units of dt
  ntotal = round(Int64, total/h)
  #voltage variables
  Ve1 = rand(Ne2) * vth
  Ve2 = rand(Ne2) * vth
  Vi1 = rand(Ni2) * vth
  Vi2 = rand(Ni2) * vth
  #voltage buffer for spike time interpolation
  Ve1_buff = Ve1
  Vi1_buff = Vi1
  Ve2_buff = Ve2
  Vi2_buff = Vi2
  #synapse variables (local)
  syn_ee_L = zeros(Ne)
  syn_ei_L = zeros(Ne)
  syn_ie_L = zeros(Ni)
  syn_ii_L = zeros(Ni)
  #synapse variables (non-local)
  syn_ee_NL = zeros(Ne)
  syn_ei_NL = zeros(Ne)
  syn_ie_NL = zeros(Ni)
  syn_ii_NL = zeros(Ni)
  #spikes times and neuron indices
  time_e1 = Float64[0]
  raster_e1 = Float64[0]
  time_i1 = Float64[0]
  raster_i1 = Float64[0]
  #second pool
  time_e2 = Float64[0]
  raster_e2 = Float64[0]
  time_i2 = Float64[0]
  raster_i2 = Float64[0]
  #calculate leak terms once in advance
  m_leak = h/tau_m
  s_leak = h/tau_s

  #start integrating

  for iter = 1:ntotal
      #update voltage variables
      #pool E1
      Ve1 .+= drive_e1 .+ (h * syn_ee_L[1:Ne2][:]) .+ (h * syn_ee_NL[1:Ne2][:]) .+ \
      (h * syn_ei_L[1:Ne2][:]) .+ (h * syn_ei_NL[1:Ne2][:]) .- (Ve1 * m_leak)
      #Pool E2
      Ve2 .+= drive_e2 .+ (h * syn_ee_L[Ne2+1:end][:]) .+ (h * syn_ee_NL[Ne2+1:end][:]) .+ \
      (h * syn_ei_L[Ne2+1:end][:]) .+ (h * syn_ei_NL[Ne2+1:end][:]) .- (Ve2 * m_leak)
      #Pool I1
      Vi1 .+= drive_i1 .+ (h * syn_ie_L[1:Ni2][:]) .+ (h * syn_ie_NL[1:Ni2][:]) .+ \
      (h * syn_ii_L[1:Ni2][:]) .+ (h * syn_ii_NL[1:Ni2][:]) .- (Vi1 * m_leak)
      #Pool I2
      Vi2 .+= drive_i2 .+ (h * syn_ie_L[Ni2+1:end][:]) .+ (h * syn_ie_NL[Ni2+1:end][:]) .+ \
      (h * syn_ii_L[Ni2+1:end][:]) .+ (h * syn_ii_NL[Ni2+1:end][:]) .- (Vi2 * m_leak)
      #leak local synapses
      syn_ee_L .-= (syn_ee_L * s_leak)
      syn_ie_L .-= (syn_ie_L * s_leak)
      syn_ei_L .-= (syn_ei_L * s_leak)
      syn_ii_L .-= (syn_ii_L * s_leak)
      #leak non-local synapses
      syn_ee_NL .-= (syn_ee_NL * s_leak)
      syn_ie_NL .-= (syn_ie_NL * s_leak)
      syn_ei_NL .-= (syn_ei_NL * s_leak)
      syn_ii_NL .-= (syn_ii_NL * s_leak)

      ###start counting spikes

      ves1 = (Ve1 .> vth)
      vesm1 = sum(ves1)
      if vesm1 > 0
        spe = find(ves1)
        for j = 1:vesm1
          js = spe[j]
          delta_h = interpolate_spike(Ve1[js], Ve1_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn_ee_L[1:Ne2] .+= W[1:Ne2, js] .* lx
          syn_ee_NL[1:Ne2] .+= W[Ne2+1:Ne, js] .* lx
          syn_ie_L[1:Ni2] .+= W[Ne+1:Ne+Ni2, js] .* lx
          syn_ie_NL[1:Ni2] .+= W[Ne+Ni2+1:N, js] .* lx
          push!(raster_e1, js)
          push!(time_e1, iter-delta_h)
        end
      end

      ves2 = (Ve2 .> vth)
      vesm2 = sum(ves2)

      if vesm2 > 0
        spe = find(ves2)
        for j = 1:vesm2
          js = spe[j]
          delta_h = interpolate_spike(Ve2[js], Ve2_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn_ee_L[Ne2+1:end] .+= W[Ne2+1:Ne, js + Ne2] .* lx
          syn_ee_NL[Ne2+1:end] .+= W[1:Ne2, js + Ne2] .* lx
          syn_ie_L[Ni2+1:end] .+= W[Ne+Ni2+1:N, js + Ne2] .* lx
          syn_ie_NL[Ni2+1:end] .+= W[Ne+1:Ne+Ni2, js + Ne2] .* lx
          push!(raster_e2, js)
          push!(time_e2, iter-delta_h)
        end
      end

      vis1 = (Vi1 .> vth)
      vism1 = sum(vis1)

      if vism1 > 0
        spi = find(vis1)
        for j = 1:vism1
          js = spi[j]
          delta_h = interpolate_spike(Vi1[js], Vi1_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn_ei_L[1:Ne2] .+= W[1:Ne2, js + Ne] .* lx
          syn_ei_NL[1:Ne2] .+= W[Ne2+1:Ne, js + Ne] .* lx
          syn_ii_L[1:Ni2] .+= W[Ne+1:Ne+Ni2, js + Ne] .* lx
          syn_ii_NL[1:Ni2] .+= W[Ne+Ni2+1:N, js + Ne] .* lx
          push!(raster_i1, js)
          push!(time_i1, iter-delta_h)
        end
      end

      vis2 = (Vi2 .> vth)
      vism2 = sum(vis2)

      if vism2 > 0
        spi = find(vis2)
        for j = 1:vism2
          js = spi[j]
          delta_h = interpolate_spike(Vi2[js], Vi2_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn_ei_L[Ne2+1:end] .+= W[Ne2+1:Ne, js + Ne + Ni2] .* lx
          syn_ei_NL[Ne2+1:end] .+= W[1:Ne2, js + Ne + Ni2] .* lx
          syn_ii_L[Ni2+1:end] .+= W[Ne+Ni2+1:N, js + Ne + Ni2] .* lx
          syn_ii_NL[Ni2+1:end] .+= W[Ne+1:Ne+Ni2, js + Ne + Ni2] .* lx
          push!(raster_i2, js)
          push!(time_i2, iter-delta_h)
        end
      end

      #reset neurons who spiked

      Ve1 .-= vth*ves1
      Ve1_buff = Ve1

      Vi1 .-= vth*vis1
      Vi1_buff = Vi1

      Ve2 .-= vth*ves2
      Ve2_buff = Ve2

      Vi2 .-= vth*vis2
      Vi2_buff = Vi2

    end
    return time_e1, raster_e1, time_i1, raster_i1, time_e2, raster_e2, time_i2, raster_i2
  end
