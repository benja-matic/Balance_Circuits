function local_random_2x2(N, IFRAC, k, Aee, Aei, Aie, Aii)

  Ni = Int64(round(N/IFRAC))
  Ne = N - Ni
  ks = sqrt(k)

  Jee = Aee/ks
  Jei = -Aei/ks
  Jie = Aie/ks
  Jii = -Aii/ks

  W1 = zeros(N,N)
  for i = 1:Ne
    e_inds = rand(1:Ne, k)
    i_inds = rand(Ne+1:N, k)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jee
      W1[i, i_inds[j]] += Jei
    end
  end

  for i = Ne+1:N
    e_inds = rand(1:Ne, k)
    i_inds = rand(Ne+1:N, k)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jie
      W1[i, i_inds[j]] += Jii
    end
  end

  return W1
end

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
  println(size(W1))
  for i = 1:Ne2
    e_inds = rand(1:Ne2, k2)
    i_inds = rand(Ne + 1:Ne + Ni2, k2)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jee
      W1[i, e_inds[j] + Ne2] += Jee
      W1[i, i_inds[j]] += Jei
      W1[i, i_inds[j] + Ni2] += Jei

      W1[i+Ne2, e_inds[j]] += Jee
      W1[i+Ne2, e_inds[j] + Ne2] += Jee
      W1[i+Ne2, i_inds[j]] += Jei
      W1[i+Ne2, i_inds[j] + Ni2] += Jei
    end
  end

  for i = Ne+1:Ne+Ni2
    e_inds = rand(1:Ne2, k2)
    i_inds = rand(Ne + 1:Ne + Ni2, k2)
    for j in eachindex(e_inds)
      W1[i, e_inds[j]] += Jie
      W1[i, e_inds[j] + Ne2] += Jie
      W1[i, i_inds[j]] += Jii
      W1[i, i_inds[j] + Ni2] += Jii

      W1[i+Ni2, e_inds[j]] += Jie
      W1[i+Ni2, e_inds[j] + Ne2] += Jie
      W1[i+Ni2, i_inds[j]] += Jii
      W1[i+Ni2, i_inds[j] + Ni2] += Jii
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

function euler_lif_2x2_CSR(h, total, N, IFRAC, W, CSR, fe1, fi1, fe2, fi2, vth, tau_m, tau_s)

  Ni = Int64(round(N/IFRAC))
  Ne = N - Ni
  Ne2 = round(Int64, Ne/2)
  Ni2 = round(Int64, Ni/2)

  drive_e = zeros(Ne)
  drive_i = zeros(Ni)

  drive_e[1:Ne2] = fe1 * h
  drive_e[Ne2+1:end] = fe2 * h
  drive_i[1:Ni2] = fi1 * h
  drive_i[Ni2+1:end] = fi2 * h

  ntotal = round(Int64, total/h)

  Ve = rand(Ne) * vth
  Vi = rand(Ni) * vth
  Ve_buff = Ve
  Vi_buff = Vi

  syn_ee = zeros(Ne)
  syn_ei = zeros(Ne)
  syn_ie = zeros(Ni)
  syn_ii = zeros(Ni)

  # SEE = zeros(100, ntotal)
  # SEI = zeros(100, ntotal)
  # SIE = zeros(100, ntotal)
  # SII = zeros(100, ntotal)

  SEE_L = zeros(40, ntotal)
  SEE_NL = zeros(40, ntotal)
  SEI_L = zeros(40, ntotal)
  SEI_NL = zeros(40, ntotal)
  SIE_L = zeros(40, ntotal)
  SIE_NL = zeros(40, ntotal)
  SII_L = zeros(40, ntotal)
  SII_NL = zeros(40, ntotal)


  time_e = Float64[0]
  raster_e = Float64[0]
  time_i = Float64[0]
  raster_i = Float64[0]

  m_leak = h/tau_m
  s_leak = h/tau_s

  kill_flag = false

  for iter = 1:ntotal

      Ve .+= drive_e .+ (h * syn_ee) .+ (h * syn_ei) .- (Ve * m_leak)
      Vi .+= drive_i .+ (h * syn_ie) .+ (h * syn_ii) .- (Vi * m_leak)

      syn_ee .-= (syn_ee * s_leak)
      syn_ie .-= (syn_ie * s_leak)
      syn_ei .-= (syn_ei * s_leak)
      syn_ii .-= (syn_ii * s_leak)

      # SEE[:,iter] = syn_ee[1:100]
      # SEI[:,iter] = syn_ei[1:100]
      # SIE[:,iter] = syn_ie[1:100]
      # SII[:,iter] = syn_ii[1:100]

      # SEE[1:50,iter] = syn_ee[1:50]#*h
      # SEE[51:100, iter] = syn_ee[Ne2+1:Ne2+50]#*h
      #
      # SEI[1:50,iter] = syn_ei[1:50]#*h
      # SEI[51:100, iter] = syn_ei[Ne2+1:Ne2+50]#*h
      #
      # SIE[1:50,iter] = syn_ie[1:50]#*h
      # SIE[51:100, iter] = syn_ie[Ni2+1:Ni2+50]#*h
      #
      # SII[1:50,iter] = syn_ii[1:50]#*h
      # SII[51:100, iter] = syn_ii[Ni2+1:Ni2+50]#*h

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
          syn_ei .+= W[1:Ne, js+Ne] .* lx
          syn_ii .+= W[Ne+1:N, js+Ne] .* lx
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
