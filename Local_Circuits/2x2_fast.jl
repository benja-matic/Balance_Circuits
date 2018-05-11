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

function euler_lif_2x2_CSR(h, total, N, IFRAC, W, CSR, fe1, fi1, fe2, fi2, vth, tau_m, tau_s)

  Ni = Int64(round(N/IFRAC))
  Ne = N - Ni
  Ne2 = round(Int64, Ne/2)
  Ni2 = round(Int64, Ni/2)
  drive = zeros(N)
  drive[1:Ne2] = fe1 * h
  drive[Ne2+1:Ne] = fe2 * h
  drive[Ne+1:Ne+Ni2] = fi1 * h
  drive[Ne+Ni2+1:N] = fi2 * h

  ntotal = round(Int64, total/h)

  V = rand(N)*vth
  V_buff = V
  syn = zeros(N)

  time = Float64[0]
  raster = Float64[0]

  m_leak = h/tau_m
  s_leak = h/tau_s

  kill_flag = false

  for iter = 1:ntotal

      V .+= drive .+ (h*syn) .- (V*m_leak)
      syn .-= (syn .* s_leak)

      vs = (V .> vth)
      vsm = sum(vs)

      if vsm > 0
        sp = find(vs)
        for j = 1:vsm
          js = sp[j]
          delta_h = interpolate_spike(V[js], V_buff[js], vth)
          lx = exp(-delta_h*h/tau_s)
          syn[CSR[js]] .+= W[CSR[js], js] .* lx
          push!(raster, js)
          push!(time, iter-delta_h)
        end
      end

      V .-= vth*vs
      V_buff = V

    end
    return time, raster
  end
