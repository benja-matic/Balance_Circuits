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

function local_random_2x2(N, p, Aee, Aei, Aie, Aii)

  Ni = Int64(round(N/5))
  Ne = N - Ni
  ke = round(Int64, Ne*p)
  ki = round(Int64, Ni*p)
  kie = round(Int64, Ne*(p/2))

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


function homogenous_4x4_weights(N, p, Aee, Aei, Aie, Aie_NL, Aii)

  N_local = Int64(round(N/2)) #divide the network into two EI circuits
  Ni_local = Int64(round(N_local/5))
  Ne_local = Int64(N_local-Ni_local)
  Ne2 = Ne_local*2
  ke_local = round(Int64, Ne_local*p)
  ki_local = round(Int64, Ni_local*p)
  kie = round(Int64, Ne_local*(p/2)) #half from your local circuit, half from the other local circuit

  Jee = Aee/ke_local
  Jei = -Aei/ki_local
  Jie = Aie/kie
  Jie_NL = Aie_NL/kie
  Jii = -Aii/ki_local

  W = zeros(N,N)

  ###EE AND EI FOR BOTH CIRCUITS
  for i = 1:Ne_local

    ee1_inds = rand(1:Ne_local, ke_local) #EE local
    ie1_inds = rand(Ne2+1:Ne2+Ni_local, ki_local) #IE local
    ee2_inds = rand(Ne_local+1:Ne2, ke_local) #other EE local
    ie2_inds = rand(N_local + Ne_local+1:N, ki_local) #other IE local

    for j in eachindex(ee1_inds) #pool1 EE
      W[i, ee1_inds[j]] += Jee
    end

    for j in eachindex(ee2_inds) #pool2 EE
      W[i+Ne_local, ee2_inds[j]] += Jee
    end

    for j in eachindex(ie1_inds) #pool 1 EI
      W[i, ie1_inds[j]] += Jei
    end

    for j in eachindex(ie2_inds) #pool 2 EI
      W[i+Ne_local, ie2_inds[j]] += Jei
    end

  end

  ###IE AND II FOR BOTH CIRCUITS
  for i = Ne2+1:Ne2+Ni_local

    ie1_inds = rand(1:Ne_local, kie)
    ii1_inds = rand(Ne2+1: Ne2+Ni_local, ki_local)
    ie2_inds = rand(Ne_local:Ne2, kie)
    ii2_inds = rand(Ne2+Ni_local:N, ki_local)

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

function sparse_row(W,N)
  flat = []
  for i = 1:N
    mi = find(W[i,:])
    push!(flat, mi)
  end
  return flat
end

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  return t
end

function euler_lif_CSR_4x4(h, total, N, W, CSR, s, vth, tau_m, tau_s, tau_a, g_a)

  ntotal = round(Int64, total/h)
  syn = zeros(N) #synapse
  A = zeros(N) #adaptation
  V = rand(N)*vth #voltage
  V_buff = V
  Input = zeros(N, ntotal)

  time = Float64[0]
  raster = Float64[0]

  m_leak = h/tau_m
  s_leak = h/tau_s
  a_leak = h/tau_a

  drive = zeros(N)
  N2 = div(N, 2)
  NeL = div(4*N2, 5)
  Ne2 = NeL*2
  Ni = div(NeL, 4)
  Ni2 = Ni*2
  drive[1:Ne2] = s*h

  syn_e = zeros(Ne2, ntotal)
  syn_i = zeros(Ni2, ntotal)

  kill_flag = false

  for iter = 1:ntotal

      incoming = drive .+ (h .* syn) .- (A .* g_a)
      Input[:, iter] = incoming
      V .+= incoming .- (V .* m_leak)

      syn .-= (syn .* s_leak)
      A .-= A*a_leak

      vs = (V .> vth)
      vsm = sum(vs)

      if vsm > 0
        spe = find(vs)
        for j = 1:vsm
          js = spe[j]
          delta_h = interpolate_spike(V[js], V_buff[js], vth)
          lx = exp(-delta_h*h/tau_m)
          syn[CSR[js]] += W[CSR[js], js] .* lx
          push!(raster, js)
          push!(time, iter-delta_h)
          if spe[j] <= Ne2
            A[spe[j]] += 1.
          end
        end
      end

      V -= vth*vs
      V_buff = V

      # if iter % 777 == 0
      #   println(iter)
      # end

    end
    return time, raster, Input
  end
