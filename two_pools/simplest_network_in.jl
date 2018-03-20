#N_in is the number of pre-synaptic neurons
#N_in*width is the fraction of those neurons synapsing onto you
#Now you have an array where the first N_in*width are populated by single-neuron strengths
#reduce to connection probability p on that domain
#now circularly shift to center the incoming connection footprint

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

#number of neurons, connection density, synaptic strengths, percentage of neurons in the pool you connect to
#width can be a max of .5 for wie before you get accidental overlap
function homogenous_weights(N, p, Aee, Aei, Aie, Aie_NL, Aii, width)

  Ni = Int64(round(N/5))
  Ne = N - Ni

  weed = top_hat_dist(Ne, width, Aee, p)
  weid = -top_hat_dist(Ni, width, Aei, p)
  wied = top_hat_dist(Ne, width, Aie, p/2.) + circshift(top_hat_dist(Ne, width, Aie_NL, p/2.), Int64(round(Ne/2)))
  wiid = -top_hat_dist(Ni, width, Aii, p)

  W = zeros(N,N)
  flat = []
  for i = 1:Ne
    W[i, 1:Ne] = top_hat_dist(Ne, width, Aee, p, i-1)
    W[i, Ne+1:end] = -top_hat_dist(Ni, width, Aei, p, div(Ni*i, Ne)-1)
  end

  for i = 1:Ni
    W[Ne+i, 1:Ne] = top_hat_dist(Ne, width, Aie, p/2., div(Ne*i, Ni)-1) + top_hat_dist(Ne, width, Aie_NL, p/2., Int64(round(Ne/2)+div(Ne*i, Ni)-1))
    W[Ne+i, Ne+1:end] = -top_hat_dist(Ni, width, Aii, p, i-1)
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

function Simple_Network_CSR(h, total, CSR, W, N, s1, s2, vth, tau_m, tau_s)

  Ni = Int64(round(N/5))
  Ne = N - Ni

  ntotal = round(Int,total/h)
  half = round(Int64, N/2)
  half_ = half + 1
  quar = round(Int64, N/4)
  quar_ = quar + 1

  V = rand(N)*vth
  V_buff = V
  drive = zeros(N)
  FPe = div(Ne,6) #each pool is 2/5 of the network
  P1s = round(Int,Ne/4-FPe)
  P1e = round(Int,Ne/4+FPe)
  P2s = round(Int,3*Ne/4-FPe)
  P2e = round(Int,3*Ne/4+FPe)
  drive[P1s:P1e] = s1*h
  drive[P2s:P2e] = s2*h
  syn = zeros(N)

  #raster
  time = Float64[0]
  raster = Float64[0]

  #Calculate drive and leak constants ahead of time
  m_leak = h/tau_m
  s_leak = h/tau_s

  println(ntotal, " beginning time loop")

  for iter = 1:ntotal

    V += (h*syn) - (V*m_leak) + drive
    syn -= (syn*s_leak)

    #check for spikes
    vs = (V.>vth)
    vsm = sum(vs)

    #update excitatory synapses and adaptation
    if vsm > 0
      spe = find(vs)
      for j = 1:vsm
        js = spe[j]
        #interpolate spike time
        delta_h = interpolate_spike(V[js], V_buff[js], vth) #time since the spike (units of h)
        lx = exp(delta_h/tau_m)
        syn[CSR[js]] += W[CSR[js], js].*lx #modify amplitude of synapse by decay since estimated time of spike
        push!(raster,js)
  	    push!(time,iter-delta_h)
      end
    end

    #reset after spike
    V -= vth*vs
    V_buff = V
  end

  return time, raster
end


function emptiness(x, funk, error_code) #checks if an empty set got populated
  if length(x) > 0
    return funk(x)
  else
    return error_code
  end
end

function CV_ISI_ALLTIME(Neurons, t, r)
  if Neurons == -5
    return -5, -5, -5
  else
  CVS = []
  for i in eachindex(Neurons)
    INT = t[find(r.==Neurons[i])]
    isi = diff(INT)
    c = cv(isi)
    if isnan(c) == false
      push!(CVS, c)
    end
  end
  CVS = convert(Array{Float64}, CVS)
  mean_cv = emptiness(CVS, mean, -5)
  median_cv = emptiness(CVS, median, -5)
  std_cv = emptiness(CVS, std, -5)
  return mean_cv, median_cv, std_cv
end
end

function score_analysis(re, N)
  #1 = top, 2 = bot
  half = div(N, 2)
  a1 = length(find(re .<= half))
  a2 = length(find(re .> half))
  biggie = max(a1, a2)
  spike_asymmetry = biggie/length(re)
  if a1 > a2
    bias = 2
  else
    bias = 1
  end
  return spike_asymmetry, bias
end

function nt_counts(t, r, bin_size, ntotal)

  ntf_bins = collect(minimum(t):bin_size:ntotal)
  ntf_bin_num = length(ntf_bins)
  ntf_nt = zeros(ntf_bin_num-1)
  nt = zeros(length(ntf_nt))

  for j = 2:ntf_bin_num
    nt[j-1] = sum(ntf_bins[j-1] .<= t .< ntf_bins[j]) #network counts
  end

  return nt
end


function network_fano(t, r, bin_size, ntotal)

  ntf_bins = collect(minimum(t):bin_size:ntotal)
  ntf_bin_num = length(ntf_bins)
  ntf_nt = zeros(ntf_bin_num-1)
  nt = zeros(length(ntf_nt))

  for j = 2:ntf_bin_num
    nt[j-1] = sum(ntf_bins[j-1] .<= t .< ntf_bins[j]) #network counts
  end

  fano = FANO(nt)
  return fano
end

function cv(isi)
  SD = std(isi)
  return SD/mean(isi)
end

function FANO(aRAY)
  return var(aRAY)/mean(aRAY)
end

function analyze_inputs(s)
  m = mean(s)
  SD = std(s)
  t = (m-1.)/SD
  return m, SD, t
end


function Neuron_finder(r, ns, mini)
  Neurons = collect(Set{Float64}(r))
  N = []
  #S = []
  for i in Neurons
    a = find(r .== i)
    if length(a) > ns
      push!(N, i)
      #push!(S, length(a))
    end
  end
  if length(N) < mini
    NF = -5
    #SF = -5
  else
    NF = convert(Array{Float64}, N)
    #SF = convert(Array{Int64}, S)
  end
  return NF
end

function Neurons_tb_ns(r, half, ns, mini)
  Neurons = collect(Set{Float64}(r))
  f1, f2 = find(Neurons .> half), find(Neurons .<= half)
  N1, N2 = Neurons[f1], Neurons[f2]
  TN, BN = [], [] #top neurons and bottom neurons; who spiked > ns times
  for i in N1
    if length(find(r .== i)) > ns
      push!(TN, i)
    end
  end

  for i in N2
    if length(find(r .== i)) > ns
      push!(BN, i)
    end
  end


  if ((length(TN) < mini) & (length(BN) < mini))
    return -5, -5
  elseif ((length(TN) < mini) & (length(BN) > mini))
    return -5, convert(Array{Float64}, BN)
  elseif ((length(BN) < mini) & (length(TN) > mini))
    return convert(Array{Float64}, TN), -5
  else
  tn = convert(Array{Float64}, TN)
  bn = convert(Array{Float64}, BN)
  return tn, bn
end
end

function count_train(bin, lt, lr, Neurons, n)
  if Neurons == -5
    return -5
  else
  bins = collect(lt[1]:bin:lt[end])
  count1 = zeros(n, length(bins)-1)
  for i in eachindex(Neurons)
    INT = lt[find(lr.==Neurons[i])]
    for j = 1:length(bins)-1
      count1[i,j] = length(find(bins[j] .<= INT .< bins[j+1]))
    end
  end
  return count1
end
end

function fano_train(Fcount, error_code)
  if Fcount == -5
    return -5, -5, -5
  else
  x = []
  for i = 1:size(Fcount)[1]
    f1 = FANO(Fcount[i,:])
    if isnan(f1) == false
      push!(x, f1)
    end
  end
  #println(x)
  a = convert(Array{Float64}, x)
  #println(a)
  meanf = emptiness(a, mean, error_code)
  medif = emptiness(a, median, error_code)
  stdvf = emptiness(a, std, error_code)
  return meanf, medif, stdvf
end
end

function rand_pair_cor(bin, lt, lr, Neurons, n)
  if Neurons == -5
    return -5
  end
  bins = collect(lt[1]:bin:lt[end])
  ya = [] #if you already looked at this neuron
  neurons = shuffle(Neurons) #
  lank = length(Neurons)
  shank = zeros(lank, length(bins)-1) #store count trains as you go
  cor_store = []
  #loop through n random pairwise correlations
  for i = 1:n
    #pick two random neurons in this pool
    x1 = rand(1:lank)
    x2 = rand(1:lank)
    if x1 == x2
      while x1 == x2
        x2 = rand(1:lank)
      end
    end
    #conditional statements just ensure that you don't waste time re-calculating count trains that you already have
    if ((x1 in ya) & (x2 in ya)) #already calculated count trains for both neurons
      c = cor(vec(shank[x1,:]), vec(shank[x2,:])) #go ahead and correlate them
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif  ((x1 in ya) & ((x2 in ya) == false)) #if you already did one but not the other
      #update shank for the new neuron
      INT2 = lt[find(lr.==neurons[x2])]
      for j = 1:length(bins)-1
        shank[x2,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x2)
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif ((x2 in ya) & ((x1 in ya) == false)) #if you already did the other but not the one
      #update shank for the one neuron
      INT2 = lt[find(lr.==neurons[x1])]
      for j = 1:length(bins)-1
        shank[x1,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x1)
      if isnan(c) == false
        push!(cor_store, c)
      end
    else #you don't have either
      #get count trains
      INT1 = lt[find(lr.==neurons[x1])]
      INT2 = lt[find(lr.==neurons[x2])]
      #add them to shank so you don't have to recalculate
      for j = 1:length(bins)-1
        shank[x1,j] = length(find(bins[j] .<= INT1 .< bins[j+1]))
        shank[x2,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      #update the list of counts you already have
      push!(ya, x1)
      push!(ya, x2)
    end
  end
  a = convert(Array{Float64}, cor_store)
  #now return the mean of these correlations
  return a
end

function moments(x)
  n = length(x)
  m = mean(x)
  skew = zeros(n)
  kurt = zeros(n)
  vary = zeros(n)
  for i =1:n
    d = x[i] - m
    skew[i] = d^3
    kurt[i] = d^4
    vary[i] = d^2
  end
  ms = mean(skew)
  mk = mean(kurt)
  mv = mean(vary)
  skewt = ms/(mv^1.5)
  KURT = mk/(mv^2)
  t = (m-1.)/sqrt(mv)
  return m, mv, t, skewt, KURT
end
