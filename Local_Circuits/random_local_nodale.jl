###Local Circuit With Random Connections

function random_weights(N, p, thresh)
    k = p*N #each cell sees on average k inputs
    ks = sqrt(k)
    W = zeros(N,N)

    for i = 1:N
        W[i,:] = -1. + 2.5*randn(N) * thresh / ks
    end

    px = rand(N,N) .< p
    W = W .* px
    # W /= (p*N/thresh) #synapses shrink with number of connections, grow with threshold

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

function LIF_delta_solve_CSC(W, CSR, S, runtime, h)

    vth = 20
    ntotal = round(Int, runtime/h)
    S = fill(S, N) .* h
    v = rand(N)*vth
    tau_m = 20.
    #dV/dt = S - v/tau_m; Euler: V = V + dV/dt * h; Euler: V += S - V/tau_m * h
    #compute h/tau_m once in advance
    m_leak = h/tau_m
    t = []
    r = []
    input = zeros(N, ntotal)
    volta = zeros(N, ntotal)

    for iter = 1:ntotal

        volta[:, iter] = v
        v += S - v*m_leak
        ves = v .> 20.
        vsm = sum(ves)

        if vsm > 0
            spe = find(ves)
            for j = 1:vsm
                js = spe[j]
                v[CSR[js]] += W[CSR[js], js]
                input[CSR[js], iter] += W[CSR[js], js]
                push!(t, iter)
                push!(r, js)
            end
        end

        v -= 20.*ves
    end

    return t, r, input, volta
end
function LIF_delta_solve(W, S, runtime, h)

    vth = 20
    ntotal = round(Int, runtime/h)
    S = fill(S, N)
    v = rand(N)*vth
    tau_m = 20.
    #dV/dt = S - v/tau_m; Euler: V = V + dV/dt * h; Euler: V += S - V/tau_m * h
    #compute h/tau_m once in advance
    m_leak = h/tau_m
    t = []
    r = []
    input = zeros(N, ntotal)
    volta = zeros(N, ntotal)

    for iter = 1:ntotal

        volta[:, iter] = v
        v += S - v*m_leak
        ves = v .> 20.
        vsm = sum(ves)

        if vsm > 0
            spe = find(ves)
            for j = 1:vsm
                v += W[:, spe[j]]
                input[:, iter] += W[:,spe[j]]
                push!(r, spe[j])
                push!(t, iter)
            end
        end

        v -= 20.*ves
    end

    return t, r, input, volta
end

function cv(isi)
  SD = std(isi)
  return SD/mean(isi)
end

function FANO(aRAY)
  return var(aRAY)/mean(aRAY)
end

function tscore(mu, sigma)
  tsc = (mu-1.)/sigma
  return tsc
end

function emptiness(x, funk, error_code) #checks if an empty set got populated
  if length(x) > 0
    return funk(x)
  else
    return error_code
  end
end

function Neuron_finder(r, ns, mini)
  Neurons = collect(Set{Float64}(r))
  N = []
  for i in Neurons
    if length(find(r .== i)) > ns
      push!(N, i)
    end
  end

  if length(N) < mini
    NF = -5
  else
    NF = convert(Array{Float64}, N)
  end
  return NF
end

function count_train_intron(bin, lt, lr, Neurons, n, shuff)
  if shuff == true
    neurons = shuffle(Neurons)[1:n]
  else
    neurons = Neurons[1:n]
  end
  bins = collect(lt[1]:bin:lt[end])
  count1 = zeros(n, length(bins)-1)
  for i in eachindex(neurons)
    INT = lt[find(lr.==neurons[i])]
    for j = 1:length(bins)-1
      count1[i,j] = length(find(bins[j] .<= INT .< bins[j+1]))
    end
  end
  return count1
end

function rand_pair_cor(bin, lt, lr, Neurons, n)
  bins = collect(lt[1]:bin:lt[end])
  count1 = zeros(n, length(bins)-1)
  ya = []
  neurons = shuffle(Neurons)
  lank = length(Neurons)
  shank = zeros(lank, length(bins)-1)
  cor_store = []
  #loop through however many sample correlations you want
  for i = 1:n
    #pick two random neurons in this pool
    x1 = rand(1:lank)
    x2 = rand(1:lank)
    #conditional statements just ensure that you don't waste time re-calculating count trains that you already have
    if ((x1 in ya) & (x2 in ya))
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif  ((x1 in ya) & ((x2 in ya) == false))
      INT2 = lt[find(lr.==neurons[x2])]
      for j = 1:length(bins)-1
        shank[x2,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x2)
      if isnan(c) == false
        push!(cor_store, c)
      end
    elseif ((x2 in ya) & ((x1 in ya) == false))
      INT2 = lt[find(lr.==neurons[x1])]
      for j = 1:length(bins)-1
        shank[x1,j] = length(find(bins[j] .<= INT2 .< bins[j+1]))
      end
      c = cor(vec(shank[x1,:]), vec(shank[x2,:]))
      push!(ya, x1)
      if isnan(c) == false
        push!(cor_store, c)
      end
    else
      #get count trains
      INT1 = lt[find(lr.==neurons[x1])]
      INT2 = lt[find(lr.==neurons[x2])]
      #add them to the counts matrix so you don't have to recalculate
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
  return mean(a)
end

function fano_train(Fcount, error_code)
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

function CV_ISI_ALLTIME(Neurons, t, r)
  CVS = []
  for i in eachindex(Neurons)
    INT = t[find(r.==Neurons[i])]
    isi = diff(INT)
    c = cv(isi)
    if isnan(c) == false
      push!(CVS, c)
    end
  end
  a = convert(Array{Float64}, CVS)
  mcv = emptiness(a, mean, -5)
  medcv = emptiness(a, median, -5)
  stdcv = emptiness(a, std, -5)
  return mcv, medcv, stdcv
end

function running_mean(signal, window, step)
  lank = div(length(signal), step)
  mod = div(window, step)
  smooth = zeros(lank-mod)
  for i in eachindex(smooth)
    ink = ((i-1)*step)+1
    smooth[i] = mean(signal[ink:ink+window])
  end
  return smooth
end

function tscore(mu, sigma)
  tsc = (mu-1.)/sigma
  return tsc
end

function SR(t, r, N)
    ST = []
    for i = 1:N
        f = find(r.==i)
        push!(ST, t[r])
    end
    return ST
end

runtime = 1000 #ms
fbinsize = 100 #ms
cbinsize = 50 #ms
h = 0.05
N = 2000
rt = runtime/1000. #convert to seconds
p = 0.1
S = 2.0
vth = 20. #sticking with 20 but be sure same in W and LIF_delta_solve

W = random_weights(N, p, vth)
CSC = sparse_rep(W,N)
@time t, r, input, volta = LIF_delta_solve_CSC(W, CSC, S, runtime, h)
println("Mean Rate: ", length(t)/(rt*N))
CVS = []
for i = 1:N
    m = find(r .== i)
    cvisi = cv(convert(Array{Float64}, diff(t[m])))
    if isnan(cvisi) == false
        push!(CVS, cvisi)
    end
end



# tic()
# Neurons = Neuron_finder(r, 10, 400)
# if ((E_Neurons == -5) | (I_Neurons == -5))
#   println("Too Little Spiking")
#   return "garbage"
# else
# Rates = [length(find(r .== i))/rt for i=1:N]
# countF = count_train_intron(fbinsize, t, r, Neurons, length(Neurons), false)
# spike_correlation = rand_pair_cor(cbinsize, t, r, Neurons, 1000)
# FANO_mean, FANO_median, FANO_std = fano_train(countF, -5)
# CV_mean, CV_median, CV_STD = CV_ISI_ALLTIME(Neurons, t, r)
# # Input_M = [mean(input[i,:]) for i=1:N]
# # Input_V = [var(input[i,:]) for i=1:N]
#
# println("FANO MEAN:",FANO_mean)
# println("CV MEAN:", E_CV_mean)
# println("COR MEAN:", spike_correlation)
# println("RATE MEAN: ", mean(Rates))
# # println("INPUT MEAN: ", mean(Input_M))
# # println("INPUT VAR: ", mean(Input_V))
# toc()
