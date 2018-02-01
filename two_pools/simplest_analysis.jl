function von_mises_dist(x, k, mu, N)
  a = exp(k*cos(x-mu))/(N*besseli(0, k))
end
function weights(Ne,Ni,kee,kei,kie_L,kii,Aee,Aei,Aie_L,Aii,pee,pei,pie,pii)
# Construct weight functions

wee = zeros(Ne,Ne)
wei = zeros(Ne,Ni)
wie = zeros(Ni,Ne)
wii = zeros(Ni,Ni)

we = 2*collect(0:Ne-1)*pi/Ne
wi = 2*collect(0:Ni-1)*pi/Ni

for i = 1:Ne
    wee[i,:]=Aee*circshift(von_mises_dist(we, kee, 0, Ne), i-1)
    wei[i,:]=Aei*circshift(von_mises_dist(wi, kei, 0, Ni),div(Ni*i,Ne)-1)
end
for i = 1:Ni
    wie[i,:]=Aie_L*circshift(von_mises_dist(we, kie_L, 0, Ne),div(Ne*i,Ni)-1)# + Aie_D*circshift(von_mises_dist(we, kie_D, pi, Ne),div(Ne*i,Ni)-1)
    wii[i,:]=Aii*circshift(von_mises_dist(wi, kii, 0, Ni),i-1)
end

# keep connections with probability pab
wee = wee.*(rand(Ne,Ne).<pee)
wei = wei.*(rand(Ne,Ni).<pei)
wie = wie.*(rand(Ni,Ne).<pie)
wii = wii.*(rand(Ni,Ni).<pii)
#println("weights complete")

return wee, -wei, wie, -wii
end

function lif(h,total,Ne,Ni,wee,wei,wie,wii,Se1,Se2)
# Leaky integrate-and-fire network

vth = 20
ntotal = round(Int,total/h)

ve = rand(Ne)*vth
vi = rand(Ni)*vth

se = zeros(Int,Ne)
si = zeros(Int,Ni)

FPe = div(Ne,5)
Se = zeros(Ne)

P1s = round(Int,Ne/4-FPe)
P1e = round(Int,Ne/4+FPe)

P2s = round(Int,3*Ne/4-FPe)
P2e = round(Int,3*Ne/4+FPe)

Se[P1s:P1e] = Se1
Se[P2s:P2e] = Se2
half = div(Ne, 2)

tau = 20
time_e = Float64[0]
raster_e = Float64[0]
time_i = Float64[0]
raster_i =Float64[0]

leak = exp(-h/tau)
kill_flag = false
for  iter = 1:ntotal
# time loop

     # administer leak
     ve *= leak
     vi *= leak

     ve += tau*Se*(1.-leak)
     # reset and record spike times
     se = (ve.>vth)
     si = (vi.>vth)

     ve -= vth*se
     vi -= vth*si

     nspikese = sum(se)
     nspikesi = sum(si)

     if nspikese > 0
     	spe = find(se)
     	for j = 1:nspikese
	    # Give synaptic kicks
  	    ve += wee[:,spe[j]]
       	vi += wie[:,spe[j]]
  	    push!(raster_e,spe[j])
  	    push!(time_e,iter)
      end
     end
     if nspikesi > 0
     	spi = find(si)
     	for j = 1:nspikesi
	    # Give synaptic kicks
  	    ve += wei[:,spi[j]]
       	vi += wii[:,spi[j]]
  	    push!(raster_i,spi[j])
  	    push!(time_i,iter)
      end
     end
     if iter == 1000
       if (length(raster_e) > 200*Ne*h*iter*.5*(1/1000)) | (length(raster_i) > 200*Ne*h*iter*.5*(1/1000))
         kill_flag = true
         break
       end
     end

end # of time loop

return time_e,raster_e,time_i,raster_i, kill_flag
end


function write_result(wta, etm, etv, ebm, ebv, cor_within, cor_between, ntf, mif, mcv, lte, id)
  println("##RESULT $(wta), $(etm), $(etv), $(ebm), $(ebv), $(cor_within), $(cor_between), $(ntf), $(mif), $(mcv), $(lte), $(id)")
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

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  ###x = dv/dt
  ###vth = v2 - xt
  ###(vth - v2)/x = t
  return t
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
  else
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
end
