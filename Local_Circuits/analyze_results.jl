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
    c = cv(convert(Array{Float64}, isi))
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
