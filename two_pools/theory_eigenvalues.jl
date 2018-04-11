function theory_pos_eigV(b, c, d)
  return b + sqrt(c+d)
end

#renormalized gain function values, which are W coefficients, set to 1 for now

function get_b(wee, wii)
  p = (wee .- wii .- 2)
  return p./2.
end

function get_c(wee, wei, wie, wii)
  p1 = (wii .+ wee) .^ 2.
  p2 = wei .* wie .* 4.
  return p1 .- p2
end

function get_d(wei, wieL)
  return 4. .* wei .* wieL
end

#renormalized gain function values, which are W coefficients, are free parameters

function get_b_g(wee, wii, ge, gi)
  p = ((ge*wee) .- (gi*wii) .- 2)
  return p./2.
end

function get_c_g(wee, wei, wie, wii, ge, gi)
  p1 = ((gi*wii) .+ (ge*wee)) .^ 2.
  p2 = wei .* wie .* 4. .* ge .* gi
  return p1 .- p2
end

function get_d_g(wei, wieL, ge, gi)
  return 4. .* wei .* wieL .* ge .* gi
end

function FI_G(S, vth, tau)
  denom = tau*log((S-(vth/tau))/S)
  return -1/denom
end

function firing_period(tau, vth, S)
  return -tau*log(1-(vth/(tau*S)))
end

function single_neuron(tau, vth, S, runtime, h)
  v0 = randn()*vth
  ntotal = Int64(div(runtime, h))
  v = v0
  m_leak = h/tau
  spikes = []
  for i = 1:ntotal
    v += (S*h) - (v * m_leak)
    if v > vth
      push!(spikes, i)
      v -= vth
    end
  end
  return spikes
end

function single_neuron_noise(tau, vth, S, runtime, h)
  v0 = randn()*vth
  ntotal = Int64(div(runtime, h))
  v = v0
  m_leak = h/tau
  spikes = []
  sh = sqrt(h)
  sigma = sqrt(S)
  ssh = sh*sigma
  eta = (randn(ntotal)*ssh) + S*h
  for i = 1:ntotal
    v = v - (v*m_leak) + eta[i]
    if v > vth
      push!(spikes, i)
      v -= vth
    end
  end
  return spikes
end

function single_neuron_noise_v(tau, vth, S, sigma, runtime, h)
  v0 = randn()*vth
  ntotal = Int64(div(runtime, h))
  v = v0
  m_leak = h/tau
  spikes = []
  sh = sqrt(h)
  ssh = sh*sigma
  eta = (randn(ntotal)*ssh) + S*h
  for i = 1:ntotal
    v = v - (v*m_leak) + eta[i]
    if v > vth
      push!(spikes, i)
      v -= vth
    end
  end
  return spikes
end

function estimated_gtiX(I)
  p = (I .-.6)*50.
  if p .< 0
    return 0
  else
    return p
  end
end

# runtime = 100000. # ms
# rts = runtime/1000. # s
# # SS = linspace(1., 2., 50)
# # ER = []
# # DR = []
# # NR = []
# # NRS = zeros(6, 100)
# # for i in SS
# #   nspikesd = (runtime/firing_period(20., 20., i))/rts
# #   push!(DR, nspikesd)
# # end
#
# SS2 = linspace(-.3, .4, 100)
# e1std = (2.13*(SS2 .^2.)) + (-1.7*SS2) + 1.83
# i1std = (-.064*(SS2 .^2.)) + (-1.34*SS2) + 1.51
#
# EZ = []
# IZ = []
#
# SD = [.5, 1., 1.5, 2., 2.5, 3.]
# for i in eachindex(SS2)
#   spikesnE = single_neuron_noise_v(20., 20., SS2[i], e1std[i], runtime, 0.1)
#   spikesnI = single_neuron_noise_v(20., 20., SS2[i], i1std[i], runtime, 0.1)
#   nspikesnE = length(spikesnE)/rts
#   nspikesnI = length(spikesnI)/rts
#   push!(EZ, nspikesnE)
#   push!(IZ, nspikesnI)
#
#   # spikesn = single_neuron_noise(20., 20., i, runtime, 0.1)
#   # nspikesn = length(spikesn)/rts
#   # spikes = single_neuron(20., 20., i, runtime, 0.1)
#   # nspikese = length(spikes)/rts
#   # push!(ER, nspikese)
#   # push!(NR, nspikesn)
#   # for j in eachindex(SD)
#   #   spikesn = single_neuron_noise_v(20., 20., SS2[i], SD[j], runtime, 0.1)
#   #   nspikesn = length(spikesn)/rts
#   #   NRS[j, i] = nspikesn
#   # end
# end
#
# SS2M = [estimated_gtiX(i) for i in SS2]
#
# # plot(SS2, ER, ".", label = "Simulation, Constant Input")
# # plot(SS, DR, ".", alpha = .5, label = "Analytical, Constant Input")
# # plot(SS2, NR, ".", label = "Simulation, Noisy")
# # plot(SS2, SS2M, ".", label = "Model g~'")
# # for i = 1:6
# #   plot(SS2, NRS[i,:][:], ".", label = "Sigma=$(SD[i])")
# # end
# plot(SS2, EZ, ".", label = "Excitatory")
# plot(SS2, IZ, ".", label = "Inhibitory")
# legend()
# xlabel("S (mV/ms)")
# ylabel("Firing rate (Hz)")
# title("Single Neuron FI Curve with Constant and Noisy Input")

N = 100

WEE_1 = ones(N)*.27
WEI_1 = ones(N)*2.2
WIE_1 = ones(N)*1.9
WIEL_1 = linspace(0, 10.5, N)
WII_1 = ones(N)*1.1
ge = .15
gi = -.1

b1 = get_b_g(WEE_1, WII_1, ge, gi)
c1 = get_c_g(WEE_1, WEI_1, WIE_1, WII_1, ge, gi)
d1 = get_d_g(WEI_1, WIEL_1, ge, gi)

eig1 = b1 + sqrt(complex(c1+d1))
eig2 = b1 + sqrt(complex(c1-d1))
eig3 = b1 - sqrt(complex(c1+d1))
eig4 = b1 - sqrt(complex(c1-d1))
figure(2)
# plot(WIEL_1, b1, label = "b")
# plot(WIEL_1, c1, label = "c")
# plot(WIEL_1, d1, label = "d")
plot(WIEL_1, eig1, label = "b + sqrt(c+d)")
plot(WIEL_1, eig2, label = "b + sqrt(c-d)")
plot(WIEL_1, eig3, label = "b - sqrt(c+d)")
plot(WIEL_1, eig4, label = "b - sqrt(c-d)")

# plot(WIEL_1, eig2, label = "b + sqrt(c-d)")
# plot(WIEL_1, eig3, label = "b - sqrt(c+d)")
# plot(WIEL_1, eig4, label = "b - sqrt(c-d)")
legend()
xlabel("WIE LONG")
title("Theoretical Eigenvalues, Components")
# theory_eig1 = theory_pos_eigV(b1, c1, d1)


#
# yd = zeros(301)
# xd = zeros(301)
#
# for i = 100:400
#   yd[i-99] = ER[end] - ER[i]
#   xd[i-99] = SS2[end] - SS2[i]
# end
