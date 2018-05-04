## Tests for weights matrix

#2x2

# Wee = sum(W[1, 1:Ne_local])
# Wei = sum(W[1, Ne_local+1:N])
# Wie = sum(W[end, 1:Ne_local])
# Wii = sum(W[end, Ne_local+1:N])
#
# ks = sqrt(k)
#
# xee = Aee*ks
# xei = -Aei*ks
# xie = Aie*ks
# xii = -Aii*ks
#
# diffs = [Wee - xee, Wei - xei, Wie - xie, Wii - xii]
# for i in diffs
# println("Theoretical - observed = $(i)")
# end

Wee1 = sum(W[1, 1:Ne_local])
Wee2 = sum(W[1, 1:N_local])

Wee3 = sum(W[Ne_local + 1, 1:Ne_local])
Wee4 = sum(W[Ne_local + 1, Ne_local + 1:N_local])


Wei1 = sum(W[1, N_local + 1:N_local + Ni_local])
Wei2 = sum(W[1, N_local + Ni_local + 1:N])

Wei3 = sum(W[Ne_local + 1, N_local + 1:N_local + Ni_local])
Wei4 = sum(W[Ne_local + 1, N_local + Ni_local + 1:N])


Wie1 = sum(W[N_local + 1, 1:Ne_local])
WieL1 = sum(W[N_local + 1, Ne_local + 1:N_local])

WieL2 = sum(W[end, 1:Ne_local])
Wie2 = sum(W[end, Ne_local + 1:N_local])


Wii1 = sum(W[N_local + 1, N_local + 1:N_local + Ni_local])
Wii2 = sum(W[N_local + 1, N_local + N_local + Ni_local + 1:N])

Wii3 = sum(W[end, N_local + 1:N_local + Ni_local])
Wii4 = sum(W[end, N_local + Ni_local + 1:N])

ks = sqrt(k)

xee = Aee*ks
xei = -Aei*ks
xie = Aie*ks
xieL = Aie_NL*ks
xii = -Aii*ks

println("For the following checks, statements should return true, and numerical values should all be 0 to machine precision")
println("")

println(Wee1==Wee2)
println("theoretical ee - observed = ", Wee1 - xee)
println(Wee3==0)
println(Wee4 - Wee1)

println("")

println("theoretical ei - observed = ", Wei1 - xei)
println(Wei2==0)
println(Wei3==0)
println(Wei4 - Wei1)

println("")

println("theoretical ie - observed = ", Wie1 - xie)
println("theoretical ie_long - observed = ", WieL1 - xieL)
println(WieL2 - WieL1)
println(Wie2 - Wie1)

println("")

println("theoretical ii - observed = ", Wii1 - xii)
println(Wii2==0)
println(Wii3==0)
println(Wii4 - Wii1)
