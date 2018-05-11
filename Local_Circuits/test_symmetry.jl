include("2x2_measure_synapse.jl")

Aee = 12.5#parse(Float64, ARGS[1]);
Aei = 50.#parse(Float64, ARGS[2]);
Aie = 20.#parse(Float64, ARGS[3]);
Aii = 50.#parse(Float64, ARGS[5]);
k = 500;

N = 5000;
IFRAC = 2;

Ni = Int64(round(N/IFRAC));
Ne = N - Ni;
W = local_random_2x2_symmetric(N, IFRAC, k, Aee, Aei, Aie, Aii);

k2 = round(Int64, k/2)
Ne2 = round(Int64, Ne/2)
Ni2 = round(Int64, Ni/2)

plot(W[1, 1:Ne2][:], W[1, Ne2 + 1:Ne][:])
plot(W[1, 1:Ne2][:], W[Ne2 + 1, 1:Ne2][:])
plot(W[1, 1:Ne2][:], W[Ne2 + 1, Ne2 + 1:Ne][:])

plot(W[1, 1 + Ne:Ne2 + Ne][:], W[1, Ne2 + 1 + Ne:Ne + Ne][:])
plot(W[1, 1 + Ne:Ne2 + Ne][:], W[Ne2 + 1, 1 + Ne:Ne2 + Ne][:])
plot(W[1, 1 + Ne:Ne2 + Ne][:], W[Ne2 + 1, Ne2 + 1 + Ne:Ne + Ne][:])

plot(W[Ne + 1, 1:Ne2][:], W[Ne + 1, Ne2 + 1:Ne][:])
plot(W[Ne + 1, 1:Ne2][:], W[Ne + Ne2 + 1, 1:Ne2][:])
plot(W[Ne + 1, 1:Ne2][:], W[Ne + Ne2 + 1, Ne2 + 1:Ne][:])

plot(W[Ne + 1, 1 + Ne:Ne2 + Ne][:], W[Ne + 1, Ne2 + 1 + Ne:Ne + Ne][:])
plot(W[Ne + 1, 1 + Ne:Ne2 + Ne][:], W[Ne + Ne2 + 1, 1 + Ne:Ne2 + Ne][:])
plot(W[Ne + 1, 1 + Ne:Ne2 + Ne][:], W[Ne + Ne2 + 1, Ne2 + 1 + Ne:Ne + Ne][:])
