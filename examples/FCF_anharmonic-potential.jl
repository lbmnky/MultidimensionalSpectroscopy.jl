#!/usr/bin/julia
using PyPlot, QuantumOptics, LinearAlgebra, Random, Combinatorics, Distributions

# https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/13%3A_Coupling_of_Electronic_and_Nuclear_Motion/13.01%3A_The_Displaced_Harmonic_Oscillator_Model
# https://www.scm.com/doc/ADF/Input/Vibrationally_resolved_el_spectra.html

# make sure to set script directory as pwd()
cd(@__DIR__)

pygui(true)

b = PositionBasis(0, 20, 200)
q0 = 7
D0 = 45.5
β0 = .035
μ0 = 4
E0 = 0

S = 1 # Huang-Rhys
α = 1 # M * ω / ħ with M reduced mass
ΔQ = sqrt(2*S/α) # ΔQ displacement
#q1 = 8.5
q1 = q0 + ΔQ
D1 = 45.5#45.5
β1 = .0375
μ1 = 4
E1 = 2.5

N = 10 # number of considered states

a = q1 - q0 # displace excited state

ω0 = sqrt(2 * D0 * β0^2 / μ0)
ω1 = sqrt(2 * D1 * β1^2 / μ1)

q = position(b)
p = momentum(b)

Hg =
    p^2 / (2 * μ0) +
    D0 * (one(b) - exp(-β0 * (dense(q) - one(b) * q0)))^2 +
    one(b) * E0

# Harmonic 
harm = false
if harm == true
    Hg = 
        p^2 / (2 * μ0) +
        .5 * μ0 * ω0^2 * (dense(q) - one(b) * q0)^2
end

energies_g, states_g = eigenstates((Hg + dagger(Hg)) / 2, N)

He =
    p^2 / (2 * μ1) +
    D1 * (one(b) - exp(-β1 * (dense(q) - one(b) * q1)))^2 +
    one(b) * E1

if harm == true
    He = 
         p^2 / (2 * μ1) +
        .5 * μ1 * ω0^2 * (dense(q) - one(b) * q1)^2 +
        one(b) * E1
end

# Alternatively to new q1 can displace with a !
#He = exp(-1im * a * p) * He * exp(1im * a * p)

energies_e, states_e = eigenstates((He + dagger(He)) / 2, N)

figure(figsize = (8, 4))
subplot(121)
xpoints = samplepoints(b)


for i = 1:length(states_g)
    plot(
        xpoints,
        abs2.(states_g[i].data) .* 20 .+ energies_g[i],
        #-real.(states_g[i].data) .+ energies_g[i],
        color = "k",
        linewidth = 1,
    )
    plot(
        xpoints,
        abs2.(states_e[i].data) .* 20 .+ energies_e[i],
        color = "g",
        linewidth = 1,
    )
end

if harm != true
    ff = D0 * (1 .- exp.(-β0 * (xpoints .- q0))) .^ 2 .+ E0
    fe = D1 * (1 .- exp.(-β1 * (xpoints .- q0 .- a))) .^ 2 .+ E1
else
    ff = .5 * μ0 * ω0^2 .* (xpoints .- q0).^2 .+ E0
    fe = .5 * μ1 * ω1^2 .* (xpoints .- q1).^2 .+ E1
end

#plot(xpoints, D0 * (1 .- exp.(-β0 * (xpoints .- q0))) .^ 2 .+ E0, alpha = 0.5)
plot(xpoints, ff, alpha = 0.5)

fill_between(
    xpoints,
    #D0 * (1 .- exp.(-β0 * (xpoints .- q0))) .^ 2 .+ E0,
    ff,
    alpha = 0.5,
)
plot(
    xpoints,
    #D1 * (1 .- exp.(-β1 * (xpoints .- q0 .- a))) .^ 2 .+ E1,
    fe,
    alpha = 0.5,
)

ylim([0, energies_e[end] + 2]);
xlim([-0, xpoints[end]]);

xlabel("Position")
ylabel("Energy")

FC = complex(zeros(1, N))

subplot(122)
# |⟨v|i⟩|² , where |v⟩ is the vibrational level on the excited electronic state and |i⟩ the initial state
for i = 1:1
    for j = 1:N
        FC[i, j] = abs2(states_g[i]' * states_e[j])^2
        #bar(
        #    energies_e[j] .- energies_g[1],
        #    real(FC[1, j]),
        #    color = "b",
        #    width = 0.1,
        #)
    end
end
γ = 0.025 # γ has an effect on the intensity envelop
ω = collect(0:0.01:10)
abs_crossSection = real(ω .* sum([γ * FC[1,i]./((i * ω0 + E1 - ω0 .- ω).^2 .- (im * γ)^2) for i in 1:length(energies_e)]))
#plot(ω,abs_crossSection./maximum(abs_crossSection)./2)
plot(ω,abs_crossSection)

FR = complex(zeros(1, N))

# ⟨f|v⟩⟨v|i⟩ ... for Raman cross-section
for i = 1:1
    for j = 1:N
        FR[i, j] = abs2(states_g[i+1]' * states_e[j] * states_e[j]' * states_g[i])
        #bar(
        #    energies_e[j] .- energies_g[1],
        #    real(FC[1, j]),
        #    color = "r",
        #    width = 0.05,
        #)
    end
end

γ = 0.025 # γ has an effect on the intensity envelop
ω = collect(0:0.01:10)
abs_crossSection = real(ω .* sum([γ * FR[1,i]./((i * ω0 + E1 - ω0 .- ω).^2 .- (im * γ)^2) for i in 1:length(energies_e)]))
#plot(ω,abs_crossSection./maximum(abs_crossSection)./2)
plot(ω,abs_crossSection)

xlabel("Energy")
#ylabel("|<ϕf|ϕi>|²")
ylabel("Intensity")

legend(["|⟨v|i⟩|²", "⟨f|v⟩⟨v|i⟩"])
xlim(2,4)
tight_layout()

# TODO convert to FockBasis and do 2D !
b_F = FockBasis(N)
test = transform(b, b_F, x0 = q0)
#...
A = test' * Hg * test
states   = eigvecs(dense(A).data)
energies = eigvals(dense(A).data)
