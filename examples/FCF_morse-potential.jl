#!/usr/bin/julia
# using PyPlot
# using QuantumOptics
# https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/13%3A_Coupling_of_Electronic_and_Nuclear_Motion/13.01%3A_The_Displaced_Harmonic_Oscillator_Model
# https://www.scm.com/doc/ADF/Input/Vibrationally_resolved_el_spectra.html

# make sure to set script directory as pwd()
cd(@__DIR__)


pygui(true)

b = PositionBasis(0,10,200)
q0 = 2.5
D0 = 45.5
β0 = .1475
μ0 = 2
E0 = 0

q1 = 3.5
D1 = 45.5
β1 = .375
μ1 = 2
E1 = 25

N = 5 # number of considered states

a = q1 - q0 # displace excited state


ω0 = sqrt(2*D0*β0^2/μ0)
ω0 = sqrt(2*D1*β1^2/μ1)


q = position(b)
p = momentum(b)

Hg = p^2/(2*μ0) + D0*(one(b) - exp(-β0 * (dense(q)-one(b)*q0) ) )^2 + one(b) * E0

energies_g, states_g = eigenstates((Hg+dagger(Hg))/2, N)

He = p^2/(2*μ1) + D1*(one(b) - exp(-β1 * (dense(q)-one(b)*q1) ) )^2 + one(b) * E1
# Alternatively to new q1 can displace with a !
#He = exp(-1im * a * p) * He * exp(1im * a * p)

energies_e, states_e = eigenstates((He+dagger(He))/2, N)

figure(figsize=(8,4))
subplot(121)
xpoints = samplepoints(b)

for i=1:length(states_g)
    plot(xpoints, abs2.(states_g[i].data).*40 .+ energies_g[i],color="k",linewidth=1)
    plot(xpoints, abs2.(states_e[i].data).*40 .+ energies_e[i],color="g",linewidth=1)
end

plot(xpoints, D0 * (1 .- exp.(-β0 * (xpoints .- q0))).^2 .+ E0,alpha=0.5)
fill_between(xpoints, D0 * (1 .- exp.(-β0 * (xpoints .- q0))).^2 .+ E0, alpha=0.5)
plot(xpoints, D1 * (1 .- exp.(-β1 * (xpoints .- q0 .- a))).^2 .+ E1,alpha=0.5)

ylim([0, energies_e[end]+5]); xlim([0,xpoints[end]])

xlabel("Position")
ylabel("Energy")

FC = complex(zeros(1,N))

subplot(122)
for i = 1:1
    for j = 1:N
        FC[i,j] = abs2(states_g[i]' * states_e[j])
        bar(energies_e[j] .- energies_g[1], real(FC[1,j]),color="b",width=0.2)
    end
end

xlabel("Energy")
ylabel("|<ϕf|ϕi>|²")
tight_layout()
show()
