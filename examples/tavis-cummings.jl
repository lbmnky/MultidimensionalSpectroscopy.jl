#!/usr/bin/julia
# using PyPlot
# using QuantumOptics

# make sure to set script directory as pwd()
cd(@__DIR__)

# include my custom cmds module
if Sys.iswindows()
    include("..\\cmds.jl")
    fn = "01_Output"
else
    include("../cmds.jl")
    fn = "01_Output"
end

pygui(true)

calc_2d = true

cmp = cmds.create_colormap("bright");




b = NLevelBasis(2)
H = diagonaloperator(b, [0, 1])
j12 = transition(b, 1, 2) # decay from 2nd into 1st level

Nb = 2
b_mb = ManyBodyBasis(b, bosonstates(b, Nb))

H_mb = manybodyoperator(b_mb, H)

j12_mb = manybodyoperator(b_mb, j12)


wc = 0.9
N = 1               # number of cavity fock states
b_fock = FockBasis(N)
a = destroy(b_fock) ⊗ one(b_mb)
at = dagger(a)

wa = 1.1
sm = one(b_fock) ⊗ j12_mb
sp = dagger(sm)
println(sm)

Psi0 = fockstate(b_fock,0) ⊗ basisstate(b_mb,[Nb,0])    # start with an excited cavity
rho0 = dm(Psi0)

rho0 = rho0 * a
g = 0.24

### this should be on the cloud

H_full = wc * at * a + wa * sp * sm + g * (at * sm + a * sp)

Γ = [0.01]
J_mb = [sm]

tlist = [0:0.1:50;]

tout, rhot = timeevolution.master(tlist,rho0,H_full,J_mb;rates=Γ)

nc = real(expect(at * a, rhot))
na = real(expect(sp * sm, rhot))

figure(figsize=(6,3))
plot(tout,nc)
plot(tout,na)
