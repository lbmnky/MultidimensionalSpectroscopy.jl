using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles, Random, Combinatorics, Distributions

# make sure to set script directory as pwd()
cd(@__DIR__)

pygui(true)

calc_2d = false

cmp = create_colormap("bright");


# https://qudev.phys.ethz.ch/static/content/science/papers/Fink2008.pdf

ωr = 4
Na = 2               # number of cavity fock states
A  = FockBasis(Na)
a  = destroy(A)
at = create(A)


ωa = 4.1
Nσ = 3              # number of levels of matter
Σ  = NLevelBasis(Nσ)
sm  = transition(Σ, 1, 2)
sp = transition(Σ, 2, 1)

sz = sp * sm

H_r = (ωr * at * a)
H_a = (ωa * sz)


"""
generate H for Nm TLSs
H_am = H_a ⊗ one(H_a) ⊗ one(H_a) ⊗ ... + one(H_a) ⊗ H_a ⊗ one(H_a) ⊗ ... +
        one(H_a) ⊗ one(H_a) ⊗ H_a ⊗ ... + ...
"""
Nm = 2

# vector with molecule indices
R = collect(1:Nm)

for i in R
    if i == 1
        H_am = H_a
        smm = sm
    else
        global H_am = H_am ⊗ one(H_a)
        global smm = smm ⊗ one(H_a)
    end
end

spp = dagger(smm)



disorder = zeros(Nm)

for j in R
    println(circshift(R,j-1))
    if j == 1
        disorder[j] = 1
        global H_amb = H_am
    else
        disorder[j] = 1 + (rand(Float64, 1)[1] - 0.5) / 20
        global H_amb = H_amb + disorder[j] * permutesystems(H_am,circshift(R,j-1))
    end
end

for i in R
    if i == 1
        println(circshift(R,i-1))
        global H_x = at ⊗ one(H_amb) * (one(A) ⊗ smm) + (one(A) ⊗ spp) * (a ⊗ one(H_amb))

    else
        global H_x = H_x + at ⊗ one(H_amb) * (one(A) ⊗ permutesystems(smm,circshift(R,i-1))) +
                    one(A) ⊗ permutesystems(spp,circshift(R,i-1)) * (a ⊗ one(H_amb))
    end
end

#H_x = at ⊗ one(H_am) * (one(A) ⊗ sm1) + one(A) ⊗ sp1 * (a ⊗ one(H_am)) +
#        at ⊗ one(H_am) * (one(A) ⊗ sm2) + one(A) ⊗ sp2 * (a ⊗ one(H_am)) +
#          at ⊗ one(H_am) * (one(A) ⊗ sm3) + one(A) ⊗ sp3 * (a ⊗ one(H_am))

#H_x = at ⊗ one(Σ) * (one(A) ⊗ sm) + one(A) ⊗ sp * (a ⊗ one(Σ))

H_TC = H_r ⊗ one(H_amb) + (one(A) ⊗ H_amb)  + 0.4 * H_x

figure(figsize=(5,4));
imshow(real(H_TC.data-H_x.data + 20*H_x.data),origin="lower",cmap="Blues")
title("H_TC with off-diagonal elements * 20")
#pcolormesh(real(H_x.data))

figure(figsize=(5,4));
imshow(real(H_TC.data),cmap="Greys")

title("H_TC")
colorbar()


num_of_exc = 0
ground = tensor(fockstate(A,0),repeat([nlevelstate(Σ,1)],Nm-num_of_exc)...)

E = eigvals(dense(H_TC).data)
xn = 0
figure(figsize=(5,4));
title("Energy levels (Eigenvalues) of H_TC")
hlines([0, ωr], -.65, -.35)
hlines([0],  .35,  .65)
hlines([ωa .* disorder[:]],  .35,  .65)
for i in 1:length(E)
    println(i)
    if i > 1 && E[i] - E[i-1] < 0.1
        hlines(E[i],-.15,.15)
        global xn += .02
    else
        hlines(E[i],-.15,.15)
        global xn = 0
    end
end
xlim([-.8, .8])

for i in R
    if i == 1
        Psi0 = fockstate(A,0) ⊗ nlevelstate(Σ,1)
    else
        global Psi0 =  Psi0 ⊗ nlevelstate(Σ,1)
    end
end

# GS density matrix
rho0 = dm(Psi0)

# excite density matrix
rho1 = (at ⊗ one(smm)) * rho0 * dagger((at ⊗ one(smm)))

# debug: view rho0 after "excitation"
# imshow(real(rho0.data),origin="lower",cmap="Blues")

tlist = [0:0.1:100;]

# [relaxation within matter, relaxtion with light field]
J = [one(at) ⊗ (smm), a ⊗ one(smm)]
Γ = [.6, .06]

# master equation time evolution
tout, rhot = timeevolution.master(tlist,rho1,H_TC,J;rates=Γ)

# visualize evolution of density matrix
# cmds.view_dm_evo(rhot,5)

# calculate population in light and matter
nc = real(expect((at ⊗ one(smm)) * (a ⊗ one(smm)), rhot))
na = real(expect(one(at) ⊗ (spp * smm), rhot))

# plot population
figure(figsize=(10,4));
plot(tout,nc)
plot(tout,na)

# calculate correlation function
corr = timecorrelations.correlation(tlist, rho1, H_TC, J, (at ⊗ one(smm)),
            (a ⊗ one(smm)); rates=Γ)

# calculate spectrum from correlation
ω, spec = timecorrelations.correlation2spectrum(tlist, corr; normalize_spec=true)

# plot both
figure(figsize=(10,4));
ax1 = subplot(211);
plot(tlist,corr)
ax2 = subplot(212);
plot(ω, spec)

# parameters for 2D
F = 0.6 * J
H = H_TC

# interaction light with cavity and TLSs
#μ12 = at ⊗ (sm*sp) + a ⊗ (sm*sp) + (a*at) ⊗ sp + (a*at) ⊗ sm
# interact light only with cavity
μ12 = at ⊗ (smm*spp) + a ⊗ (smm*spp)
# interact with both (?) for higher states
#μ23 = at ⊗ (spp*smm) + a ⊗ (spp*smm) + (at*a) ⊗ spp + (at*a) ⊗ smm
# or only the atom ?
#μ23 = (at*a) ⊗ spp + (at*a) ⊗ smm
# or only the atom ?
μ23 = at ⊗ (spp*smm) + a ⊗ (spp*smm)

# GS density matrix (reset it you idiot)
rho0 = dm(Psi0)

if calc_2d
    ## calculate (complex) 3rd order corr function (with T=0)
    zp = 11 # zeropad up to 2^zp

    ## calculate 2D spectra at
    T = [0] #fs

    spectra2d = Array{out2d}(undef, length(T))
    for i = 1:length(T)
        spectra2d[i] = make2Dspectra(tlist,[rho0, rho0],[H, H],[F, F],[μ12, μ12],[μ23, μ23],T[i],
                                            "lindblad";debug=true,zp=zp)
    end

    ## crop 2D data and increase dw
    spectra2d = [crop2d(spectra2d[i],1;w_max=20,step=2) for i = 1:length(T)]

    ## assign ω-axis from output
    ω = spectra2d[1].ω
end

if calc_2d
    ## plot 2D spectra for each(?) T
    # what to plot
    rep="absorptive"
    scal="lin"

    ## make  subplot layout
    nplots = length(T); ncols = ceil(sqrt(nplots)); nrows = ceil(nplots / ncols)

    ## create figure
    figure(figsize=(ncols*3,nrows*3));
    suptitle(rep * " 2D spectrum (" * scal * ". scaling)")
    for i = 1:nplots
        subplot(nrows,ncols,i,aspect="equal")
        plot2d(ω,round.(spectra2d[i].full2d,digits=1);repr=rep,scaling=scal)
        title("2D spectrum at $(T[i]) fs")
        #xlim([0, E₁+E₂]); ylim([0, E₁+E₂]);
    end
end

# The end
#########################################################################

#H_TC = H_r ⊗ one(Σ) ⊗ one(Σ) + one(A) ⊗ H_a ⊗ one(Σ)

"""
b = NLevelBasis(2)
H = diagonaloperator(b, [0, 1])
j12 = transition(b, 1, 2) # decay from 2nd into 1st level

Nb = 3
b_mb = ManyBodyBasis(b, bosonstates(b, [1, 1, 1]))
H_mb = manybodyoperator(b_mb, H)
j12_mb = manybodyoperator(b_mb, j12)


wc = 1
N = 1               # number of cavity fock states
b_fock = FockBasis(N)
a = destroy(b_fock) ⊗ one(b_mb)
at = dagger(a)

wa = 1
sm = one(b_fock) ⊗ j12_mb
sp = dagger(sm)
println(dense(sm))

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