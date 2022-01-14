using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles

## some theory
# chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/13%3A_Coupling_of_Electronic_and_Nuclear_Motion/13.01%3A_The_Displaced_Harmonic_Oscillator_Model
# www.scm.com/doc/ADF/Input/Vibrationally_resolved_el_spectra.html

# make sure to set script directory as pwd()
cd(@__DIR__)

pygui(true)

cmp = create_colormap("bright");

calc_2d = true

## Huang-Rhys factor
S = .3
α = 1 # M * ω / ħ with M reduced mass
ΔQ = sqrt(2*S/α) # ΔQ displacement
d = ΔQ# sqrt(D_HR)

#S = M * ω / ħ * ΔQ^2
# bond length change between C-C and C=C is  0.19 Å

b_tls = NLevelBasis(2)          # Hilbert-space of system                   Basis: {|ui⟩}
b_vib = FockBasis(5)            # Hilbert-space of oscillator               Basis: {|vi⟩}
b     = b_tls ⊗ b_vib          # combined basis

j21 = transition(b_tls,2,1)                     # |e⟩⟨g|
j12 = dagger(j21)                               # |g⟩⟨e|
at  = create(b_vib)                             # ...
a   = dagger(at)
D   = displace(b_vib,d)

# ground states
Eg = 0
ωg = .15
Hg_el = Eg * j12 * j21     # == HG = g * Eg ⊗ dagger(g)
#Hg_vib = ωg * (at*a + a*at * 1/2)
Hg_vib = ωg * (number(b_vib) + identityoperator(b_vib) * 1/2)
Hg = Hg_el ⊗ one(b_vib) + j12*j21 ⊗ (Hg_vib)
#Hg = one(Hg_el) ⊗ Hg_vib
#Hg = Hg_el ⊗ Hg_vib

# HAVE TO created excited states w/o displaced potential and then use displacement
# operator to act on transition dipole moment ... or calculate FC factors and use
# this as transition dipole operator ... #CHECK

# excited state
Ee = 1.65
ωe = 0.15
He_el = Ee * j21 * j12
#He_vib = dagger(D) * Hg_vib * D # issue with D and dagger(D)
#He_vib = D * Hg_vib * dagger(D) # issue with D and dagger(D)
#He_vib = Hg_vib  # TODO figure out what's going on here
He_vib  = ωe * (number(b_vib) + identityoperator(b_vib) * 1/2)
He      = He_el ⊗ one(b_vib) + j21*j12 ⊗ (He_vib)
#He = He_el ⊗ one(b_vib) + j21*j12 ⊗ He_vib
#He = He_el ⊗ He_vib

## full system Hamiltonian
#H = Hg + (one(b_tls) ⊗ D) * He * (one(b_tls) ⊗ D')
H = Hg + He
#H = (one(b_tls) ⊗  D ) * H * (one(b_tls) ⊗  D )'
## diagonalized Hamiltonian
states   = eigvecs(dense(H).data)
energies = eigvals(dense(H).data)

# and GS
energies_g, states_g = eigenstates(dense(Hg_vib))
# and ES
energies_e, states_e = eigenstates(dense(He_vib))

## plot energy levels
N = length(energies)
println(Int8(ceil(N/2)))
figure(figsize=(6,4))
subplot(121)
scatter(zeros(Int8(ceil(N/2))), energies[1:Int8(ceil(N/2))])
x = [-1:0.1:2*d;]

# TODO: transform to position basis to plot PECs
#xx = PositionBasis(0,10,200)
#xxx = position(xx)

m = 50                        # arbitrary value for m
f = m * ωg^2
# GS
plot(x,f*x.^2)
# ES
plot(x,f * (x.-d).^2 .+ Ee)
scatter(zeros(Int8(floor(N/2))) .+ d, energies[end-Int8(floor(N/2))+1:end])
ylim(0, 6)

#for i=1:length(states_g)
#    plot(xpoints, abs2.(states_g[i].data).*40 .+ energies_g[i],color="k",linewidth=1)
#    plot(xpoints, abs2.(states_e[i].data).*40 .+ energies_e[i],color="g",linewidth=1)
#end

## calculate Franck-Condon factors
FC = complex(zeros(length(energies_g),length(energies_e)))
for ii = 1:length(energies_g)
    for jj = 1:length(energies_e)
        FC[ii,jj] = abs(conj(states_g[1]' * (D*states_e[jj])) * states_g[1]' * (D*states_e[ii]))
        #FC[ii,jj] = abs(conj(states_g[1]' * states_e[jj]) * (states_g[1]' * states_e[ii]))
    end
end

# the previous algorith corresponds to calculating the FC-factors
l = complex(zeros(length(energies_g),length(energies_e)))
for ii = 1:length(energies_g)
    for jj = 1:length(energies_e)
        l[ii,jj] = exp(-d^2) * d^(2*jj) / factorial(jj)
    end
end

ln = zeros(5,)
for n = 1:5
    ln[n] = exp(-d^2) * d^(2*n) / factorial(n)
end

μ = j12  ⊗ D + j21 ⊗ D'
#μ12 = μ ⊗ (D * D')
#μ12 = μ ⊗ one(b_vib)# + (j21) ⊗  (D')# + j21 ⊗ D
#μ12 = (one(b_tls) ⊗ D') * μ
#μ12 = μ ⊗ one(b_vib)
μ12 = μ
μ23 = copy(μ)
fill!(μ23.data,0)

# initialize in groud state
Psi0 = nlevelstate(b_tls,1) ⊗ fockstate(b_vib,0)
rho0 = Psi0 ⊗ dagger(Psi0)  # alt.: rho0 = dm(Psi0)

Psi1 = nlevelstate(b_tls,2) ⊗ fockstate(b_vib,0)
rho1 = dm(Psi1)

## TODO used thermally populated ground state instead
#T = 0.01                               # Temperature T in what unit ?
#rho0 = thermalstate(H,T)

#c_ops = [sqrt(0.05) * one(b_tls) ⊗ (a+at), sqrt(0.075) * (j21+j12) ⊗ one(b_vib)]
Γ = [sqrt(0.01), sqrt(0.001), sqrt(0.04), sqrt(0.04)]
L = Γ .* [one(b_tls) ⊗ a, j12 ⊗ one(b_vib), one(j21 * j12) ⊗ (at * a), (j21 * j12) ⊗ one(at * a)]

tlist = [0:1:150;]

corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)

#tout, rhot = timeevolution.master(tlist,μ12*rho0,H,L)
#corr = expect(μ12,rhot)

subplot(222)
plot(tlist,real.(corr))

zp = 12
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

ω_abs, spec_abs = timecorrelations.correlation2spectrum(tnew, corr;
                                                        normalize_spec=true)


corr = []

corr = timecorrelations.correlation(tlist, rho1, H, L, μ12, μ12)

tout, rhot = timeevolution.master(tlist,μ12 * rho1,H,L)
corr = expect(μ12,rhot)

zp = 11
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

ω_em, spec_em = timecorrelations.correlation2spectrum(tnew, corr;
                                                        normalize_spec=true)

subplot(224)
plot(ω_abs,spec_abs)
plot(-ω_em ,spec_em,"r:")

legend(["absorption", "emission"])

γ = 0.02 # γ has an effect on the intensity envelop
abs_crossSection = real(ω_abs .* sum([γ * FC[1,i]./((i * ωe + Ee - ωg .- ω_abs).^2 .- (im * γ)^2) for i in 1:length(energies_e)]))
#plot(-ω_abs,abs_crossSection./maximum(abs_crossSection))

"""
temp = FC
temp = temp ./ maximum(real(temp))
for jj = 1:length(energies_g)
    plot(-([Ee-Eg, Ee-Eg] .+ 2*(jj-1)*ωg), [0, real(temp[jj,jj])], color="r",linewidth=1)
end
"""

tight_layout()

tout, rhot = timeevolution.master(tlist,μ12 * rho0 * μ12,H,L)


vecs   = eigvecs(dense(H).data)
eivecs = [Ket(H.basis_l,vecs[:,i]) for i in 1:length(vecs[1,:])]

figure()
title("eigen/coupled basis") # levels: 1 -> GS, 2 -> LP, 3 - N-1 -> ..., N -> UP
for i in 1:length(eivecs)
    es_n = real(expect(dm(eivecs[i]),rhot))
    plot(tout,es_n,linestyle="dashed", label=string(i)); #append!(leg,[string(i)])
    ylim((-0.1, 1))
end
legend()

F = L # works with method "lindblad"

if calc_2d

    ## zeropad up to 10^zp
    zp = 11

    ## calculate 2d spectra at
    T = [0, 100]
    #T = collect([0:5:100;])

    spectra2d = Array{out2d}(undef, length(T))

    # use multithreding...
    Threads.@threads for i = 1:length(T)
    # ... or not
    #for i = 1:length(T)
        spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                            "lindblad";debug=true,use_sub=false,
                                                zp=zp, t2coh="kin");
    end

    ## crop 2D data and increase dw
    spectra2d = [crop2d(spectra2d[i],1.4;w_max=2.2,step=1) for i = 1:length(T)]
end

if calc_2d
    ## simulate effect of laser spectrum
    laserSpec =  .3 * exp.(-(spectra2d[1].ω.-1.35).^2/(2*(.19)^2)) +
                  exp.(-(spectra2d[1].ω.-1.6).^2/(2*(.11)^2)) +
                   .3 * exp.(-(spectra2d[1].ω.-1.85).^2/(2*(.05)^2))

    ## plot laser spectrum
    #figure()
    #plot(out2d[1].ω,laserSpec)

    ## 1D convoltion (in ω_excitation)
    IRF = (laserSpec.^2*(ones(length(laserSpec)))')
    ## 2D convolution
    IRF = (laserSpec.^2*(laserSpec.^2)')

    #for i in 1:length(out2d)
        conv2d = spectra2d[1].full2d .* IRF
    #end

    ## plot 2D spectra for each(?) T
    # what to plot
    rep  = "absorptive"
    scal = "lin"

    ## make  subplot layout
    nplots = length(T);
    ncols  = Int32(ceil(sqrt(nplots)));
    nrows  = Int32(ceil(nplots / ncols));

    # determine maximum value in dataset out2d[:].full2d[:,:]
    maxi = maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)])

    # plot 2D spectra
    fig, ax = subplots(nrows,ncols,sharex=true,sharey=true,figsize=(ncols*3.2,nrows*3))
    if nplots == 1
        ax = [ax]
    end

    fig.suptitle(rep * " 2D spectrum (" * scal * ". scaling)")
    k = 0
    for i = 1:ncols
        for j = 1:nrows
            global k += 1
            if k > nplots
                continue
            end
            sca(ax[i,j])
            ax[i,j].set_aspect="equal"
            plot2d(spectra2d[k].ω,round.(spectra2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi)
            #plot2d(spectra2d[k].ω,conv2d;repr=rep,scaling=scal,norm=0)
            title("2D spectrum at $(T[k]) fs")
        end
    end
    tight_layout()
    subplots_adjust(top=0.8)

    ω = spectra2d[1].ω

    #DELETE or improve
    #plot([Ee-Eg, Ee-Eg] .+ 0*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([Ee-Eg, Ee-Eg] .+ 2*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([Ee-Eg, Ee-Eg] .+ 4*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 0*ωg,"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 2*ωg,"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 4*ωg,"k--",linewidth=0.5)

end
