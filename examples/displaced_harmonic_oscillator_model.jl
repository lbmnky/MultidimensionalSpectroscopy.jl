#!/usr/bin/julia
using PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles

## some theory
# chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/13%3A_Coupling_of_Electronic_and_Nuclear_Motion/13.01%3A_The_Displaced_Harmonic_Oscillator_Model
# www.scm.com/doc/ADF/Input/Vibrationally_resolved_el_spectra.html

# make sure to set script directory as pwd()
cd(@__DIR__)

pygui(true)

# include my custom cmds module
if Sys.iswindows()
    include("..\\cmds.jl")
    fn = "01_Output"
else
    include("../cmds.jl")
    fn = "01_Output"
end

cmp = cmds.create_colormap("bright");

calc_2d = true

## Huang-Rhys factor
D_HR = .55
d = sqrt(D_HR)

b_tls = NLevelBasis(2)          # Hilbert-space of system                   Basis: {|ui⟩}
b_vib = FockBasis(1)            # Hilbert-space of oscillator               Basis: {|vi⟩}
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
Hg = Hg_el ⊗ one(b_vib) + j12*j21 ⊗ Hg_vib
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
He_vib = ωe * (number(b_vib) + identityoperator(b_vib) * 1/2)
He = He_el ⊗ one(b_vib) + j21*j12 ⊗ He_vib
#He = He_el ⊗ He_vib

## full system Hamiltonian
H = Hg + He

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

### ????? which one ??? #CHECK
D.data = FC

μ = transition(b_tls,1,2) + transition(b_tls,2,1)
μ = μ ⊗ D
μ12 = μ
μ23 = copy(μ)
fill!(μ23.data,0)

# initialize in groud state
Psi0 = nlevelstate(b_tls,1) ⊗ fockstate(b_vib,0)
rho0 = Psi0 ⊗ dagger(Psi0)  # alt.: rho0 = dm(Psi0)

## TODO used thermally populated ground state instead
#T = 0.01                               # Temperature T in what unit ? 
#rho0 = thermalstate(H,T)

#c_ops = [sqrt(0.05) * one(b_tls) ⊗ (a+at), sqrt(0.075) * (j21+j12) ⊗ one(b_vib)]
Γ = [sqrt(0.05), sqrt(0.075)]
L = Γ .* [one(b_tls) ⊗ a, j12 ⊗ one(b_vib)]

tlist = [0:0.05:15;]*2*π

corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)

#tout, rhot = timeevolution.master(tlist,μ12*rho0,H,c_ops)
#corr = expect(μ23,rhot)

subplot(222)
plot(tlist,real.(corr))

zp = 11
corr = cmds.zeropad(corr,zp)
tnew, ~ = cmds.interpt(tlist,zp)

# TO DO: correctly get absorption and emission!
ω_abs, spec_abs = timecorrelations.correlation2spectrum(tnew, corr;
                                                        normalize_spec=true)

corr = []
corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)
ω_em, spec_em = timecorrelations.correlation2spectrum(tlist, corr;
                                                        normalize_spec=true)

subplot(224)
plot(ω_abs,spec_abs)

temp = FC
temp = temp ./ maximum(real(temp))
for jj = 1:length(energies_g)
    plot(-([Ee-Eg, Ee-Eg] .+ 2*(jj-1)*ωg), [0, real(temp[jj,jj])], color="r",linewidth=1)
end

tight_layout()

F = L # works with method "lindblad"

if calc_2d

    ## zeropad up to 10^zp
    zp = 11

    ## calculate 2d spectra at
    T = [0.]

    out2d = Array{cmds.out2d}(undef, length(T))

    # use multithreding...
    Threads.@threads for i = 1:length(T)
    # ... or not
    #for i = 1:length(T)
        out2d[i] = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                            "lindblad";debug=true,use_sub=false,
                                                zp=zp, t2coh=false);
    end

    ## crop 2D data and increase dw
    out2d = [cmds.crop2d(out2d[i],1.3;w_max=1.89,step=1) for i = 1:length(T)]

    ## simulate effect of laser spectrum
    laserSpec =  .3 * exp.(-(out2d[1].ω.-1.35).^2/(2*(.19)^2)) +
                  exp.(-(out2d[1].ω.-1.6).^2/(2*(.11)^2)) +
                   .3 * exp.(-(out2d[1].ω.-1.85).^2/(2*(.05)^2))
    
    ## plot laser spectrum
    figure()
    plot(out2d[1].ω,laserSpec)

    ## 1D convoltion (in ω_excitation)
    IRF = (laserSpec.^2*(ones(length(laserSpec)))')
    ## 2D convolution
    IRF = (laserSpec.^2*(laserSpec.^2)')

    #for i in 1:length(out2d)
        conv2d = out2d[1].full2d .* IRF
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
    maxi = maximum([maximum(real(out2d[i].full2d)) for i in 1:length(out2d)])

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
            #cmds.plot2d(out2d[k].ω,round.(out2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi)
            cmds.plot2d(out2d[k].ω,conv2d;repr=rep,scaling=scal,norm=0)
            title("2D spectrum at $(T[k]) fs")
        end
    end
    tight_layout()
    subplots_adjust(top=0.8)

    ω = out2d[1].ω

    #DELETE or improve
    #plot([Ee-Eg, Ee-Eg] .+ 0*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([Ee-Eg, Ee-Eg] .+ 2*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([Ee-Eg, Ee-Eg] .+ 4*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 0*ωg,"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 2*ωg,"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 4*ωg,"k--",linewidth=0.5)

end
