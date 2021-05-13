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
D_HR = .25
d = sqrt(D_HR)
D_HR2 = .55
d2 = sqrt(D_HR2)

b_tls = NLevelBasis(2)          # Hilbert-space of system                   Basis: {|ui⟩}
b_vib = FockBasis(3)            # Hilbert-space of oscillator               Basis: {|vi⟩}
b     = b_tls ⊗ b_vib          # combined basis

j21 = transition(b_tls,2,1)                     # |e⟩⟨g|
j12 = dagger(j21)                               # |g⟩⟨e|
at  = create(b_vib)                             # ...
a   = dagger(at)
D   = displace(b_vib,d)
D2  = displace(b_vib,d2)


# ground states
Eg  = 0
ωg  = 0.10415
ωg2 = 0.16862
Hg_el = Eg * j12 * j21     # == HG = g * Eg ⊗ dagger(g)
#Hg_vib = ωg * (at*a + a*at * 1/2)
Hg_vib = ωg * (number(b_vib) + identityoperator(b_vib) * 1/2) 
Hg_vib2 = ωg2 * (number(b_vib) + identityoperator(b_vib) * 1/2)
#Hg_vib = ωg * (at*a + identityoperator(b_vib) * 1/2) 
#Hg_vib2 = ωg2 * (at*a + identityoperator(b_vib) * 1/2)
Hg = Hg_el ⊗ one(b_vib) ⊗ one(b_vib) + j12*j21 ⊗ (1 * Hg_vib) ⊗ one(b_vib) + j12*j21 ⊗ one(b_vib) ⊗ (1 * Hg_vib2)
#Hg = Hg_el ⊗ Hg_vib

# HAVE TO created excited states w/o displaced potential and then use displacement
# operator to act on transition dipole moment ... or calculate FC factors and use
# this as transition dipole operator ... #CHECK

# excited state
Ee  = 1.65
ωe  = .9 * 0.10415
ωe2 = .9 * 0.16862
He_el = Ee * j21 * j12
#He_vib = dagger(D) * Hg_vib * D # issue with D and dagger(D)
#He_vib = D' * Hg_vib * D # issue with D and dagger(D)
#He_vib2 = D' * Hg_vib * D
#He_vib = Hg_vib  # TODO figure out what's going on here
He_vib  = ( ωe  * (number(b_vib) + identityoperator(b_vib) * 1/2))
He_vib2 = ( ωe2 * (number(b_vib) + identityoperator(b_vib) * 1/2))
He = He_el ⊗ one(b_vib) ⊗ one(b_vib) + j21*j12 ⊗ (1 * He_vib) ⊗ one(b_vib) + j21*j12 ⊗ one(b_vib) ⊗ (1 * He_vib2)
#He = He_el ⊗ He_vib

## full system Hamiltonian
Het = (one(b_tls) ⊗ D ⊗ one(b_vib))' * He * (one(b_tls) ⊗ D ⊗ one(b_vib)) +
      (one(b_tls) ⊗ one(b_vib) ⊗ D2)' * He * (one(b_tls) ⊗ one(b_vib) ⊗ D2)

#He = (one(b_tls) ⊗ D ⊗ D2)' * He * (one(b_tls) ⊗ D ⊗ D2) 

H = Hg + He
H = (H + dagger(H)) / 2

## diagonalized Hamiltonian
states   = eigvecs(dense(H).data)
energies = eigvals(dense(H).data)

energies, states = eigenstates(dense(H))

# and GS
energies_g, states_g = eigenstates(dense(Hg_vib))
energies_g2, states_g2 = eigenstates(dense(Hg_vib2))
# and ES
energies_e, states_e = eigenstates(dense(D' * He_vib * D))
energies_e2, states_e2 = eigenstates(dense(D2' * He_vib2 * D2))


## plot energy levels
N = length(energies)
println(Int8(ceil(N/2)))
figure(figsize=(6,4))
subplot(121)
scatter(zeros(length(energies_g)), energies_g)
scatter(zeros(length(energies_g2)), energies_g2)

x = [-1:0.01:2*d;]

# TODO: transform to position basis to plot PECs
#xx = PositionBasis(0,10,200)
#xxx = position(xx)

m = 50                        # arbitrary value for m
f = m * ωg^2
f2 = m * ωg2^2
# GS
plot(x,f*x.^2)
plot(x,f2*x.^2)
# ES
f = m * ωe^2
f2 = m * ωe2^2
scatter(zeros(length(energies_e)) .+ d, energies_e .+ Ee)
scatter(zeros(length(energies_e2)) .+ d2, energies_e2 .+ Ee)
plot(x,f * (x.-d).^2 .+ Ee)
plot(x,f2 * (x.-d2).^2 .+ Ee)
ylim(0, 3.5)

#for i=1:length(states_g)
#    plot(xpoints, abs2.(states_g[i].data).*40 .+ energies_g[i],color="k",linewidth=1)
#    plot(xpoints, abs2.(states_e[i].data).*40 .+ energies_e[i],color="g",linewidth=1)
#end

## calculate Franck-Condon factors
# |⟨v|i⟩² , where |v⟩ is the vibrational level on the excited electronic state and |i⟩ the initial state
FC = complex(zeros(length(energies_g),length(energies_e)))
for ii = 1:length(energies_g)
    for jj = 1:length(energies_e)
        #FC[ii,jj] = abs(conj(states_g[1]' * (D*states_e[jj])) * states_g[1]' * (D*states_e[ii]))
        #FC[ii,jj] = abs(conj((D*states_e[jj])' * states_g[1]) * (D*states_e[ii])' * states_g[1])
        FC[ii,jj] = abs(states_e[jj]' * states_g[1])^2
        #FC[ii,jj] = abs(conj(states_g[1]' * states_e[jj]) * (states_g[1]' * states_e[ii]))
    end
end

FC2 = complex(zeros(length(energies_g2),length(energies_e2)))
for ii = 1:length(energies_g2)  # only for lowest vib level on ground state ... length(energies_g2)
    for jj = 1:length(energies_e2)
        #FC2[ii,jj] = abs(conj(states_g2[1]' * (D2*states_e2[jj])) * states_g2[1]' * (D2*states_e2[ii]))
        #FC2[ii,jj] = abs(conj((D2*states_e2[jj])' * states_g2[1]) * (D2*states_e2[ii])' * states_g2[1])
        FC2[ii,jj] = abs(states_e2[jj]' * states_g2[1])^2
    end
end

## RAMAN overlaps
## calculate Franck-Condon factors
# ⟨f|v⟩⟨v|i⟩
FR = complex(zeros(length(energies_g),length(energies_e)))
for ii = 1:length(energies_g)
    for jj = 1:length(energies_e)
        FR[ii,jj] = abs(conj(states_g[2]' * D*states_e[jj]) * ( (D*states_e[ii])' * states_g[1]))^2
    end
end

## RAMAN overlaps
## calculate Franck-Condon factors
# ⟨f|v⟩⟨v|i⟩
FR2 = complex(zeros(length(energies_g),length(energies_e)))
for ii = 1:length(energies_g)
    for jj = 1:length(energies_e)
        FR2[ii,jj] = abs(conj(states_g2[2]' * D*states_e2[jj]) * ( (D*states_e2[ii])' * states_g2[1]))^2
    end
end


# the previous algorithm corresponds to calculating the FC-factors
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

μ12 = j12 ⊗ D ⊗ D2# + j21 ⊗ D' ⊗ D2'
#μ12 = j12 ⊗ D ⊗ one(b_vib) + j12 ⊗ one(b_vib) ⊗ D2
μ12 = μ12 + dagger(μ12)
#μ12 = μ ⊗ D ⊗ one(b_vib) + μ ⊗ one(b_vib) ⊗ D2
#μ12 = μ ⊗ (a*at) ⊗ (a*at)  # ???
#μ12 = rho0 * μ12 + μ12 * rho0
#μ12 = rho0 * (μ ⊗ D ⊗ D2)
#μ12 = μ12 + μ12'
#μ12 = μ ⊗ (a + at) ⊗ (a + at)
μ23 = copy(μ12)
fill!(μ23.data,0)

# initialize in groud state
Psi0 = nlevelstate(b_tls,1) ⊗ fockstate(b_vib,0) ⊗ fockstate(b_vib,0)
rho0 = Psi0 ⊗ dagger(Psi0)  # alt.: rho0 = dm(Psi0)

Psi1 = nlevelstate(b_tls,2) ⊗ fockstate(b_vib,0) ⊗ fockstate(b_vib,0) 
rho1 = dm(Psi1)

 ## TODO used thermally populated ground state instead
#T = 0.01                               # Temperature T in what unit ? 
#rho0 = thermalstate(H,T)

Γ = [sqrt(0.0001), sqrt(0.0001), sqrt(0.045)]
L = Γ .* [one(b_tls) ⊗ a ⊗ one(b_vib), one(b_tls) ⊗ one(b_vib) ⊗ a, j21*j12 ⊗ one(b_vib) ⊗ one(b_vib)]

tlist = [0:0.1:30;]*2*π

corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)

# calculate ⟨f|i(t)⟩ from https://pubs.acs.org/doi/pdf/10.1021/ar960240c eq.10
#corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, j12⊗a⊗one(b_vib))
#figure()
#plot(tlist,corr)


subplot(222)
plot(tlist,real.(corr))

zp = 11
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

ω_abs, spec_abs = timecorrelations.correlation2spectrum(tnew, corr;
                                                        normalize_spec=true)

corr = []
corr = timecorrelations.correlation(tlist, rho1, H, L, μ12, μ12)
#tout, rhot = timeevolution.master(tlist,rho1*μ12,H,L)
#corr = real(expect(μ12,rhot))
zp = 11
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

ω_em, spec_em = timecorrelations.correlation2spectrum(tnew, corr;
                                                        normalize_spec=true)

subplot(224)
plot(ω_abs,spec_abs)
plot(-ω_em,spec_em,"r:")

legend(["absorption", "emission"])


γ = 0.02 # γ has an effect on the intensity envelop
abs_crossSection = real(ω_abs .* sum([γ * FC[1,i]./((i * ωe + Ee - ωg .- ω_abs).^2 .- (im * γ)^2) for i in 1:length(energies_e)]))
plot(-ω_abs,abs_crossSection./maximum(abs_crossSection)./2)

abs_crossSection2 = real(ω_abs .* sum([γ * FC2[1,i]./((i * ωe2 + Ee - ωg2 .- ω_abs).^2 .- (im * γ)^2) for i in 1:length(energies_e)]))
plot(-ω_abs,abs_crossSection2./maximum(abs_crossSection2)./2)

abs_crossSec = (abs_crossSection + abs_crossSection2) ./maximum(abs_crossSection + abs_crossSection2)

plot(-ω_abs,abs_crossSec)

#temp = FC
#temp = temp ./ maximum(real(temp))
#for jj = 1:length(energies_g)
#    plot(-([Ee-Eg, Ee-Eg] .+ (jj-1)*ωg), [0, real(temp[jj,jj])], color="r",linewidth=2)
#end

#temp = FC2
#temp = temp ./ maximum(real(temp))
#for jj = 1:length(energies_g)
#    plot(-([Ee-Eg, Ee-Eg] .+ (jj-1)*ωg2), [0, real(temp[jj,jj])], color="b",linewidth=1)
#end
#=
temp = FR
temp = temp ./ maximum(real(temp))
for jj = 1:length(energies_g)
    plot(-([Ee-Eg, Ee-Eg] .+ (jj-1)*ωg) .+ 0.005, [0, real(temp[jj,jj])], color="y",linewidth=2,linestyle=":")
end

temp = FR2
temp = temp ./ maximum(real(temp))
for jj = 1:length(energies_g)
    plot(-([Ee-Eg, Ee-Eg] .+ (jj-1)*ωg2) .+ 0.005, [0, real(temp[jj,jj])], color="g",linewidth=2,linestyle=":")
end
=#
tight_layout()

F = L # works with method "lindblad"

if calc_2d

    ## zeropad up to 10^zp
    zp = 11

    ## calculate 2d spectra at
    T = [0.]

    spectra2d = Array{out2d}(undef, length(T))

    # use multithreding...
    Threads.@threads for i = 1:length(T)
    # ... or not
    #for i = 1:length(T)
        spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                            "lindblad";debug=true,use_sub=false,
                                                zp=zp, t2coh=false);
    end

    ## crop 2D data and increase dw
    spectra2d = [crop2d(spectra2d[i],1;w_max=4,step=1) for i = 1:length(T)]

    ## simulate effect of laser spectrum
    laserSpec =  .3 * exp.(-(spectra2d[1].ω.-1.35).^2/(2*(.19)^2)) +
                  exp.(-(spectra2d[1].ω.-1.6).^2/(2*(.11)^2)) +
                   .3 * exp.(-(spectra2d[1].ω.-1.85).^2/(2*(.05)^2))
    
    ## plot laser spectrum
    #figure()
    #plot(spectra2d[1].ω,laserSpec)

    ## 1D convoltion (in ω_excitation)
    IRF = (laserSpec.^2*(ones(length(laserSpec)))')
    ## 2D convolution
    IRF = (laserSpec.^2*(laserSpec.^2)')

    #for i in 1:length(out2d)
        conv2d = spectra2d[1].full2d .* IRF
    #end
end

if calc_2d

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
