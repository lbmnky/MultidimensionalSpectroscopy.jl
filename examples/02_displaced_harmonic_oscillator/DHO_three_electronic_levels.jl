using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles

## good source of information 
# chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/13%3A_Coupling_of_Electronic_and_Nuclear_Motion/13.01%3A_The_Displaced_Harmonic_Oscillator_Model
# www.scm.com/doc/ADF/Input/Vibrationally_resolved_el_spectra.html

# make sure to set script directory as pwd()
cd(@__DIR__)

# let pyplot make popup figures
pygui(true)

# custom colormap (neg-0-pos : purple-blue-white-green-yellow-red)
cmp = create_colormap("bright");

# do you want to calculate 2D spectra ?
calc_2d = true

# Huang-Rhys factor (excited state; w.r.t. ground state)
D_HR = .3
# displacement
d = sqrt(D_HR)

# Huang-Rhys factor (second excited state for ESA)
D_HRf = .01
df = sqrt(D_HRf)

# create Hilbert spaces
b_3ls = NLevelBasis(3)   # Hilbert-space of system  (3 electronic levels)                Basis: {|ui⟩}
b_vib = FockBasis(1)     # Hilbert-space of oscillator (1, 2, ... vibrational states)    Basis: {|vi⟩}
b = b_3ls ⊗ b_vib       # combined basis

# operators of 3 level system
j21 = transition(b_3ls,2,1)                     # |e⟩⟨g|
j12 = dagger(j21)                               # |g⟩⟨e|
j32 = transition(b_3ls,3,2)                     # |e⟩⟨f|
j23 = dagger(j32)                               # |f⟩⟨e|

# operators of vibrational mode
at = create(b_vib)                              # ...
a  = dagger(at)
D  = displace(b_vib,d)                          # displacement operator for excited state
Df = displace(b_vib,df)                         # displacement operator for second excited state

# ground states Hamiltonian
Eg = 0                                          # energy of electronic level 
ωg = 0.15                                       # energy of vibration
Hg_el = Eg * j12 * j21                          # electronic part of Hamiltonian        == HG = g * Eg ⊗ dagger(g)
#Hg_vib = ωg/2 * (at*a + one(b_vib))
Hg_vib = ωg * (number(b_vib) + identityoperator(b_vib) * 1/2)   # vibrational part of H

#Hg = Hg_el ⊗ one(b_vib) + one(b_3ls) ⊗ Hg_vib
Hg = Hg_el ⊗ one(b_vib) + j12*j21 ⊗ Hg_vib    # combined H for ground state -- H = |0⟩⟨0| + #TODO fill in
#Hg = Hg_el ⊗ Hg_vib

# #TODO Check: HAVE TO created excited states w/o displaced potential and the use displacement
# operator to act on transition dipole moment ... or calculate FC factors and use
# this as transition dipole operator ...

## excited state Hamiltonian
Ee = 1.65                                       # energy of excited electronic level
ωe = 0.15                                       # energy of vibration
He_el = Ee * j21 * j12                          # electronic part of H
#He_vib = dagger(D) * Hg_vib * D # issue with D and dagger(D)
#He_vib = D * Hg_vib * dagger(D) # issue with D and dagger(D)
He_vib = ωe * (number(b_vib) + identityoperator(b_vib) * 1/2)    #vibrational part of H
#He = He_el ⊗ one(b_vib) + one(b_3ls) ⊗ He_vib
He = He_el ⊗ one(b_vib) + j21*j12 ⊗ He_vib
#He = He_el ⊗ He_vib

## 2nd excited state Hamiltonian
Ef = 2 * Ee + 0.3 #4.45                                                       # energy of second exc. el. level
ωf = 0.15                                                               # energy of vibration
Hf_el = Ef * j32 * j23                                                  # el. part of H
Hf_vib = ωf * (number(b_vib) + identityoperator(b_vib) * 1/2)           # vib. part of H
Hf = Hf_el ⊗ one(b_vib) + j32*j23 ⊗ Hf_vib

## construct system Hamiltonian
H = Hg + He + Hf
#println("Display Hamiltonian")
#display(dense(H).data[1:10,1:4])

## diagonalize system Hamiltonian
states   = eigvecs(dense(H).data)
energies = eigvals(dense(H).data)

# same for ground state
energies_g, states_g = eigenstates(dense(Hg_vib))
# and excited state
energies_e, states_e = eigenstates(dense(He_vib))

## number of eigentstates
N = length(energies)

## plot energy levels
figure(figsize=(6,4))
subplot(121)
# visualize ground state
scatter(zeros(Int8(ceil(N/3))), energies[1:Int8(ceil(N/3))])
# and arbitrary potential energy surface (harmonic)
x = [-1:0.01:2*d;]
m = 100                                         # arbitrary
f = 1/2 * m * ωg^2
plot(x,f*x.^2,"b")

#TODO use transformation into position basis instead
#xx = PositionBasis(0,10,200)
#xxx = position(xx)
#### go to position basis ?? ?
#xpoints = samplepoints(b)

# visualize excited state
m = 100                                          # arbitrary
f = 1/2 * m * ωe^2
scatter(zeros(Int8(floor(N/3))) .+ d, energies[Int8(floor(N/3))+1:end-Int8(floor(N/3))])
plot(x,f * (x.-d).^2 .+ Ee,"y")

# visualize second exc. state
m = 100                                         # arbitrary
f = 1/2 * m * ωf^2
scatter(zeros(Int8(floor(N/3))) .+ df, energies[end-Int8(floor(N/3))+1:end])
plot(x,f * (x.-df).^2 .+ Ef,"g"); ylim(0, 6)

## ????
#for i=1:length(states_g)
#    plot(xpoints, abs2.(states_g[i].data).*40 .+ energies_g[i],color="k",linewidth=1)
#    plot(xpoints, abs2.(states_e[i].data).*40 .+ energies_e[i],color="g",linewidth=1)
#end

## calculate Franck-Condon factors between ground and excited state
FC = complex(zeros(length(energies_g),length(energies_e)))
for ii = 1:length(energies_g)
    for jj = 1:length(energies_e)
        FC[ii,jj] = abs(conj(states_g[1]' * (D*states_e[jj])) * states_g[1]' * (D*states_e[ii]))
        #FC[ii,jj] = abs(conj(states_g[1]' * states_e[jj]) * (states_g[1]' * states_e[ii]))                 # this with D somewhere else ? 
    end
end

## ?????? #TODO CHECK!!! the previous algorith corresponds to calculating the FC-factors
l = complex(zeros(length(energies_g),length(energies_e)))
for ii = 1:length(energies_g)
    for jj = 1:length(energies_e)
        l[ii,jj] = exp(-d^2) * d^(2*jj) / factorial(jj)
    end
end

## ??????
ln = zeros(5,)
for n = 1:5
    ln[n] = exp(-d^2) * d^(2*n) / factorial(n)
end

## # TODO what's going on here ?
D.data = FC

## TDM operator between ground and excited state in NLevelBasis only
μ12 = j12 + j21
μ12 = μ12 ⊗ D

## TDM operator between excited and second excited state in NLevelBasis only
μ23 = j23 + j32
μ23 = μ23 ⊗ Df

## initialize system in the ground state
Psi0 = nlevelstate(b_3ls,1) ⊗ fockstate(b_vib,0)
rho0 = Psi0 ⊗ dagger(Psi0)                                                                  # alt.: rho0 = dm(Psi0)

## TODO: used thermally populated ground state instead
#T = 0.01                                                               # Temperature T in what unit ? 
#rho0 = thermalstate(H,T)

## Lindblad operators with rates Γ
Γ = [sqrt(0.05), sqrt(0.065)]
L = Γ .* [one(b_3ls) ⊗ (a), (j12) ⊗ one(b_vib)]

## time dimension
tlist = [0:0.12:20;]*2*π

## calculate correlation function for ground/excited state transition
corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)

## could do the same in the following way
#tout, rhot = timeevolution.master(tlist,μ12*rho0,H,L)
#corr = expect(μ23,rhot)

## plot correlation function
subplot(222)
plot(tlist,real.(corr))

## zeropad to smooth data
zp = 11
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

#TODO : correctly get absorption and emission! (with correct frequency sign)
## get absorption spectrum
ω_abs, spec_abs = timecorrelations.correlation2spectrum(tnew, corr;
                                                        normalize_spec=true)

## get emission spectrum (just sign switch so far ... calc relaxation in excited state ?)                                                    
corr = []
corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)
ω_em, spec_em = timecorrelations.correlation2spectrum(tlist, corr;
                                                        normalize_spec=true)


## plot absorption spectrum                                                    
subplot(224)
plot(ω_abs,spec_abs)

## plot stick spectrum of Franck-Condon factors
temp = FC
temp = temp ./ maximum(real(temp))
for jj = 1:length(energies_g)
    plot(-([Ee-Eg, Ee-Eg] .+ 2*(jj-1)*ωg), [0, real(temp[jj,jj])], color="r",linewidth=1)
end

tight_layout()

## calculate 2D spectra
if calc_2d
    F = L                                           # works with method "lindblad"
    zp = 10                                         # zeropad data up to 10^zp values

    ## calculate 2d spectra at
    T = [0.]

    ## initialize output array/structure    
    spectra2d = Array{out2d}(undef, length(T))

    ## use multithreading, if you can:
    Threads.@threads for i = 1:length(T)
    # ... or else:
    #for i = 1:length(T)
        spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                            "lindblad";debug=true,use_sub=false,
                                                zp=zp, t2coh=false);
    end

    ## crop 2D data and increase dw (step>1) if necessary
    spectra2d = [crop2d(spectra2d[i],1.4;w_max=2.2,step=1) for i = 1:length(T)]

end

if calc_2d
    
    ## plot 2D spectra for each(?) T
    # decide what to plot
    rep  = "absorptive"                             # "absorptive", "dispersive", "absolute"
    scal = "lin"                                    # "lin", "tan", "cube", "asinh"

    ## make subplot layout
    nplots = length(T);                             # total number of 2D spectra
    ncols  = Int32(ceil(sqrt(nplots)));
    nrows  = Int32(ceil(nplots / ncols));

    ## determine maximum value in dataset out2d[:].full2d[:,:]
    maxi = maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)])

    ## plot 2D spectra
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
            title("2D spectrum at $(T[k]) fs")
        end
    end
    tight_layout()
    subplots_adjust(top=0.85)

    ## #TODO: remove the following of make more universal
    ω = spectra2d[1].ω

    #plot([Ee-Eg, Ee-Eg] .+ 0*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([Ee-Eg, Ee-Eg] .+ 2*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([Ee-Eg, Ee-Eg] .+ 4*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 0*ωg,"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 2*ωg,"k--",linewidth=0.5)
    #plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 4*ωg,"k--",linewidth=0.5)

end
