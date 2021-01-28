#!/usr/bin/julia
using PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles
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

D_HR = .5
d = sqrt(D_HR)


D_HRf = .05
df = sqrt(D_HRf)


b_3ls = NLevelBasis(3)  # Hilbert-space of system                   Basis: {|ui⟩}
b_vib = FockBasis(1)    # Hilbert-space of oscillator               Basis: {|vi⟩}
b = b_3ls ⊗ b_vib       # combined basis

#=
# operators direct in b
# transition operators system
j21 = transition(b_tls,2,1) ⊗ one(b_vib)
j12 = dagger(j21)

# creation, annihilation operator vib mode
at = one(b_tls) ⊗ create(b_vib)
a = dagger(at)

# displacement is in b_vib only
D = one(b_tls) ⊗ displace(b_vib,d)
=#
j21 = transition(b_3ls,2,1)                     # |e⟩⟨g|
j12 = dagger(j21)                               # |g⟩⟨e|
j32 = transition(b_3ls,3,2)                     # |e⟩⟨f|
j23 = dagger(j32)                               # |f⟩⟨e|
at = create(b_vib)                              # ...
a = dagger(at)
D = displace(b_vib,d)
Df = displace(b_vib,df)

# ground states
Eg = 0
ωg = .2
Hg_el = Eg * j12 * j21     # == HG = g * Eg ⊗ dagger(g)
#Hg_vib = ωg/2 * (at*a + one(b_vib))
Hg_vib = ωg/2 * (2*at*a + a*at)

#Hg = Hg_el ⊗ one(b_vib) + one(b_3ls) ⊗ Hg_vib
Hg = Hg_el ⊗ one(b_vib) + j12*j21 ⊗ Hg_vib
#Hg = Hg_el ⊗ Hg_vib

# HAVE TO created excited states w/o displaced potential and the use displacement
# operator to act on transition dipole moment ... or calculate FC factors and use
# this as transition dipole operator ...

# excited state
Ee = 2.05
ωe = .2
He_el = Ee * j21 * j12
#He_vib = dagger(D) * Hg_vib * D # issue with D and dagger(D)
#He_vib = D * Hg_vib * dagger(D) # issue with D and dagger(D)
He_vib = ωe/2 * (2*at*a + a*at)
#He = He_el ⊗ one(b_vib) + one(b_3ls) ⊗ He_vib
He = He_el ⊗ one(b_vib) + j21*j12 ⊗ He_vib
#He = He_el ⊗ He_vib

# 2nd excited state
Ef = 4.6
ωf = .3
Hf_el = Ef * j32 * j23
#Hf_vib = dagger(D) * He_vib * D
Hf_vib = ωf/2 * (2*at*a + a*at)
#Hf = Hf_el ⊗ one(b_vib) + one(b_3ls) ⊗ Hf_vib
Hf = Hf_el ⊗ one(b_vib) + j32*j23 ⊗ Hf_vib


H = Hg + He + Hf
println("Display Hamiltonian")
#display(dense(H).data[1:10,1:4])

# diagonalized Hamiltonian
states = eigvecs(dense(H).data)
energies = eigvals(dense(H).data)

energies_g, states_g = eigenstates(dense(Hg_vib))

energies_e, states_e = eigenstates(dense(He_vib))

# Display
#display((states))
#display(energies)
N = length(energies)
figure(figsize=(6,4))
subplot(121)
scatter(zeros(Int8(ceil(N/3))), energies[1:Int8(ceil(N/3))])
x = [-1:0.1:2*d;]
m = 100            # arbitrary
f = 1/2 * m * ωg^2
plot(x,f*x.^2,"b")

#xx = PositionBasis(0,10,200)
#xxx = position(xx)

m = 100            # arbitrary
f = 1/2 * m * ωe^2
scatter(zeros(Int8(floor(N/3))) .+ d, energies[Int8(floor(N/3))+1:end-Int8(floor(N/3))])
plot(x,f * (x.-d).^2 .+ Ee,"y")
ylim(0, 6)



# f
m = 100            # arbitrary
f = 1/2 * m * ωf^2
scatter(zeros(Int8(floor(N/3))) .+ df, energies[end-Int8(floor(N/3))+1:end])
plot(x,f * (x.-df).^2 .+ Ef,"g")
ylim(0, 6)

#### go to position basis ?? ?
#xpoints = samplepoints(b)

#for i=1:length(states_g)
#    plot(xpoints, abs2.(states_g[i].data).*40 .+ energies_g[i],color="k",linewidth=1)
#    plot(xpoints, abs2.(states_e[i].data).*40 .+ energies_e[i],color="g",linewidth=1)
#end

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

### ????? which one ???
D.data = FC

#μ12 = transition(b_3ls,1,2) + transition(b_3ls,2,1)
μ12 = j12 + j21
μ12 = μ12 ⊗ D

#μ12 = transition(b_3ls,1,2) + transition(b_3ls,2,1)
μ23 = j23 + j32
μ23 = μ23 ⊗ Df

# initialize in groud state
Psi0 = nlevelstate(b_3ls,1) ⊗ fockstate(b_vib,0)
rho0 = Psi0 ⊗ dagger(Psi0)  # Or do: rho0 = dm(Psi0)

# thermal density matrix ??? Must only populate vibrations
#T = 0.01
#rho0 = thermalstate(Hg,T)

#c_ops = [sqrt(0.05) * one(b_tls) ⊗ (a+at), sqrt(0.075) * (j21+j12) ⊗ one(b_vib)]
Γ = [sqrt(0.05), sqrt(0.175)]
c_ops = [Γ[1] * one(b_3ls) ⊗ (a), Γ[2] * (j12) ⊗ one(b_vib)]

tlist = [0:0.05:15;]*2*π

corr = timecorrelations.correlation(tlist, rho0, H, c_ops, μ12, μ12)

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
corr = timecorrelations.correlation(tlist, rho0, H, c_ops, μ12, μ12)
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
show()


if calc_2d
    F = c_ops # works with method "lindblad"
    zp = 10

    ## calculate 2d spectra at
    T = [0.]
    #out2d = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T,"lindblad";debug=false,zp=zp);
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
    out2d = [cmds.crop2d(out2d[i],1.5;w_max=3,step=2) for i = 1:length(T)]

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
            cmds.plot2d(out2d[k].ω,round.(out2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi)
            #cmds.plot2d(out2d[k].ω,abs.(out2d[k].full2d_nr);repr=rep,scaling=scal,norm=0)
            title("2D spectrum at $(T[k]) fs")
            colorbar()
        end
    end
    tight_layout()
    subplots_adjust(top=0.8)

    ω = out2d[1].ω

    plot([Ee-Eg, Ee-Eg] .+ 0*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    plot([Ee-Eg, Ee-Eg] .+ 2*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    plot([Ee-Eg, Ee-Eg] .+ 4*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 0*ωg,"k--",linewidth=0.5)
    plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 2*ωg,"k--",linewidth=0.5)
    plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 4*ωg,"k--",linewidth=0.5)

end
