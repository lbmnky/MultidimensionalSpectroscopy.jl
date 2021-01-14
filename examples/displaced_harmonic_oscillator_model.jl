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

D_HR = .55
d = sqrt(D_HR)

b_tls = NLevelBasis(2)  # Hilbert-space of system
b_vib = FockBasis(3)    # Hilbert-space of oscillator
b = b_tls ⊗ b_vib       # combined basis

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
j21 = transition(b_tls,2,1)
j12 = dagger(j21)
at = create(b_vib)
a = dagger(at)
D = displace(b_vib,d)

# ground states
Eg = 0
ωg = 0.2
Hg_el = Eg * j12 * j21     # == HG = g * Eg ⊗ dagger(g)
Hg_vib = ωg * (at*a + one(b_vib) * 1/2)
Hg = Hg_el ⊗ one(b_vib) + one(b_tls) ⊗ Hg_vib
#Hg = Hg_el ⊗ Hg_vib

# HAVE TO created excited states w/o displaced potential and the use displacement
# operator to act on transition dipole moment ... or calculate FC factors and use
# this as transition dipole operator ...

# excited state
Ee = 4
He_el = Ee * j21 * j12
He_vib = dagger(D) * Hg_vib * D # issue with D and dagger(D)
#He_vib = D * Hg_vib * dagger(D) # issue with D and dagger(D)
He_vib = Hg_vib
He = He_el ⊗ one(b_vib) + one(b_tls) ⊗ He_vib
#He = He_el ⊗ He_vib
H = Hg + He
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
println(Int8(N/2))
figure(figsize=(6,4))
subplot(121)
scatter(zeros(Int8(N/2)), energies[1:Int8(N/2)])
x = [-1:0.1:2*d;]

#xx = PositionBasis(0,10,200)
#xxx = position(xx)

f = sqrt(1/(2*ωg))
plot(x,f*x.^2)
scatter(zeros(Int8(N/2)) .+ d, energies[Int8(N/2)+1:end])
plot(x,f * (x.-d).^2 .+ Ee)
#plot(xxx,xxx)

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

μ = transition(b_tls,1,2) + transition(b_tls,2,1)
μ = μ ⊗ D
μ12 = μ
μ23 = μ

# initialize in groud state
Psi0 = nlevelstate(b_tls,1) ⊗ fockstate(b_vib,0)
rho0 = Psi0 ⊗ dagger(Psi0)  # Or do: rho0 = dm(Psi0)

# thermal density matrix ??? Must only populate vibrations
#T = 0.02
#rho0 = thermalstate(Hg,T)

c_ops = [sqrt(0.05) * one(b_tls) ⊗ (a+at), sqrt(0.075) * (j21+j12) ⊗ one(b_vib)]
#c_ops = [sqrt(0.025) * one(b_tls) ⊗ (a+at)]
#c_ops = []

tlist = [0:0.05:15;]*2*π

corr = timecorrelations.correlation(tlist, rho0, H, c_ops, μ12, μ23)

tout, rhot = timeevolution.master(tlist,μ12*rho0,H,c_ops)
corr = expect(μ23,rhot)

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
    T = [0.,10,20]
    #out2d = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T,"lindblad";debug=false,zp=zp);
    out2d = Array{cmds.out2d}(undef, length(T))
    for i = 1:length(T)
        out2d[i] = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                            "lindblad";debug=true,use_sub=false,
                                                zp=zp);
    end

    ## crop 2D data and increase dw
        out2d = [cmds.crop2d(out2d[i],3.5;w_max=5.5,step=2) for i = 1:length(T)]

    ## assign ω-axis from output
    ω = out2d[1].ω

    ## plot 2D spectra for each(?) T
    # what to plot
    rep="absorptive"
    scal="lin"

    ## make  subplot layout
    nplots = length(T); ncols = ceil(sqrt(nplots)); nrows = ceil(nplots / ncols)

    ## create figure
    figure(figsize=(ncols*3+0.25,nrows*3+0.5));
    suptitle(rep * " 2D spectrum (" * scal * ". scaling)")
    for i = 1:nplots
            subplot(nrows,ncols,i,aspect="equal")
            cmds.plot2d(ω,out2d[i].full2d;repr=rep,scaling=scal)
            title("2D spectrum at $(T[i]) fs")
            #xlim([0, E₁+E₂]); ylim([0, E₁+E₂]);
    end
    tight_layout()

    plot([Ee-Eg, Ee-Eg] .+ 0*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    plot([Ee-Eg, Ee-Eg] .+ 2*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    plot([Ee-Eg, Ee-Eg] .+ 4*ωg,[ω[1], ω[end]],"k--",linewidth=0.5)
    plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 0*ωg,"k--",linewidth=0.5)
    plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 2*ωg,"k--",linewidth=0.5)
    plot([ω[1], ω[end]],[Ee-Eg, Ee-Eg] .+ 4*ωg,"k--",linewidth=0.5)

    # choose one of the following plotting ranges
    #xlim([ω[1], ω[end]]); ylim([ω[1], ω[end]]);
    #xlim([0, ω[end]]); ylim([0, ω[end]]);
    #xlim([0, 2*(Eg+Ee)]); ylim([0, 2*(Eg+Ee)]);
end

show()
