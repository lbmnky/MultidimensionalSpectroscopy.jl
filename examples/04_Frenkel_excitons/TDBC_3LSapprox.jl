#!/usr/bin/julia
using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles

using Bridge, Random

# make sure to set script directory as pwd()
cd(@__DIR__)

# to use PyPlot in GUI mode under Linux
pygui(true)
# run once with calc_2d = false to initialize functions
calc_2d = true

cmp = create_colormap("bright");

Δ = 0.05

C = 2416 / 8065.5
# https://www.pks.mpg.de/fileadmin/user_upload/MPIPKS/group_pages/QuantumAggregates/PDF/pubilcations/Eisfeld_Briggs_2002_CP_61.pdf
# https://www.pks.mpg.de/~eisfeld/wwww/pdf/Eisfeld_Briggs_2006_CP_376.pdf

Emon = 2.41532 # 513 nm ... https://aip.scitation.org/doi/pdf/10.1063/1.469393
E = 2.08
E = Emon - C # Eisfeld
E2 = 2*E + Δ

# or analytically https://aip.scitation.org/doi/pdf/10.1063/1.469393
V = -2505 / 2
Eg_k1 = Emon + 2*V * cos(pi/(15+1)) / 8065
E = Eg_k1

Ek1_k1k2 = Emon + 2*V * cos(2*pi/(15+1)) / 8065

E2 = Eg_k1 + Ek1_k1k2



# make Hilbert-space
b_TLS = NLevelBasis(3)  # Hilbert-space of monomer

# create transition operators in monomer basis
j12 = transition(b_TLS,1,2)     # == g_e, or a / σ₋ (annihilation)
j21 = transition(b_TLS,2,1)     # == e_g, or a† /σ⁺ (creation)

j13 = transition(b_TLS,1,3)
j31 = transition(b_TLS,3,1)
j23 = transition(b_TLS,2,3)
j32 = transition(b_TLS,3,2)
# g_e / e_g notation from Constantin script

# define Hamiltonian of the dimer
# H₁ and H₂ are the Hamiltonians of the individual dipoles
H = E * j21 * j12 + E2 * j31 * j13 #+ .0 * (j12 + j21)


# get eigenstates and eigenvalues of H
#  using Julia LinearAlgebra functions instead of QuantumOptics requires to use H.data
states   = eigvecs(dense(H).data)
energies = eigvals(dense(H).data)

#states = eigenstates(dense(H))

# plot site basis and exciton basis energy diagrams
figure(figsize=(8,5))
ax2 = subplot(241)
#hlines(energies[:],-0.25, 0.25)
#xticks([-1, 0, 1], ["TLS"])
plot_levels(H,0)
hlines(.5*E2,-.45,.45, color="red",linestyle="dashed")
title("TLS energy levels")

# make initial density matrix
Psi0 = nlevelstate(b_TLS,1)
rho0 = dm(Psi0)
println("rho0")
display(rho0)

# time list for evaluation of master equation
Tmax = 550
dt = 1.
tlist = [0:dt:Tmax;]

#     relaxation        dephasing
#    _____________    _______________
L = [j12, j13, j23,  j21*j12, j31*j13]
Γ = [.01,.00 ,.05,    .13,     .01]
L = Γ .* L

μ12 = j12 + j21 + j13 + j31


display(dense(μ12))
# ... and with the doubly excited state:
μ23 = j23 + j32

rho1 = μ12 * rho0 * μ12
μ12  = μ12 / sqrt(tr(rho1))
rho1 = μ12 * rho0 * μ12
rho2 = μ23 * rho1 * μ23
#μ23  = μ23 / sqrt(tr(rho2)) * sqrt(2)
rho2 = μ23 * rho1 * μ23




N = 1
corr = sum(timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12) for i in 1:N) ./ N
#

D = .05
tc = 100
gt = D^2 * tc^2 * (exp.(-collect(-150:0.2:150)./tc) .+ collect(-150:0.2:150)./tc .- 1)
gt = D^2 * tc^2 * (exp.(-tlist./tc) .+ tlist./tc .- 1)

subplot(244)
plot(tlist,exp.(-gt))
##
#corr = corr .* exp.(-gt)

tout, rhot = timeevolution.master(tlist,μ12 * rho0,H,L)

# zeropad corr and extrapolate tlist
zp = 0
corr = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
ω,spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)

# add the results to previous figure
ax4 = subplot(244)
plot(tnew,real(corr),"k")
title("correlation function")
xlabel("time"); ylabel("Corr. func. ⟨μ(t) μ₀₁ ρ₀⟩")

ax5 = subplot(248)
plot(-[E, E], [0, .5])
plot(ω,spec,"k")
plot(-[energies[2], energies[2]], [0, .5],"k--")
title("TLS abs. spectrum")
xlabel("ω"); ylabel("Absorption")

tout, rhot = timeevolution.master(tlist,μ12 * rho0 * μ12,H,L)

vecs   = eigvecs(dense(H).data)
eivecs = [Ket(H.basis_l,vecs[:,i]) for i in 1:length(vecs[1,:])]

subplot(223)
title("eigen/coupled basis") # levels: 1 -> GS, 2 -> LP, 3 - N-1 -> ..., N -> UP
for i in 1:length(eivecs)
    es_n = real(expect(dm(eivecs[i]),rhot))
    plot(tout,es_n,linestyle="dashed", label=string(i)); #append!(leg,[string(i)])
    ylim((-0.1, 1))
end
legend()

try     # need to delete method first when working on it
    m = @which spectral_density(1)
    Base.delete_method(m)
catch
end
function spectral_density(ω) # power spectral density, thermal bath noise spectrum
    #w = [(ω .- ωq_TDBC[i]) ./ (.5 * width) for i in 1:length(ωq_TDBC)]
    b = 1 / 0.0259 # eV k_B*T at T=300 K
    f = .5 # 0 ... 1: relative importance of non-radiation vs. dephasing ... f = 0 > no dephasing
    γ = 0.7 # from absorption linewidth γ = γ_non-rad + γ_dephasing
    η = 1
    η = f * γ / (2 * pi * 0.0259)
    ω_cut = .1
    if ω == 0
        Jω = .01
    #elseif !isless(real(ω),0)
    #else
    #    w = [(ω .- ωq_TDBC[i]) ./ (.5 * width) for i in 1:length(ωq_TDBC)]
    #    Jω = sum([ampl_TDBC[i] ./ (1 .+ w[i].^2) for i in 1:length(ωq_TDBC)])  # probably most realistic ... ???
    #    #Jω = 1
    end
    if  !isless(real(ω),0) && ω != 0
        w = [(ω .- ωq_TDBC[i]) ./ (.5 * width) for i in 1:length(ωq_TDBC)]
        Jω = sum([ampl_TDBC[i] ./ (1 .+ w[i].^2) for i in 1:length(ωq_TDBC)])
        Jω = η * ω * exp(-(ω/ω_cut)^2) + Jω
        Jω = Jω * (1/exp(b *  ω - 1) + 1)
    elseif isless(real(ω),0)
        w = [(ω .- ωq_TDBC[i]) ./ (.5 * width) for i in 1:length(ωq_TDBC)]
        Jωa = sum([ampl_TDBC[i] ./ (1 .+ w[i].^2) for i in 1:length(ωq_TDBC)])
        Jω = η * -ω * exp(-(-ω/ω_cut)^2)
        Jω = Jω * (1/exp(b * -ω - 1)) + Jωa
    end
    if ω == 0
        Jω = 0.0000
    end
    return Jω

end

width = 0.0005
# from ... http://www.rsc.org/suppdata/c8/sc/c8sc00171e/c8sc00171e1.pdf
ωq_TDBC = [40, 80, 120, 150, 185, 197] ./ 8065 #.* 0.02998 # cm-1 in Hz #./ 8065 # in eV instead
ωq_TDBC = [0]
# => [0.00496 0.00992 0.01488 0.01860  0.02294 0.02443]
ampl_TDBC = [14, 18,  25,  43,  42, 67, 60] ./ 1000
ampl_TDBC = [0]


Γ = 50
a_ops = [Γ * j21 * j12, spectral_density, Γ/2 * j31 * j13, spectral_density]
#a_ops = [Γ * (j21+j31+j32) * (j12+j13+j23), spectral_density]
a_ops = [Γ * j21 * j12, spectral_density]

R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops]; J=L)

tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0,R,H)

#corr = expect(μ12,rhot)
#corr = corr .* exp.(-gt)

# calculate and plot spectrum
zp = 11
corr = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
ω, spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true);
subplot(248)
plot(ω,spec)
tight_layout()

subplot(244)
plot(tnew,corr)
#figure()
subplot(247)
plot(collect(-1:.001:1),(spectral_density.(collect(-1:.001:1))))

# load TA
TDBC_TA_t0 = readdlm("TDBC_TA_t0.dat",',')                  # https://aip.scitation.org/doi/pdf/10.1063/1.469393
TDBC_agg_abs = readdlm("TDBCaggr_abs.dat",',')
calc_ω = ω
calc_abs = spec

method = "redfield"
if method == "redfield"
    F       = R
    use_sub = false
elseif method == "lindblad"
    F       = L
    use_sub = false
end

#calc_2d = false
if calc_2d

        ## calculate (complex) 3rd order corr function (with T=0)
        zp = 11 # zeropad up to 2^zp

        ## calculate 2D spectra at
        T = [0] #0.66 fs

        spectra2d = Array{out2d}(undef, length(T))

        N = 50
        #E0 = [init_E() for i in 1:N]
        #H0 = Ham.(0,E0)
        #H0 = [Ham(50,0.1) for i in 1:N]
        #Threads.@threads for i = 1:length(T)
        for i = 1:length(T)
            #HH = [H0, Ham.(T[i],E0)]
            spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                                method;debug=false,use_sub=use_sub,zp=zp);
        end

        ## crop 2D data and increase dw
        spectra2d = [crop2d(spectra2d[i],1.5;w_max=3,step=1) for i = 1:length(T)]
end

if calc_2d
        ## plot 2D spectra for each(?) T
        # what to plot
        rep="abst"
        scal="lin"

        ## make  subplot layout
        nplots = length(T);
        ncols = Int32(ceil(sqrt(nplots)));
        nrows = Int32(ceil(nplots / ncols));

        # determine maximum value in dataset spectra2d[:].full2d[:,:]
        maxi = maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)])

        # plot 2D spectra
        fig, ax = subplots(nrows,ncols,sharex=true,sharey=true,figsize=(ncols*3.5,nrows*3))
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
        subplots_adjust(top=0.9)

        ω = spectra2d[1].ω
        ## plot additional things, like energy levels of states
        #plot([E, E], [ω[1], ω[end]],"k--",linewidth=1,alpha=0.25)
        #plot([ω[1], ω[end]],[E, E],"k--",linewidth=1,alpha=0.25)
        #plot([energies[2], energies[2]], [ω[1], ω[end]],"g--",linewidth=1,alpha=0.25)
        #plot([ω[1], ω[end]],[energies[2], energies[2]],"g--",linewidth=1,alpha=0.25)

        ## Save data
        #save_2d([round.(real(spectra2d[i].full2d),digits=1) for i = 1:length(T)],T,fn)

        ## plot TA (summed 2D spectrum)
        figure()
        ta = [sum(spectra2d[i].full2d,dims=1) for i in 1:length(spectra2d)] ./ length(spectra2d)
        plot(ω,vcat(ta...)')
        plot(ω,zeros(size(ω)),linestyle = "dashed")
        plot(1239.8 ./ TDBC_TA_t0[:,1], -1 * TDBC_TA_t0[:,2]/maximum(-TDBC_TA_t0[:,2]) * maximum(real(ta[1])),"o-")
        plot(TDBC_agg_abs[:,1],TDBC_agg_abs[:,2] ./ maximum(TDBC_agg_abs[:,2]) * maximum(real(ta[1])))
        plot(-calc_ω,calc_abs./maximum(calc_abs[-calc_ω.>0]).*maximum(real(ta[1])),"k--")
        xlabel("Energy/frequency")
        ylabel("Difference absorption")
        tight_layout()

end
