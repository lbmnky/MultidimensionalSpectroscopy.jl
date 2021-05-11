# https://www.nature.com/articles/s41598-017-16178-8.pdf?proof=t
using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles, Random, Combinatorics, Distributions

QuantumOpticsBase.set_printing(standard_order=false,rounding_tol=1e-3)

# make sure to set script directory as pwd()
cd(@__DIR__)

mod_name = "A"
fn_log =  "logs/" * mod_name * ".log"

# erase .log file
try
    rm("logs/" * mod_name * ".log")
    rm(fn_log)
catch
end

### -- select model parameters --- ###
mod_name = "parameters/" * mod_name
include( mod_name * ".params")





logger(String(read(mod_name * ".params")), fn_log)

pygui(true)

calc_2d = true

logger("\nDo 2D calculation: ", calc_2d, fn_log)

cmp = create_colormap("bright");

# use 2 level system to approximate Frenkel excitons ... seems to coarse, but well
Nlev = 3

b = NLevelBasis(Nlev)

# energy of TDBC Frenkel exciton bright state #TODO: add some details, such as J angles, etc.
Emon = E_monomer
Emonb = 2 * E_monomer + 0.1
# not important here, but in general for TDBC aggregates
# Δ = .299                          # shift monomer -> aggregate 
# J = -Δ / (2 * cos(pi /(N+1)))     # nearest-neighbour coupling constant

# transitions and monomer Hamiltonian
σ⁺ = transition(b,2,1) 
σ31 = transition(b,3,1) 
σ32 = transition(b,3,2)
σ⁻ = transition(b,1,2) 
σ13 = transition(b,1,3) 
σ23 = transition(b,2,3)
H_mon = Emon * σ⁺ * σ⁻ + Emonb * σ31*σ13

σ⁺ = σ⁺ + σ32
σ⁻ = σ⁻ + σ23

logger("\nMonomer Hamiltonian:\n\n", dense(H_mon), fn_log)

# use N uncoupled TLSs
N = N_monomers      # from .params file
B = b^N             # basis for matter part of H

# random (normal dist) list of energy shifts (factors) scaled arbitrarily by scale. Used to account for static disorder of exciton transition energy
dE_n = energy_disorder  # from .params file

# keep only factors that increase E_mon, since Frenkel exciton absorption spectrum is not Gaussian but favors high energy side. This could be explained by 
# a maximum delocalisation length resulting in the lowest energy bright state beyond which no longer delocalisation is possible. Shorter delocalisation lengths,
# due to structural breaks in the exciton or other effecsts, lead to the lowest/bright transition being at higher energies
# idea1
#dE_n[dE_n.<=1] = 2 .- dE_n[dE_n.<=1] #BUG #CHECK This shifts polariton spectrum to higher E
# idea2, better since keeps all values
#dE_n = 1 .+ abs.(dE_n .- 1)

# can deactive random energy shifts
dE_n = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] .+ 1

# functions/methods to build matter H and transition operators, by embedding monomer H and σ at site N
H_n(n)  = embed(B,n,H_mon * dE_n[n]) # Hamiltonian at site n => Hₙ with random, normally distributed energy offset dE
Σ⁺_n(n) = embed(B,n,σ⁺)              # raising operator at site n 
Σ⁻_n(n) = embed(B,n,σ⁻)              # lowering operator at site n

# need this workaroud to deal with a single TLS... #BUG
if N != 1
    Σ⁺   = sum(Σ⁺_n.(1:N))
    Σ⁻   = sum(Σ⁻_n.(1:N))
    Σ⁺n  = Σ⁺_n.(1:N)
    Σ⁻n  = Σ⁻_n.(1:N)
    Hexc = sum(H_n.(1:N))                   # Hamiltonian of disorder TLSs ... NOT coupled to each other
else
    Σ⁺   = σ⁺
    Σ⁻   = σ⁻
    Σ⁺n   = σ⁺
    Σ⁻n   = σ⁻
    Hexc = H_mon                             # Hamiltonian of disorder TLSs ... NOT coupled to each other
end

logger("\nMatter Hamiltonian:\n\n", dense(Hexc), fn_log)

# create cavity basis
bcav = FockBasis(N_photonStates)  #FACT: #CHECK need second FockState for ESA !?
#bcav = FockBasis(1)
#BUG#TODO: Having 2 photon states messes up transforming between eigen and site basis in subspace
#IDEA MAYBE I need two 1 photon Fock states ? 
Ecav = E_cavity         # slightly red-tuned from excitonic transition, as main system in paper
a    = destroy(bcav) # ⊗ one(bcav) + one(bcav) ⊗ destroy(bcav)
at   = create(bcav)  #⊗ one(bcav) + one(bcav) ⊗ create(bcav)

# embed operators into full basis
A  = a ⊗ one(Hexc)
At = at ⊗ one(Hexc)
Σ⁺ = one(a) ⊗ Σ⁺
Σ⁻ = one(a) ⊗ Σ⁻
Σ⁺n = [one(a) ⊗ Σ⁺n[i] for i in 1:N]
Σ⁻n = [one(a) ⊗ Σ⁻n[i] for i in 1:N]
 
# create cavity H
Hcav = Ecav * at * a

logger("\nCavity Hamiltonian:\n\n", dense(Hcav), fn_log)

# create basis of complete system
Bfull = bcav ⊗ B

# create Jaynes/Tavis-Cummings Hamiltonian for light-matter coupling g 
g = 0.15 / sqrt(N)                               # should make splliting independent of N
#g = 0
H = embed(Bfull,1,Hcav) + one(Hcav) ⊗ Hexc + g * (Σ⁺ * A + Σ⁻ * At)

logger("\nFull Hamiltonian:\n\n", dense(H), fn_log)

# sort by increasing values of Hamiltonian (rough #CHECK ) ... sort into excitation sectors


idx = sortperm(real(diag((H).data)))
H.data = H.data[idx,idx]
A.data = A.data[idx,idx]
At.data = At.data[idx,idx]
Σ⁺.data = Σ⁺.data[idx,idx]
Σ⁻.data = Σ⁻.data[idx,idx]
for i in 1:N
    Σ⁺n[i].data = Σ⁺n[i].data[idx,idx]
    Σ⁻n[i].data = Σ⁻n[i].data[idx,idx]
end

logger("\nOrdered Hamiltonian:\n\n", dense(H), fn_log)

# transition operator
#μ = Σ⁺ + Σ⁻    # transitions between excitonic states only, #FACT: not possible induced by ext. EM field in cavity, every via optical mode
μ = At + A #+ At^2 + A^2                                     #FACT: transitions induced by external field 

# density matrices
rho0 = dm(fockstate(bcav,0) ⊗ tensor([nlevelstate(b,1) for i in 1:N]))
rho1 = dm(fockstate(bcav,0) ⊗ tensor([nlevelstate(b,2) for i in 1:N]))

# transition operator from ground to first ES
μ12 = rho0 * μ + μ * rho0

# sort (see above)


μ12.data  = μ12.data[idx,idx]
rho0.data = rho0.data[idx,idx]
rho1.data = rho1.data[idx,idx]

# transition from first ES to higher ES #FACT: should connect excitation on exciton site and excitation of exciton + cavity
μ23 = μ - μ12

# normalization of μ12 and μ23 ... seems not really necessary #DELETE
rho1 = μ12 * rho0 * μ12
#μ12 = μ12 / sqrt(tr(rho1))
#rho1 = μ12 * rho0 * μ12
rho2 = μ23 * rho1 * μ23
#μ23 = μ23 / sqrt(tr(rho2))
#μ23 = μ23 
rho2 = μ23 * rho1 * μ23

# Lindblad dissipation operators
L1 = A           # FACT: #CHECK only the relaxation channel from single exc. sector to GS is relevant ...
L2 = Σ⁻          # FACT: actually, having the elements that connect single and double exc. sector seems to ...
                        # FACT: broaden the ESA peaks (the only signal that considers single and double exc. sectors)                        
L3 = .5 * ((At*A) - (A*At)) # not sure if this one is relevant/meaningful/... #CHECK

#L3 = A  - L1
#L4 = Σ⁻ - L2

# rates and list of Lindblad dissipation operators
Γ =      [.25, .000001]         # Γ[1] was 0.2 before
L = Γ .* [L1,  L2]

# plot energy levels
figure(); title("Energy level diagram")
plot_levels(H,0)
plot_levels(Hcav,-1)
plot_levels(Hexc,1)
plot_levels(H_mon,2)
xticks([-1, 0, 1, 2], ["cavity", "full", "matter", "monomer"])

##

# create subspace (single excitation, double excitation, ...) in eigenbasis/polariton basis
H_si, transf_op_si, P_si, L_si, rho0_si, μ12_si, μ23_si, Σ⁺_si, Σ⁻_si, At_si, A_si = create_subspace([H],"si", L, rho0, μ12, μ23, Σ⁺, Σ⁻, At, A)
H_si = H_si[1]
H, transf_op, P, L, rho0, μ12, μ23, Σ⁺, Σ⁻, At, A = create_subspace([H],"bi", L, rho0, μ12, μ23, Σ⁺, Σ⁻, At, A)
H = H[1]
#FACT: Using "bi_lowest" does not give ESA that eliminates above diagonal cross-peak
#FACT: Neither does using "bi_polariton", which takes lowest and highest levels of double excitation sector
#FACT: only using "bi" adds ESA signal at above diagonal cross-peak that cancels out the GSB above diagonal cross-peak after some time T



Σ⁺n_si = [P_si * Σ⁺n[i] * P_si' for i in 1:N]
Σ⁻n_si = [P_si * Σ⁻n[i] * P_si' for i in 1:N]

Σ⁺n = [P * Σ⁺n[i] * P' for i in 1:N]
Σ⁻n = [P * Σ⁻n[i] * P' for i in 1:N]

logger("\nCHANGE TO SUBSPACE", fn_log)
logger("\n0th, 1st, 2nd Manifold", fn_log)
logger("\nCoupled (eigen) basis", fn_log)
logger("\nHamiltonian:\n\n", dense(H), fn_log)
logger("\nTransition operator gs <-> es:\n\n", dense(μ12), fn_log)
logger("\nTransition operator es <-> d-es:\n\n", dense(μ23), fn_log)
logger("\n0th, 1st, 2nd Manifold", fn_log)
logger("\nUncoupled (site) basis", fn_log)
logger("\nHamiltonian:\n\n", transf_op * H * transf_op', fn_log)
logger("\nTransition operator gs <-> es:\n\n", transf_op * μ12 * transf_op', fn_log)
logger("\nTransition operator es <-> d-es:\n\n", transf_op * μ23 * transf_op', fn_log)
logger("\n0th, 1st Manifold", fn_log)
logger("\nCoupled (eigen) basis", fn_log)
logger("\nHamiltonian:\n\n", dense(H_si), fn_log)
logger("\nTransition operator gs <-> es:\n\n", μ12_si, fn_log)
logger("\nTransition operator es <-> d-es:\n\n", μ23_si, fn_log)
logger("\nUncoupled (site) basis", fn_log)
logger("\nHamiltonian (o, 1, 2 Manifold; site):\n\n", transf_op_si * H_si * transf_op_si', fn_log)
logger("\nTransition operator gs <-> es:\n\n", transf_op_si * μ12_si * transf_op_si', fn_log)
logger("\nTransition operator es <-> d-es:\n\n", transf_op_si * μ23_si * transf_op_si', fn_log)
plot_levels(H_si,-0.05,col="r",ls="dashed")
plot_levels(H   , 0.05,col="b",ls="dotted")


# transf_op can be used to convert between polariton and site basis #FACT: required prior sorting of operators into blocks of excitation sectors
# H_site    = transf_op     * H     * transf_op'
# H_site_si = transf_op_si  * H_si  * transf_op_si'

#TODO: vary g instead of E_mon ... try to see effect
g = g .* dE_n
#Hk(i) = (Ecav * At_si * A_si) + (Emon * Σ⁺_si * Σ⁻_si) + g[1] * (Σ⁺_si * A_si + Σ⁻_si * At_si)

# time list
tlist = [0:1.1:200;]

logger("\nTime list:\n\n", tlist, fn_log)

# Lindblad time evolution
corr    = timecorrelations.correlation(tlist, rho0_si, H_si, L_si, μ12_si, μ12_si)
# #TODO: test with varying g ... 
#corr = sum([timecorrelations.correlation(tlist, rho0_si, Hk(i), L_si, μ12_si, μ12_si) for i in 1:1]) ./ 1

# zeropadding for smoother data only
zp = 10
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

# calculate and plot spectrum
ω, spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true);
figure()
subplot(211)
plot(ω,spec)

# load and plot experimental absorption spectrum
TDBC_cav   = readdlm("TDBC_cav.dat")
TDBC_cav_E = TDBC_cav[:,2] ./ 8065      # convert from cm^-1 to eV
TDBC_cav_I = TDBC_cav[:,3] ./ 0.176     # normalize to 1 (in the visible)

plot(-TDBC_cav_E,TDBC_cav_I,"r:")

# TODO: spectral density ... so far only trial and error ... Seems like it needs to bridge the E gap between P+ and P- for relaxation to occur
width = .005
#ωq_TDBC     = [0.12 0.336]
#ampl_TDBC   = [0.12 3] / (N)

# from ... http://www.rsc.org/suppdata/c8/sc/c8sc00171e/c8sc00171e1.pdf
#ωq_TDBC = [40, 80, 120, 150, 185, 197] ./ 8065 #.* 0.02998 # cm-1 in Hz #./ 8065 # in eV instead
# => [0.00496 0.00992 0.01488 0.01860  0.02294 0.02443]
#ampl_TDBC = [14, 18,  25,  43,  42, 67, 60]  

#ωq_TDBC     = [0.0084, 0.0120, 0.1204, 0.1289, 0.1304, 0.1424, 0.303] #FACT: relaxation from UP to LP when ωq_TDBC bridges energy gap between UP and LP. However, no relaxation to dark states ... DUE TO LACK OF COUPLING ????
ωq_TDBC     = ωq        # from .params file
ampl_TDBC   = Aq        # from .params file


logger("Spectral frequencies: ", ωq_TDBC, fn_log)
logger("Spectral amplitudes: ", ampl_TDBC, fn_log)
logger("Spetral peak width: ", width, fn_log)


try     # need to delete method first when working on it
    local m = @which spectral_density(1)
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
        Jω = .0
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
    return Jω

end

#plot spectral density 
#figure()
subplot(212)
plot(collect(-1:.001:1),(spectral_density.(collect(-1:.001:1))))

# Rate and operator for redfield process (GS -> single exc sector)
Γ_R = .23

# this means collective mode of all TLSs dephases ...
a_ops_si = [Γ_R * Σ⁺_si * Σ⁻_si, spectral_density]
# ... and this: each TLS dephases individually
a_ops_si = vcat([[Σ⁺n_si[i] * Σ⁻n_si[i], spectral_density] for i in 1:N]...)


R_si, ekets = timeevolution.bloch_redfield_tensor(H_si, [a_ops_si]; J=1 .* L_si)

# time evolution, Redfield
tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12_si*rho0_si,R_si,H_si)

corr = expect(μ12_si,rhot)
# calculate and plot spectrum
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
ω, spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true);
subplot(211)
plot(ω,spec)
#FACT: J(0) > 0 induces pure dephasing and leads to broadening of ALL transitions
#FACT: J(ω) > 0 for ω > 0 leads to relaxation between states with Ei - Ef = ω ... 
#CHECK what role does the operator play ? 

figure(figsize=(14,4));
subplot(131)
plot(tout,real(corr[1:length(tout)]))
tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12_si*rho0_si*μ12_si,R_si,H_si)
rho1test = copy(rho0_si)
rho1test.data[1,1] = 0
rho1test.data[2,2] = 1
#tout, rhot = timeevolution.master_bloch_redfield(tlist,rho1test,R_si,H_si)
#corr_gs = expect(Σ⁺_si*Σ⁻_si,rhot)
corr_gs = expect(rho0_si,rhot)
#corr_es = expect(Σ⁻_si*Σ⁺_si,rhot)
#corr_es = expect(μ12_si * rho0_si * μ12_si,rhot)
corr_es_mat = expect(Σ⁺_si * rho0_si * Σ⁻_si / N,rhot)
corr_es_cav = expect(At_si * rho0_si * A_si, rhot)
plot(tout,real(corr_gs))
#plot(tout,corr_es) # == corr_es_cav
plot(tout,real(corr_es_cav))
plot(tout,real(corr_es_mat))
legend(["corr. function", "ground state", "cavity", "matter"])

rhot_site = [transf_op_si * i  * transf_op_si' for i in rhot]
corr_es1_mat = expect(rho0_si, rhot_site)
#plot(tout,corr_es1_mat)

vecs   = eigvecs(dense(H_si).data)
eivecs = [Ket(H_si.basis_l,vecs[:,i]) for i in 1:length(vecs[1,:])]

subplot(132);
title("eigen/coupled basis") # levels: 1 -> GS, 2 -> LP, 3 - N-1 -> ..., N -> UP
for i in 1:H_si.basis_l.shape[1]
    local es_n = real(expect(dm(eivecs[i]),rhot))
    plot(tout,es_n,linestyle="dashed", label=string(i)); #append!(leg,[string(i)])
    ylim((-0.1, 1))
end
legend()

es_mat = zeros(length(tout))
subplot(133); title("uncoupled basis") # levels: 1 -> GS, 2 -> cav., 3 - N -> TLSs
for i in 1:H_si.basis_l.shape[1]
    local es_n = real(expect(dm(eivecs[i]),rhot_site))
    plot(tout,es_n,linestyle="dashed", label=string(i)); #append!(leg,[string(i)])
    if i > 2
        global es_mat = es_mat .+ es_n
    end
end
plot(tout, es_mat, label="Σ(TLSs)")
legend()


# Rate and operator for redfield process (GS up to double exc sector)
a_ops = [Γ_R * Σ⁺ * Σ⁻, spectral_density]
a_ops = vcat([[Σ⁺n[i] * Σ⁻n[i], spectral_density] for i in 1:N]...)

R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops]; J=1 .* L)

# choose method to use for simulation of 2D spectra
method = "redfield"

# change parameters accordingly #TODO: simplify a bit
if method == "lindblad"
    F = L
    F_si = L_si
    use_sub = true
    H = H[1]
elseif method == "redfield"
    F = R
    F_si = R_si
    use_sub = true
end

if calc_2d
    ## calculate (complex) 3rd order corr function (with T=0)
    zp = 11 # zeropad up to 2^zp

    ## calculate 2D spectra at
    T = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 100, 200, 300] # energy eV -> t => [hbar / eV] ~ 0.66 fs
    T = [0, 30]
    spectra2d = Array{out2d}(undef, length(T))

    # cannot plot inside cmds.jl when multithreading
    #for i = 1:length(T)
    # multithreading (run several T steps in parallel)!
    Threads.@threads for i = 1:length(T)
        spectra2d[i] = make2Dspectra(tlist,[rho0, rho0_si],[H, H_si],[F, F_si],[μ12, μ12_si],[μ23, μ23_si],T[i],
                                            method;debug=true,use_sub=use_sub,
                                                t2coh="kin",zp=zp)
    end

    ## crop 2D data and increase dw
    spectra2d = crop2d.(spectra2d,1.5;w_max=2.6,step=1)
    spectra2d = round2d.(spectra2d)
end

if calc_2d
    ## plot 2D spectra for each(?) T
    # what to plot
    rep="absl"
    scal="lin"

    ## make  subplot layout
    nplots = length(T);
    ncols = Int32(ceil(sqrt(nplots)));
    nrows = Int32(ceil(nplots / ncols));

    # determine maximum value in dataset out2d[:].full2d[:,:]
    maxi = Float64(maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)]))

    # plot 2D spectra
    fig, ax = subplots(nrows,ncols,sharex=true,sharey=true,figsize=(ncols*3.5,nrows*3))
        if nplots == 1
                ax = [ax]
        end
    fig.suptitle(rep * " 2D spectrum (" * scal * ". scaling)")
        k = 0
        for i = 1:nrows
            for j = 1:ncols
                global k += 1
                if k > nplots
                    continue
                end
                sca(ax[j,i])
                ax[j,i].set_aspect="equal"
                plot2d(spectra2d[k].ω,round.(spectra2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi)
                title("2D spectrum at $(T[k]) fs")
            end
        end
        tight_layout()
        subplots_adjust(top=0.825)
  
    ## plot TA (summed 2D spectrum)
    figure()
    ta = [sum(spectra2d[i].full2d,dims=1) for i in 1:length(spectra2d)]
    #[plot(out2d[1].ω,ta[i]') for i in 1:length(ta)]
    [plot(3e8 * 4.136e-15 ./ spectra2d[1].ω ./ 1e-9,-ta[i]') for i in 1:length(ta)]
    xlabel("Detection (ω₃)")
    legend([string(i) * " fs" for i in T])
    tight_layout()
end


# to plot time traces 
# cmds.plot_timeTrace([real(i.full2d) for i in out2d],T,out2d[1].ω,[1.925, 1.925, 2.24, 2.24],[1.925, 2.24, 2.24, 1.925])

