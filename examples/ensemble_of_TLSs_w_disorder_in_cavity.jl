using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles, Random, Combinatorics, Distributions

# make sure to set script directory as pwd()
cd(@__DIR__)

# plot in gui instead of side panel
pygui(true)

# set various options
calc_2d = true                                          # might be costly
use_sub = true                                          # can use single and double excitation subspace with lindblad master equation

# set colormap from cmds module
cmp = create_colormap("bright");

## start actual problem
# CAVITY MODE PARAMETERS
ωr  = 5                                                 # energy
Na  = 1                                                 # number of cavity fock states
A   = FockBasis(Na)                                     # basis
a   = destroy(A)                                        # annihilation operator
at  = create(A)                                         # creation operator 
H_r = (ωr * at * a)                                     # cavity Hamiltonian

# SINGLE TLS PARAMETERS
ω_TLS  = 5                              # Energy
N_TLS  = 2                              # number of states (GS, 1st ES, ...)
b_TLS  = NLevelBasis(N_TLS)             # basis
sm     = transition(b_TLS, 1, 2)        # σ- transition 2-->1 ; lowering operator
sp     = transition(b_TLS, 2, 1)        # σ+ transition 2<--1 ; raising operator
sz     = sp * sm                        # number operator
H_TLS = (ω_TLS * sz)                    # Hamiltonian for single TLS

# theory background for embedding operation: https://docs.qojulia.org/tutorial/ search embed()

# set number of two-level systems to be used
num_of_TLSs = 2

# generate composite basis of TLSs
B_TLSs = b_TLS^num_of_TLSs               # creates a composite basis (alt.: CompositeBasis(repeat([b_TLS],num_of_TLSs)...))

# generate basis of full system
B_full = B_TLSs ⊗ A

# #CHECK ???  !! NEED TO GET RID OF ELEMENTS WITH HIGHER EX SECTORS
cccc  = diagonaloperator(b_TLS,[1,1/2])

# define σ- and σ+ in composite basis #TODO call them Σ for consistency
#Sm = +([tensor(circshift(append!([sm],repeat([one(sm)],num_of_TLSs-1)),i-1)...) for i in 1:num_of_TLSs]...)
Sm = +([tensor(circshift(append!([sm],repeat([cccc],num_of_TLSs-1)),i-1)...) for i in 1:num_of_TLSs]...)
Sp = dagger(Sm)

# TODO figure out how to construct j_ops correctly for different situations
# jump operators for Lindblad time evolution
# j_ops needs to be "list" of objects for relaxation from the individual excited states

# TODO what's this useful for / what is correct ? 
temp = one(sm)
temp.data[2,2] = 1/sqrt(2)
#j_ops = [tensor(circshift(append!([sm],repeat([one(sm)],num_of_TLSs-1)),i-1)...) for i in 1:num_of_TLSs]
j_ops = 1 * [tensor(circshift(append!([sm],repeat([temp],num_of_TLSs-1)),i-1)...) for i in 1:num_of_TLSs]
#j_ops = [tensor(circshift(append!([sm],[j*one(sm) for j in 1:num_of_TLSs-1 ]),i-1)...) for i in 1:num_of_TLSs]

j_ops_b = 1. * [tensor(circshift(append!([sp,sm],repeat([temp],num_of_TLSs-2)),i-1)...) for i in 1:num_of_TLSs]
j_ops_c = 1 * [tensor(circshift(append!([sm,sp],repeat([temp],num_of_TLSs-2)),i-1)...) for i in 1:num_of_TLSs]

pop!(j_ops_b)
pop!(j_ops_c)

append!(j_ops,j_ops_b)
#append!(j_ops,j_ops_c)

#j_ops = [j_ops[i] + dagger(j_ops[i]) for i in 1:length(j_ops)]

## define transition dipole operators
# transitions from GS (should only go to cavity mode #CHECK)
μ12   = +([tensor(circshift(push!([sm*sp],repeat([sm*sp],num_of_TLSs-1)...),i-1)...) for i in 1:num_of_TLSs]...)
#μ12   = tensor(circshift(push!([sm],repeat([sm*sp],num_of_TLSs-1)...),1-1)...)

## add Fock space
#μ12 = μ12 ⊗ one(at) + one(μ12) ⊗ (at + a)
μ12 = μ12 ⊗ (at + a)

# transitions from single excited state to double excited state (energy reside in matter part = cavity can be excited again #CHECK )
elements = collect(permutations(push!([sp*sm],repeat([one(b_TLS)],num_of_TLSs-1)...)))
μ23 = +([tensor(elements[i]...) for i in 1:length(elements)]...)

# just one ! 
#μ23 = tensor(one(b_TLS) , sp*sm )
#μ23 = (μ23 + dagger(μ23))
# normalize later
#μ23 = μ23 ⊗ one(at) + one(μ23) ⊗ (at + a)
μ23 = μ23 ⊗ (at + a)

## create static disorder by normally distributing the site energies
nvals = 1000
disorder = sort(rand(Normal(0,.035),nvals+1))[1:Int32(floor(nvals/(num_of_TLSs-1))):nvals+1]
disorder[1] = sort(rand(Normal(0,.035),nvals+1))[50]
disorder[end] = sort(rand(Normal(0,.035),nvals+1))[end-50]
# scale disorder
disorder ./= 2
# + 1, because I multiply with this values rather than add it to the energy 
disorder .+= 1

# turn disorder off
disorder .= 1

## create diagonal, composite Hamiltonian from TLS Hamiltonian and include disorder of energy levels
# no coupling, isolated TLSs
H_en = +([disorder[i] * embed(B_TLSs,i,H_TLS) for i in 1:num_of_TLSs]...);

## #IDEA: can possibly consider nearest-neighbour interactions for Hamiltonian
#H_int = +([tensor(circshift(push!([sm, sp], repeat([one(sm)],num_of_TLSs-2)...),i-1)...) for i in 1:num_of_TLSs]...)
#H_int = H_int + dagger(H_int)

## construct full composite Hamiltonian
g0 = 0.1 
g = g0 * sqrt(num_of_TLSs)
H = H_en ⊗ one(H_r) + one(H_en) ⊗ H_r + g * (Sm ⊗ at + Sp ⊗ a)

## get energies and eigenvectors
E      = eigvals(dense(H).data)
eivecs = eigvecs(dense(H).data)
eivecs = [Ket(B_full,eivecs[:,i]) for i in 1:length(E)]
Psi0   = eivecs[1]
rho0   = dm(Psi0)

## print eigenenergies
xn = 0
figure(figsize=(5,4));
title("Energy levels (Eigenvalues) of H_FH")
#hlines([0, ωr], -.65, -.35)
hlines([0],  .35,  .65)
hlines([ω_TLS .* disorder[:]],  .35,  .65)
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


## normalize μ12, such that tr(rho1) = 1
μ12   = μ12 / binomial(num_of_TLSs,1)
# excite density matrix to normalize μ23
rho1 = μ12 * (rho0 * μ12)

#TODO: correctly normalize μ23
rho2 =  μ23 *  rho1 * μ23
# μ23 = μ23 / binomial(num_of_TLSs,0) ### somehow in this way ??? 

## j_ops in full system
Γ     = [.3] #.1 * [.1,.6]
j_ops = Γ .* [tensor(j_ops[i],one(a)) for i in 1:3]


# #DEBUG: view rho1 after "excitation"
figure(figsize=(4,4));
#imshow(real(rho1.data),origin="lower",cmap="Blues")
pcolormesh(real(rho1.data),edgecolor="k",cmap="Blues",linewidth=".5")

## time list
tlist = [0:0.3:150;]

## master equation time evolution starting after population transfer into excited state (cavity mode)
tout, rhot = timeevolution.master(tlist,rho1,H,j_ops);#;rates=Γ);
view_dm_evo(rhot,5)
# ground and excited state population, time evolution
n_gs = real(expect(rho0, rhot))
n_es = real(expect(rho1, rhot))

## plot population in ground and excited state
leg = []                          # legend
figure(figsize=(10,4));
subplot(121)
plot(tout,n_gs); leg = ["ngs"]
plot(tout,n_es); append!(leg,["nes"])
legend(leg)
leg = []
# also plot partial population of each TLS and the cavity
subplot(122)
for i in 1:length(E) #TODO: distinguish between eigenstates
    local n_es = real(expect(dm(eivecs[i]),rhot))
    plot(tout,n_es); append!(leg,[string(i)])
end
legend(leg)

#DELETE
rhotest = tri(tri(rho1,"U"),"L")
rhotest = rho1
rhotest.data .= 0
rhotest.data[2,2] = 1

## calculate correlation function
corr    = timecorrelations.correlation(tlist, rho0, H, j_ops, μ12, μ12)             # for ground state absorption
#corr    = timecorrelations.correlation(tlist, rho1, H, j_ops, μ23, μ23);           # for excited state (rho1) absorption #CHECK: shouldn't work like this yet

## calculate spectrum from correlation function
ω, spec = timecorrelations.correlation2spectrum(tlist, corr; normalize_spec=true);

## plot correlation function and spectrum
figure(figsize=(10,4));
ax1 = subplot(211);
plot(tlist,corr)
ax2 = subplot(212);
plot(ω, spec)

## create Redfield tensor #CHECK #TODO
ωC = 1; S = 0.1; γ1 = 10
function noise_power(ω)
    if ω == 0.0
        return 0
    else
        return 1 / π * γ1 * ωC ./ (ω.^2 .+ ωC^2)
    end
end
test = Sm*Sp ⊗ one(a)
test.data = tril(triu(test.data))
a_ops = [sqrt(.5)*test, noise_power]
R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops], J=j_ops)

## select methods for time evolution during 2D calculation
method = "lindblad"                                                     # working
#method = "redfield"                                                    # not working

## change parameters
if method == "lindblad"
    F = j_ops
elseif method == "redfield"
    F = R
end

## calculate (complex) 3rd order corr function (with T=0)
if calc_2d
    
    zp = 10 # zeropad up to 2^zp

    ## calculate 2D spectra at
    T = [5] #fs

    # init output array
    spectra2d = Array{out2d}(undef, length(T))

    # cannot plot inside cmds.jl when multithreading
    #for i = 1:length(T)
    # multithreading (run several T steps in parallel)!
    Threads.@threads for i = 1:length(T)
        spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                            method;debug=true,use_sub=false, #TODO! use_sub true and false give different results. For false SE and ESA have sign changes, true looks better
                                                zp=zp)
    end

    ## crop 2D data and increase dw
    spectra2d = [crop2d(spectra2d[i],4;w_max=6,step=1) for i = 1:length(T)]

    ## plot 2D spectra for each(?) T
    # what to plot
    rep="absorptive"
    scal="lin"

    ## make  subplot layout
    nplots = length(T);
    ncols = Int32(ceil(sqrt(nplots)));
    nrows = Int32(ceil(nplots / ncols));

    ## determine maximum value in dataset out2d[:].full2d[:,:]
    maxi = maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)])

    ## plot 2D spectra
    fig, ax = subplots(nrows,ncols,sharex=true,sharey=true,figsize=(ncols*3.2,nrows*3))
    if length(T) == 1
        ax = [ax]
    end
    fig.suptitle(rep * " 2D spectrum (" * scal * ". scaling)")
    k = 0
    for i = 1:ncols
        for j = 1:nrows
            global k += 1
            sca(ax[i,j])
            if k > nplots
                continue
            end
            ax[i,j].set_aspect="equal"
            plot2d(spectra2d[k].ω,round.(spectra2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi)
            title("2D spectrum at $(T[k]) fs")
        end
    end
    tight_layout()
    subplots_adjust(top=0.8)

    ## plot TA (summed 2D spectrum) #TODO convolute with narrow (pump) laser spectrum
    figure()
    ta = [sum(real(spectra2d[i].full2d),dims=1) for i in 1:length(spectra2d)]
    [plot(spectra2d[1].ω,ta[i]') for i in 1:length(ta)]

end



## #DEBUG How to check execution time (run twice after function assignment!):

#=
function test_timing()
    @time H = +([disorder[i] * embed(B_TLSs,i,H_TLS) for i in 1:num_of_TLSs]...);  # works with +, -, *, /
    #@time H = last(accumulate(+, [disorder[i] * embed(B_TLSs,i,H_TLS) for i in 1:num_of_TLSs]));  # works with +, -, *, /
    return H
end
test_timing()
=#

