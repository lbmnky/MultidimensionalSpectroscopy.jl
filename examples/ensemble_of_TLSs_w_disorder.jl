using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles, Random, Combinatorics, Distributions

# make sure to set script directory as pwd()
cd(@__DIR__)

# plot in gui instead of side panel
pygui(true)

# set different options
calc_2d = true
use_sub = true

# set colormap from cmds module
cmp = create_colormap("bright");

## start actual problem
# CAVITY MODE (not used so far)
ωr  = 4
Na  = 2               # number of cavity fock states
A   = FockBasis(Na)
a   = destroy(A)
at  = create(A)
H_r = (ωr * at * a)

# set number of two-level systems to be used
num_of_TLSs = 3
# SINGLE TLS parameters
ω_TLS  = 5                              # Energy
N_TLS  = 2                              # number of states (GS, 1st ES, ...)
b_TLS  = NLevelBasis(N_TLS)             # basis
sm     = transition(b_TLS, 1, 2)        # σ- transition 2-->1
sp     = transition(b_TLS, 2, 1)        # σ+ transition 2<--1
sz     = sp * sm                        # "number" operator
H_TLS = (ω_TLS * sz)                    # Hamiltonian for single TLS

# background for embedding operation: https://docs.qojulia.org/tutorial/ search embed()

# generate composite basis
B_TLSs = b_TLS^num_of_TLSs # creates a composite basis (alt.: CompositeBasis(repeat([b_TLS],num_of_TLSs)...))

# σ- and σ+ in composite basis
# Sm = sm ⊗ one(sm) ⊗ one(sm) ⊗ ... + one(sm) ⊗ sm ⊗ one(sm) ⊗ ... + one(sm) ⊗ one(sm) ⊗ sm ⊗ ... + ...
Sm = +([tensor(circshift(append!([sm],repeat([one(sm)],num_of_TLSs-1)),i-1)...) for i in 1:num_of_TLSs]...)
Sp = dagger(Sm)

# TODO figure out how to construct j_ops correctly for different situations
# jump operators for Lindblad time evolution
# j_ops needs to be "list" of objects for relaxation from the individual excited states
temp = one(sm)
#temp.data[2,2] = 1/sqrt(2)

# Lindblad operators for relaxation of each TLS to ground state
# Lᵢ = Σ{j=1}{N²-1} uⱼᵢAⱼ ... what's what ?
# N-dimensional system: N=2 molecules -> 4-1=3
L = .1 * [tensor(circshift(append!([sm],repeat([temp],num_of_TLSs-1)),i-1)...) for i in 1:num_of_TLSs]
# L[1] = [sm ⊗ temp ⊗ temp ⊗ ...]
# L[2] = [temp ⊗ sm ⊗ temp ⊗ ...]
# L[3] = ...

# coupling/decay between TLS i and i+1
# K[1] = [sm ⊗ sp ⊗ temp ⊗ ...]
# K[2] = [temp ⊗ sm ⊗ sp ⊗ ...]
K = .4 * [tensor(circshift(append!([sm,sp],repeat([one(temp)],num_of_TLSs-2)),i-1)...) for i in 1:num_of_TLSs+1]

# coupling between TLS i and i+2
#append!(L,K)
#j_ops = [tensor(circshift(append!([sm],[j*one(sm) for j in 1:num_of_TLSs-1 ]),i-1)...) for i in 1:num_of_TLSs]

#j_ops_b = 1. * [tensor(circshift(append!([sp,sm],repeat([temp],num_of_TLSs-2)),i-1)...) for i in 1:num_of_TLSs]
#j_ops_c = 1. * [tensor(circshift(append!([sm,sp],repeat([temp],num_of_TLSs-2)),i-1)...) for i in 1:num_of_TLSs]
#j_ops_d = 1. * [tensor(circshift(append!([sp,temp,sm],repeat([temp],num_of_TLSs-3)),i-1)...) for i in 1:num_of_TLSs]
#j_ops_e = 1. * [tensor(circshift(append!([sm,temp,sp],repeat([temp],num_of_TLSs-3)),i-1)...) for i in 1:num_of_TLSs]
#pop!(j_ops_b)
#pop!(j_ops_c)


#append!(j_ops,j_ops_b)
#append!(j_ops,j_ops_c)

#j_ops = [j_ops[i] + dagger(j_ops[i]) for i in 1:length(j_ops)]

# transition dipole operators
# transitions from GS
# μ12 = sm ⊗ (sm*sp) ⊗ (sm*sp) ⊗ ... + (sm*sp) ⊗ sm ⊗ (sm*sp) ⊗ ... + ...
# μ12 = Σᵢ σ⁻ᵢ ??????
μ12   = +([tensor(circshift(push!([sm],repeat([sm*sp],num_of_TLSs-1)...),i-1)...) for i in 1:num_of_TLSs]...)
# normalize μ12, such that tr(rho1) = 1
μ12   = (μ12 + dagger(μ12)) / sqrt(num_of_TLSs)

# transitions from single excited state to double excited state
elements = collect(permutations(push!([sm,sp*sm],repeat([sm*sp],num_of_TLSs-2)...)))
μ23 = +([tensor(elements[i]...) for i in 1:length(elements)]...)
μ23 = (μ23 + dagger(μ23))
# normalize later

# create static disorder by normally distributing the site energies
nvals = 1000
disorder = sort(rand(Normal(0,.035),nvals+1))[1:Int32(floor(nvals/(num_of_TLSs-1))):nvals+1]
disorder[1] = sort(rand(Normal(0,.035),nvals+1))[50]
disorder[end] = sort(rand(Normal(0,.035),nvals+1))[end-50]

disorder .+= 1

# create diagonal composite Hamiltonian from TLS Hamiltonian
H_en = +([disorder[i] * embed(B_TLSs,i,H_TLS) for i in 1:num_of_TLSs]...);

# consider nearest-neighbours only for interaction Hamiltonian
H_int = +([tensor(circshift(push!([sm, sp], repeat([one(sm)],num_of_TLSs-2)...),i-1)...) for i in 1:num_of_TLSs]...)
H_int = H_int + dagger(H_int)

# construct full composite Hamiltonian
H = H_en + .0 * H_int

# test to generate subspace
eivecs      = eigvecs(dense(H).data)
n_sub       = 1 + length(H.basis_l.bases) + binomial(length(H.basis_l.bases),2)
eivecs_sub  = [Ket(B_TLSs,eivecs[:,i]) for i in 1:n_sub]
Psi0 = eivecs_sub[1]
rho0 = dm(Psi0)
b_sub = SubspaceBasis(H.basis_l,eivecs_sub)
println("dim(Hilbert space): ", length(B_TLSs))
P = sparse(projector(b_sub, H.basis_l))
Pt = dagger(P)
smm = embed(B_TLSs,1,sm)
"""
if use_sub
    println("dim(subspace): ", length(b_sub))
    H = P * H * Pt
    μ12 = P * μ12 * Pt
    μ23 = P * μ23 * Pt
    rho0 = P * rho0 * Pt
    Sm = P * Sm * Pt
    Sp = P * Sp * Pt
    smm = P * smm * Pt
    j_ops = [P * j_ops[i] * Pt for i in 1:num_of_TLSs]
end
"""

# get energies and eigenvectors
E      = eigvals(dense(H).data)
vecs   = eigvecs(dense(H).data)
eivecs = [Ket(B_TLSs,vecs[:,i]) for i in 1:n_sub]

# print eigenenergies
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

# excite density matrix to normalize μ23
rho1 = (μ12 * (rho0 * μ12))
rho2 =  μ23 *  rho1 * μ23
μ23  =  sqrt(2) * μ23 / sqrt(tr(rho2))
rho2 =  μ23 *  rho1 * μ23

# debug: view rho1 after "excitation"
figure(figsize=(4,4));
#imshow(real(rho1.data),origin="lower",cmap="Blues")
pcolormesh(real(rho1.data),edgecolor="k",cmap="Blues",linewidth=".5")
tlist = [0:0.3:150;]

Γ     = [.3] #.1 * [.1,.6]
#L = Γ .* L

# master equation time evolution
tout, rhot = timeevolution.master(tlist,rho1,H,L);#;rates=Γ);

view_dm_evo(rhot,5)

n_gs = real(expect(rho0, rhot))
n_es = real(expect(rho1, rhot))

leg = []
# plot population
figure(figsize=(10,4));
subplot(121)
plot(tout,n_gs); leg = ["ngs"]
plot(tout,n_es); append!(leg,["nes"])
legend(leg)
leg = []
subplot(122)
for i in 1:num_of_TLSs
    local n_es = real(expect(dm(eivecs[i]),rhot))
    plot(tout,n_es); append!(leg,[string(i)])
end
legend(leg)

rhotest = tri(tri(rho1,"U"),"L")
rhotest = rho1
rhotest.data .= 0
rhotest.data[2,2] = 1

# calculate correlation function
# for ground state absorption
corr    = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12);
# for excited state (rho1) absorption
#corr    = timecorrelations.correlation(tlist, rho1, H, j_ops, μ23, μ23);

# calculate spectrum from correlation
ω, spec = timecorrelations.correlation2spectrum(tlist, corr; normalize_spec=true);

# plot correlation function and spectrum
figure(figsize=(10,4));
ax1 = subplot(211);
plot(tlist,corr)
ax2 = subplot(212);
plot(ω, spec)

# Redfield

τc = 60
Λ = 1/τc
λ = 1

#  Ohmic  exponentially  cutoffspectral density
ωc = 0.01860 * 100
λ  = 0.00434 * 100
ω0 = .3
width = .1

G = 1
LL = .2


try
    m = @which spectral_density(1)
    Base.delete_method(m)   
catch
end
function spectral_density(ω)
    if ω == 0.0
        return 0
    else
        w = (ω .- ω0) ./ (.5 * width)
        return abs.(λ * 2 * Λ * ω ./ (ω.^2 .+ Λ^2)) .+ .1 ./ (1 .+ (w).^2)
        #return abs.(λ * 2 * Λ * ω ./ (ω.^2 .+ Λ^2))
        #return .1 ./ (1 .+ w.^2)
        #return abs.(λ .* ω .+ .1 ./ (1 .+ (w).^2))
        #return abs.(λ / ωc * ω .* exp.(-ω ./ ωc))       #  Ohmic  exponentially  cutoffspectral density: λ reorganization energy, ωc cutoff frequency 
        #return abs.(LL .* ω .* G) #./ (G^2 .+ w.^2)
    end
end
test = Sm*Sp
test.data = tril(triu(test.data))
#a_ops = [sqrt(.5)*test, noise_power];

#S = append!(j_ops_b,j_ops_c)

#a_ops = [[10 * j_ops_b[i] * j_ops_b[i]', noise_power] for i in 1:length(j_ops_b)];
#a_ops = []

#a_ops = [1/3 * tensor(circshift(append!([sm*sp],repeat([one(b_TLS)],num_of_TLSs-1)),i)...) for i in 1:num_of_TLSs]
a_ops = [Sp * Sm, spectral_density]

temp2 = []
for i = 1:length(K)
    println(i)
    append!(temp2,[.3 * (K[i] + dagger(K[i]))])
    append!(temp2,[spectral_density])
end

S = temp2

w = collect(0:.01:10)
figure()
plot(w, spectral_density(w))
## Here, the TLSs are not coupled to each other. S describes the energy transfer
#  between TLSs and L describes the decay of each individual TLS to the ground
#  state. Thus S must be an off-diagonal element connecting two diagonal
#  populations (between two TLSs) and L is off-diagonal between each TLS and the
#  ground state.
#R, ekets = timeevolution.bloch_redfield_tensor(H, a_ops, J=0.1 .* append!(j_ops_b,j_ops_c))
R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops], J= .03 .* L)

tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0,R,H);
corr = expect(μ12,rhot);

zp = 11
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tout,zp)

ω,spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)

sca(ax1)
plot(tnew,real(corr),"g",linewidth=1,label="RF")
sca(ax2)
plot(ω,spec,"g",linewidth=1,label="RF")

rho1 = μ12 * rho0 * μ12
tout, rhot = timeevolution.master_bloch_redfield(tlist,rho1,R,H);

n_gs = real(expect(rho0, rhot))
n_es = real(expect(rho1, rhot))

leg = []
# plot population
figure(figsize=(10,4));
subplot(121)
plot(tout,n_gs); leg = ["ngs"]
plot(tout,n_es); append!(leg,["nes"])
legend(leg)
leg = []
subplot(122)
for i in 1:length(eivecs)
    local n_es = real(expect(dm(eivecs[i]),rhot))
    plot(tout,n_es); append!(leg,[string(i)])
    ylim(0, 1)
end
legend(leg)

# select methods for time evolution during 2D calculation
method = "lindblad"
#method = "redfield"

# change parameters
if method == "lindblad"
    F = L
    use_sub = true
elseif method == "redfield"
    F = R
    use_sub = false
end

if calc_2d
    ## calculate (complex) 3rd order corr function (with T=0)
    zp = 10 # zeropad up to 2^zp

    ## calculate 2D spectra at
    T = [0, 15] #fs

    spectra2d = Array{out2d}(undef, length(T))

    # cannot plot inside cmds.jl when multithreading
    #for i = 1:length(T)
    # multithreading (run several T steps in parallel)!
    Threads.@threads for i = 1:length(T)
        spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                            method;debug=true,use_sub=use_sub,
                                                zp=zp)
    end

    ## crop 2D data and increase dw
    spectra2d = [crop2d(spectra2d[i],4;w_max=6,step=1) for i = 1:length(T)]

    ## assign ω-axis from output
    #ω = out2d[1].ω
end

if calc_2d

    ## plot 2D spectra for each(?) T
    # what to plot
    rep="absorptive"
    scal="lin"

    ## make  subplot layout
    nplots = length(T);
    ncols = Int32(ceil(sqrt(nplots)));
    nrows = Int32(ceil(nplots / ncols));

    ## create figure
    #figure(figsize=(ncols*3+0.2,nrows*3+0.2));
    #suptitle(rep * " 2D spectrum (" * scal * ". scaling)")
    #for i = 1:nplots
    #    subplot(nrows,ncols,i,aspect="equal")
    #    cmds.plot2d(out2d[i].ω,round.(out2d[i].full2d,digits=1);repr=rep,scaling=scal)
    #    title("2D spectrum at $(T[i]) fs")
        #xlim([0, E₁+E₂]); ylim([0, E₁+E₂]);
    #end

    # determine maximum value in dataset out2d[:].full2d[:,:]
    maxi = maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)])

    # plot 2D spectra
    fig, ax = subplots(nrows,ncols,sharex=true,sharey=true,figsize=(ncols*3.2,nrows*3))
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

    ## plot TA (summed 2D spectrum)
    figure()
    ta = [sum(spectra2d[i].full2d,dims=1) for i in 1:length(spectra2d)]
    [plot(spectra2d[1].ω,ta[i]') for i in 1:length(ta)]

end



## How to check execution time (run twice after function assignment!):

"""
function test_timing()
    @time H = +([disorder[i] * (B_TLSs,i,H_TLS) for i in 1:num_of_TLSs]...);  # works with +, -, *, /
    #@time H = last(accumulate(+, [disorder[i] * embed(B_TLSs,i,H_TLS) for i in 1:num_of_TLSs]));  # works with +, -, *, /
    return H
end
test_timing()


# How to construct full Hamiltonian
# http://hitoshi.berkeley.edu/221A/tensorproduct.pdf

H = H_a ⊗ one(H_a) + one(H_a) ⊗ H_a + 0.1 * ((sm⊗sp)+(sp⊗sm))

#H = manybodyoperator(ManyBodyBasis(Σ,bosonstates(Σ,[1,1])),H_a)


H = H_a ⊗ one(H_a) ⊗ one(H_a) + one(H_a) ⊗ H_a ⊗ one(H_a) + one(H_a) ⊗ one(H_a) ⊗  H_a + 0.1 * ((sm⊗sp⊗one(sm))+(sp⊗sm⊗one(sm)) + (sm⊗one(sm)⊗sp)+(sp⊗one(sm)⊗sm) + (one(sm)⊗sm⊗sp)+(one(sm)⊗sp⊗sm))

Nm = 3


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
        disorder[j] = 1 #+ (rand(Float64, 1)[1] - 0.5) / 20
        global H_amb = H_amb + disorder[j] * permutesystems(H_am,circshift(R,j-1))
    end
end

for i in R
    if i == 1
        println(circshift(R,i-1))
        #global H_x = (one(spp) ⊗ permutesystems(smm,circshift(R,i-1))) +
        #            one(smm) ⊗ permutesystems(spp,circshift(R,i-1))
        #global H_x = (sp ⊗ sm) + (sm ⊗ sp)
                        #one(smm) ⊗ permutesystems(spp,circshift(R,i-1))
    else
        #global H_x = H_x + (one(spp) ⊗ permutesystems(smm,circshift(R,i-1))) +
        #            one(smm) ⊗ permutesystems(spp,circshift(R,i-1))
        #global H_x = H_x
    end
end

H = H_amb + 0.1*H_x



E = eigvals(dense(H).data)
xn = 0
figure(figsize=(5,4));
title("Energy levels (Eigenvalues) of H_TC")
#hlines([0, ωr], -.65, -.35)
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
"""



""" WAY TO CONSTRUCT THEM ... TOO DIFFICULT
num_of_exc = 0
elms = push!(repeat([nlevelstate(b_TLS,2)],num_of_exc),repeat([nlevelstate(b_TLS,1)],num_of_TLSs-num_of_exc)...)
ground = tensor(elms)

num_of_exc = 1
elms = push!(repeat([nlevelstate(b_TLS,2)],num_of_exc),repeat([nlevelstate(b_TLS,1)],num_of_TLSs-num_of_exc)...)
singl = push!([tensor(circshift(elms,i-1)...) for i in 1:num_of_TLSs])

num_of_exc = 2
elms = push!(repeat([nlevelstate(b_TLS,2)],num_of_exc),repeat([nlevelstate(b_TLS,1)],num_of_TLSs-num_of_exc)...)
#doubl = push!([tensor(circshift(elms,i-1)...) for i in 1:num_of_TLSs])
doubl = push!([tensor(unique(permutations(elms))[i]...) for i in 1:binomial(num_of_TLSs,2)])

test = Ket[]
test = push!(test,ground)
test = push!(test,singl...)
test = push!(test,doubl...)
"""
