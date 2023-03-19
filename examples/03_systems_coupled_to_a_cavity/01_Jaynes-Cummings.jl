using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles, Random

# make sure to set script directory as pwd()
cd(@__DIR__)

pygui(true)

calc_2d = true

cmp = create_colormap("bright");

Δ = -0.2
wc_k(θ) = (2.0 + Δ + sin(θ))  * 2 * pi # cavity frequency
wc = wc_k(0)
wa = 2.0  * 2 * pi  # atom frequency
g  = 0.15 * 2 * pi   # coupling strength
κ = 0.5         # cavity dissipation rate
γ = 0.00005         # atom dissipation rate
N = 1               # number of cavity fock states
n_th_a = 0.0        # avg number of thermal bath excitation ???
#use_rwa = True

tlist = [0:0.05:50;]
b_fock = FockBasis(N)
b_atom = NLevelBasis(2)
Psi0 = fockstate(b_fock,1) ⊗ nlevelstate(b_atom,1)    # start with an excited cavity

a = destroy(b_fock) ⊗ one(b_atom)
at = dagger(a)
sm = one(b_fock) ⊗ transition(b_atom,1,2)
sp = dagger(sm)

H = wc * at * a + wa * sp * sm + g * (at * sm + a * sp)

H_k(θ) = wc_k(θ) * at * a + wa * sp * sm + g * (at * sm + a * sp)

println(H)

#J = [sqrt(γ)*sm,sqrt(κ)*a]
# the same ?
J     = [a,sm]
rates = [sqrt(κ), sqrt(γ)]
#rates  = [1, 1]
tout, rhot = timeevolution.master(tlist,Psi0,H,rates .* J)
nc = real(expect(at * a, rhot))
na = real(expect(sp * sm, rhot))

figure(figsize=(6,6))
subplot(211)
plot(tout,nc)
plot(tout,na)

#ω, spec = timecorrelations.correlation2spectrum(tlist, nc; normalize_spec=true)
#figure(figsize=(6,3))
#plot(ω, spec)

rho0 = dm(Psi0)

corr = timecorrelations.correlation(tlist, rho0, H, J, at, a; rates=rates)
ω, spec = timecorrelations.correlation2spectrum(tlist, corr; normalize_spec=true)
plot(tlist,corr)

subplot(212)
plot(ω, spec)



width = .5
ωq_TDBC     = [H.data[3,3] - H.data[2,2]]
ampl_TDBC   = [1] / (N)
gam = .1

try
    local m = @which spectral_density(1)
    Base.delete_method(m)
catch
end
function spectral_density(ω) # power spectral density, thermal bath noise spectrum
    w = [(ω .- ωq_TDBC[i]) ./ (.5 * width) for i in 1:length(ωq_TDBC)]
    if ω == 0
        #Jω = 0
        Jω = gam
    elseif !isless(real(ω),0)
        Jω = sum([ampl_TDBC[i] ./ (1 .+ w[i].^2) for i in 1:length(ωq_TDBC)])  # probably most realistic ... ???
        #Jω = gam / 2 * (ω / (2 * pi)) * (ω > 0.0)
    elseif  isless(real(ω),0)
        Jω = 0
        #Jω = gam / 2 * (-ω / (2 * pi)) * (ω < 0.0)
    end
    return Jω
end

Γ_R = .01
a_ops = [Γ_R * sp * sm, spectral_density]# plot correlation function and spectrum
a_ops = [Γ_R * at * a,  spectral_density]# plot correlation function and spectrum

R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops]; J=rates .* J)


μ12 = at * sm*sp + a * sm*sp
μ23 = at * a * sp + at * a * sm

tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0,R,H)
corr = expect(μ12,rhot)

zp = 10

corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
ωnew, spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true);
plot(ωnew,spec)


rhot = nothing
tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0*μ12,R,H)

vecs   = eigvecs(dense(H).data)
eivecs = [Ket(H.basis_l,vecs[:,i]) for i in 1:length(vecs[1,:])]

#Operator(H.basis_l, H.basis_l, eivecs)

transf = copy(H)
transf.data = vecs

rhot_coupl = [transf' * i * transf for i in rhot]

H_coupl = transf' * H * transf

tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0*μ12,R,H)

vecs   = eigvecs(dense(H_coupl).data)
eivecs = [Ket(H_coupl.basis_l,vecs[:,i]) for i in 1:length(vecs[1,:])]

figure()
subplot(132);
title("eigen/coupled basis") # levels: 1 -> GS, 2 -> LP, 3 - N-1 -> ..., N -> UP
for i in 1:length(H.basis_l)
    es_n = real(expect(dm(eivecs[i]),rhot))
    plot(tout,es_n,linestyle="dashed", label=string(i)); #append!(leg,[string(i)])
    ylim((-0.1, 1))
end
legend()

vecs   = eigvecs(dense(H).data)
eivecs = [Ket(H.basis_l,vecs[:,i]) for i in 1:length(vecs[1,:])]

es_mat = zeros(length(tout))
subplot(133); title("uncoupled basis") # levels: 1 -> GS, 2 -> cav., 3 - N -> TLSs
for i in 1:length(H.basis_l)
    es_n = real(expect(dm(eivecs[i]),rhot))
    plot(tout,es_n,linestyle="dashed", label=string(i)); #append!(leg,[string(i)])
    if i > 2
        global es_mat = es_mat .+ es_n
    end
end
#plot(tout, es_mat, label="Σ(TLSs)")
legend()

corr = []
spec = []
angles = collect(-pi/4:0.001*pi:pi/2)
for i in angles
    #temp1 = timecorrelations.correlation(tlist, rho0, H_k(i), J, at, a; rates=rates)
    local R, ekets = timeevolution.bloch_redfield_tensor(H_k(i), [a_ops]; J=rates .* J)
    local tout, temp = timeevolution.master_bloch_redfield(tlist,μ12*rho0,R,H_k(i))
    temp1 = expect(μ12,temp)
    append!(corr, [temp1])
    local ω, temp2 = timecorrelations.correlation2spectrum(tlist, temp1; normalize_spec=false)
    append!(spec,[temp2])
end

spec = hcat(spec...)
figure()
contourf(angles ./ (2*pi) .* 360,ω, spec)

Psi0 = fockstate(b_fock,0) ⊗ nlevelstate(b_atom,1)    # start with system in gs
rho0 = dm(Psi0)

use_sub = false

if calc_2d
        ## calculate (complex) 3rd order corr function (with T=0)
        zp = 11 # zeropad up to 2^zp

        F = R

        ## calculate 2d spectra at
        T = [0., 10., 100.]

        spectra2d = Array{out2d}(undef, length(T))
        for i = 1:length(T)
            spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                                "redfield";debug=false,use_sub=use_sub,zp=zp);
        end

        ## crop 2D data and increase dw
            spectra2d = [crop2d(spectra2d[i],3.5;w_max=20,step=2) for i = 1:length(T)]

        ## assign ω-axis from output
        ω = spectra2d[1].ω

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
                plot2d(ω,spectra2d[i].full2d;repr=rep,scaling=scal)
                title("2D spectrum at $(T[i]) fs")
                #xlim([0, E₁+E₂]); ylim([0, E₁+E₂]);
        end
        tight_layout()

        #plot([E₁, E₁], [ω[1], ω[end]],"k--",linewidth=1,alpha=0.25)
        #plot([E₂, E₂], [ω[1], ω[end]],"k--",linewidth=1,alpha=0.25)
        #plot([ω[1], ω[end]],[E₁, E₁],"k--",linewidth=1,alpha=0.25)
        #plot([ω[1], ω[end]],[E₂, E₂],"k--",linewidth=1,alpha=0.25)
        #plot([energies[2:3], energies[2:3]], [ω[1], ω[end]],"g--",linewidth=1,alpha=0.25)
        #plot([ω[1], ω[end]],[energies[2:3], energies[2:3]],"g--",linewidth=1,alpha=0.25)


        # choose one of the following plotting ranges
        #xlim([ω[1], ω[end]]); ylim([ω[1], ω[end]]);
        #xlim([0, ω[end]]); ylim([0, ω[end]]);
        #xlim([0, wa+wc]); ylim([0, wa+wc]);
end
