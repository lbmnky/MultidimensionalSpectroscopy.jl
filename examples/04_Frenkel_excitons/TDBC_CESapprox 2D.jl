#=
Coherent exciton scattering (CES) approximation (Roden, Eisfeld, Briggs)

* mean-field type approximation
* xact monomer Green function (electronic and vibrational DOF) is replaced by its average in the vibrational ground state
* reproduces absorption and circ. dichroism of H- and J-aggregates

=#

using PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles, MultidimensionalSpectroscopy

# make sure to set script directory as pwd()
cd(@__DIR__)

# to use PyPlot in GUI mode under Linux
pygui(true)
# run once with calc_2d = false to initialize functions
calc_2d = false   # !! use less values for k if running 2D 

Ω(k)      = Emon + 2 * V * cos( pi * k / (N + 1) )
Ω2(k1,k2) = Emon + 2 * V * cos( pi * k1 / (N + 1) ) + Emon + 2 * V * cos( pi * k2 / (N + 1) )
δ(l,m) = l == m
# sum over 2 indices that are not equal
[i+j for i in 1:2 for j in 1:2 if i != j]
# https://aip.scitation.org/doi/pdf/10.1063/1.469393
μμ(k)        = μ_mon * sqrt( 2 / (N+1)) *             (1 - (-1)^k ) / 2 * cot( pi * k  / (2 * (N+1)))
μμ2(k1,k2,k) = μ_mon * sqrt( 2 / (N+1)) * ( δ(k2,k) * (1 + (-1)^k1) / 2 * cot( pi * k1 / (2 * (N+1))) - 
                                            δ(k1,k) * (1 + (-1)^k2) / 2 * cot( pi * k2 / (2 * (N+1))) +
                                            1/2 * ( δ(k1-k2-k,0) - δ(k1-k2+k,0) ) * 
                                                  ( cot(pi * k1 / (2 * (N+1))) + cot(pi * k2 / (2 * (N+1))) ) +
                                            1/2 * ( δ(k1+k2+k,2*(N+1)) - δ(k1+k2-k,0) ) *
                                                  ( cot(pi * k1 / (2 * (N+1))) - cot(pi * k2 / (2 * (N+1))) ) )

Sp(i) = transition(B_TDBC,1+i,i) 
Sm(i) = transition(B_TDBC,i,1+i)

k  = [1, 2, 3, 5, 7, 9] 
k  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] 
k1 = [1, 2]#, 3, 4, 5, 6, 7]
k2 = [1, 2]#, 3, 4, 5, 6, 7]
K  = [(2,1),(3,2), (4,1), (6,1), (5,2)]

V = -2455 / 2 / 8065
μ_mon = 1
Emon = 2.41532 # 513 nm ... https://aip.scitation.org/doi/pdf/10.1063/1.469393
Emon = 2.41032 # 513 nm ... https://aip.scitation.org/doi/pdf/10.1063/1.469393

B_TDBC = NLevelBasis(1 + length(k) + length(K))

try     # need to delete method first when working on it
    m = @which make_Ops(1)
    Base.delete_method(m)   
catch
end
function make_Ops(num)
    global N = num
    ##
    ωs_ge  = [Ω(i) for i in k] 
    #ωs_ef  = [Ω2(i,j) for i in k for j in k2 if i > j]
    ωs_ef  = [Ω2(i...) for i in K]
    ωs     = vcat(ωs_ge,ωs_ef)

    μs_ge = [μμ(i) for i in k]
    #μs_ge = sqrt.(complex(μs_ge))
    μs_ef = zeros(length(k),length(K))
    for l in 1:length(K)
        for i in 1:length(k)
            μs_ef[i,l] = μμ2(K[l]...,k[i]) 
        end
    end

    
    H   = sum(ωs[i] * transition(B_TDBC,1+i,1) * transition(B_TDBC,1,1+i) for i in 1:length(ωs))
    μ12 = sum(μs_ge[i] * transition(B_TDBC,1+i,1) for i in 1:length(k))
    μ23 = sum(μs_ef[i,j] * transition(B_TDBC,1+length(k)+j,1+i) for j in 1:length(K) for i in 1:length(k))

    
    
    ## 

    Σ⁺ge = μ12
    Σ⁺ef = μ23
    Σ⁺ = μ12 + μ23
    #Σ⁺ = sum(transition(testb,1+i,1) for i in 1:length(k))
    

    #Σ⁺ = sum(Sp(i) for i in k)
    Σ⁻ = dagger(Σ⁺);
    μ12 = (μ12) + dagger(μ12);
    #μ12.data = sqrt.(μ12.data)
    μ23 = μ23 + dagger(μ23)
    #μ23.data = sqrt.(μ23.data)

    # make initial density matrix
    Psi0 = nlevelstate(B_TDBC,1)
    rho0 = dm(Psi0)

    # make collapse operator
    #ccc = diagonaloperator(b_TLS,[1, 1/4, 1])

    #     relaxation        dephasing
    #    _____________    _______________   
    #L = [j12, j13, j23,  j21*j12, j31*j13]
    #Γ = [.01,.00 ,.05,    .18,     .01]
    #L = Γ .* L
    #...redo for new system
    #L = [.05 * cmds.tri(μ12,"L") * cmds.tri(μ12,"U")]
    L1 = 0.01 .* [transition(B_TDBC,1,1+i) for i in 1:length(k)]
    L2 = 0.09 .* [transition(B_TDBC,2,2+i) for i in 1:length(k)]

    L_dephasing     = [.04 * Sp(i) * Sm(i) for i in 1:length(k)]
    L_dephasingESA  = [.06 * Sp(i) * Sm(i) for i in length(k):length(k)+length(K)]

    L = vcat(L1,L2,L_dephasing)
    L = vcat(L1,L2)

    ##

    rho1 = μ12 * rho0 * μ12
    μ12  = μ12 / sqrt(tr(rho1))
    μ23  = μ23 / sqrt(tr(rho1))
    rho1 = μ12 * rho0 * μ12
    rho2 = μ23 * rho1 * μ23
    #μ23  = μ23 / sqrt(tr(rho2)) 
    rho2 = μ23 * rho1 * μ23

    return H, μ12, μ23, Σ⁺ge, Σ⁺ef, Σ⁺, Σ⁻, μs_ge, μs_ef, ωs_ge, ωs_ef, rho0, rho1, rho2, L

end

H, μ12, μ23, Σ⁺ge, Σ⁺ef, Σ⁺, Σ⁻, μs_ge, μs_ef, ωs_ge, ωs_ef, rho0, rho1, rho2, L = make_Ops(14)



##
#figure()
#bar(ωs_ge,μs_ge.^2,width=0.003)
##



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

#V = -2505 / 2
#Eg_k1 = Emon + 2*V * cos(pi/(20+1)) / 8065
#E = Eg_k1

#Ek1_k1k2 = Emon + 2*V * cos(2*pi/(20+1)) / 8065

#E2 = Eg_k1 + Ek1_k1k2

# get eigenstates and eigenvalues of H
#  using Julia LinearAlgebra functions instead of QuantumOptics requires to use H.data
#states   = eigvecs(dense(H).data)
#energies = eigvals(dense(H).data)

# set up figure
fig  = figure(figsize=(15,9));
ax1  = subplot(241); title("Energy levels")
ax1a = subplot(242); title("Noise-power spectrum S(ω)")
ax2  = subplot(243); title("correlation function LB")
ax2a = subplot(244); title("correlation function RF")
ax3  = subplot(223); title("Kinetics in eigen/coupled basis") 
ax4  = subplot(247); title("Abs. spectrum LB")
ax4a = subplot(248); title("Abs. spectrum RF")



sca(ax1)
plot_levels(H,0)
hlines(.5*E2,-.45,.45, color="red",linestyle="dashed")
hlines(Emon,-.45,.45, color="green",linestyle="dashed")
hlines(2*Emon,-.45,.45, color="green",linestyle="dashed")


# time list for evaluation of master equation
Tmax = 400
dt = 1.
tlist = [0:dt:Tmax;]

#M = 1
#corr = sum(timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12) for i in 1:M) ./ M
##
D = .1
tc = .9
#gt = D^2 * tc^2 * (exp.(-collect(-150:0.2:150)./tc) .+ collect(-150:0.2:150)./tc .- 1)
gt = D^2 * tc^2 * (exp.(-tlist./tc) .+ tlist./tc .- 1)

Δ = .01 # root-mean-square amplitude of ⟨δω(t)|δω(0)⟩
g_slow = Δ^2 .* tlist.^2 ./2 # static limit
Λ = 100 # inverse correlation time (tc)⁻¹
g_fast = Δ^2 .* tlist ./ Λ   # Markovian limit

# or ...
#Δ_slow = 10 * 0.00414
#Λ_slow = 1 
#Δ_fast = 54 * 0.00414
#Λ_fast = 5
#
#g_slow = Δ_slow^2 / Λ_slow^2 .* (exp.(-Λ_slow .* tlist) .+ Λ_slow .* tlist .-1)
#g_fast = Δ_fast^2 / Λ_fast^2 .* (exp.(-Λ_fast .* tlist) .+ Λ_fast .* tlist .-1)

# from MUKAMEL Chapter 8

#gt = 2 * λ * kb * T / (ħ * Λ^2) * (exp.(-Λ .* tlist) .+ Λ .* tlist .- 1) .
#        - 1im * (λ/Λ) .* (exp.(-Λ .* tlist) .+ Λ .* tlist .- 1)

κ = 3
Λ = .02

λ = 0.2 # half the Stokes shift in eV !?
kB = 8.617333262145 * 10^-5 # in eV * K^-1
T = 300 # Kelvin
ħ = 6.582119569 * 10^-16 # in eV s
#Δ = sqrt(2 * λ * kB * T / ħ)


Δ = Λ / κ # unit 1/s

gt = Δ^2 / Λ^2 * (exp.(-Λ .* tlist) .+ Λ .* tlist .- 1) 

#gt = g_slow + g_fast


corr = zeros(size(tlist))
list = [12, 13, 14, 15, 16, 17, 18]
list = collect(14:.5:16)
for s in list
    local H, μ12, μ23, Σ⁺ge, Σ⁺ef, Σ⁺, Σ⁻, μs_ge, μs_ef, ωs_ge, ωs_ef, rho0, rho1, rho2, L = make_Ops(s);
    global corr = exp(-(s-15)^2) .* (timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12) + corr) .* exp.(-gt) ./ 2
end

corr = corr ./ maximum(real(corr))


sca(ax2)
##
corr = corr .* exp.(-gt)

tout, rhot = timeevolution.master(tlist,μ12 * rho0,H,L)

# zeropad corr and extrapolate tlist
zp = 11
corr = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
ω,spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)

#figure(figsize=(6,3));
# add the results to previous figure
sca(ax2)
plot(tnew,real(corr))
plot(tlist,exp.(-gt))

xlabel("time"); ylabel("Corr. func. ⟨μ(t) μ₀₁ ρ₀⟩")

sca(ax4)
plot(-ω,spec)
xlabel("ω"); ylabel("Absorption")

tout, rhot = timeevolution.master(tlist,μ23 * rho1,H,L)
corr = sum(timecorrelations.correlation(tlist, rho1, H, L, μ23, μ23) for i in 1:N) ./ N
corr = corr .* exp.(-gt)

corr = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
ω,specESA = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)
#plot(-[E, E], [0, .5])
plot(-ω,-specESA)
xlabel("ω"); ylabel("Absorption")




tout, rhot = timeevolution.master(tlist,μ12 * rho0 * μ12,H,L)

vecs   = eigvecs(dense(H).data)
eivecs = [Ket(H.basis_l,vecs[:,i]) for i in 1:length(vecs[1,:])]

sca(ax3)
# levels: 1 -> GS, 2 -> LP, 3 - N-1 -> ..., N -> UP
for i in 1:length(eivecs)
    es_n = real(expect(dm(eivecs[i]),rhot))
    plot(tout,es_n,linestyle="dashed", label=string(i)); #append!(leg,[string(i)])
    ylim((-0.1, 1))
end
#legend()



#

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
        Jω = .02
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
        Jω = 0.00000
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


Γ = .2
#a_ops = [Γ * (j21+j31+j32) * (j12+j13+j23), spectral_density]
dat = tri(tri(Σ⁺ * Σ⁻,"U"),"L")
a_ops = [Γ * Σ⁺ * Σ⁻, spectral_density]
a_ops = [Γ * Σ⁺ge * dagger(Σ⁺), spectral_density]#, Γ * Σ⁺ef * dagger(Σ⁺ef), spectral_density]
#a_ops = [dat, spectral_density]

R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops]; J=L)

# linear absorption
tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0,R,H)
corr = expect(μ12,rhot)
#corr = corr .* exp.(-gt)
#plot(tlist,corr .* exp.(-gt))

# calculate and plot spectrum
corr = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
sca(ax2a)
plot(tnew,real(corr))
ω, spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true);
sca(ax4a)
plot(-ω,spec)

# ESA
tout, rhot = timeevolution.master_bloch_redfield(tlist,μ23*rho1,R,H)
corr = expect(μ23,rhot)
#corr = corr .* exp.(-gt)
#plot(tlist,corr .* exp.(-gt))

# calculate and plot spectrum
corr = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)
sca(ax2a)
plot(tnew,real(corr))
ω, specESA = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true);
sca(ax4a)
plot(-ω,-specESA)



#figure()
sca(ax1a)
plot(collect(-1:.001:1),(spectral_density.(collect(-1:.001:1))))

# load TA
TDBC_TA_t0 = readdlm("TDBC_TA_t0.dat",',') # https://aip.scitation.org/doi/pdf/10.1063/1.469393
TDBC_agg_abs = readdlm("TDBCaggr_abs2.dat",',') 
calc_ω = ω[ω.<0]
calc_abs = spec[ω.<0]

#convert to energy axis if it is in wavelength
if TDBC_agg_abs[1,1] > 100
    TDBC_agg_abs[:,1] = 6.582119569 * 10^-16 * 3e8 ./ (TDBC_agg_abs[:,1] * 1e-9) * 2 * pi
end

sca(ax4a)
plot(TDBC_agg_abs[:,1],TDBC_agg_abs[:,2],"k--")
bar(ωs_ge,(μs_ge./maximum(μs_ge)).^2,width=0.01,color="k")
#bar(ωs_ef,-(μs_ef[1]./maximum(μs_ef[1])).^2,width=0.01,color="k")

# scale ESA tdm s by ge tdm
μs_ef = μs_ef
temp = [μs_ef[:,i] ./ maximum(μs_ef[:,i]) .* (μs_ge./maximum(μs_ge))[i] for i in 1:length(ωs_ef)]
ωtemp = [[i .- j for j in ωs_ge] for i in ωs_ef]

for i in 1:length(temp)
    #bar(ωs_ef[i] .- ωs_ef[1] .+ ωs_ge,-temp[i], width=0.01)
    bar(ωtemp[i],-temp[i], width=0.01)
end

sca(ax4)
plot(TDBC_agg_abs[:,1],TDBC_agg_abs[:,2],"k--")
bar(ωs_ge,(μs_ge./maximum(μs_ge)).^2,width=0.01,color="k")
for i in 1:length(temp)
    bar(ωtemp[i],-temp[i], width=0.01)
end

tight_layout()

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
        T = [0, 50, 100, 200, 500, 1000] #0.66 fs

        spectra2d = Array{out2d}(undef, length(T))

        N = 50
        Threads.@threads for i = 1:length(T)
        #for i = 1:length(T)
            spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                                method;debug=false,use_sub=use_sub,t2coh="kin",zp=zp);
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

        # determine maximum value in dataset out2d[:].full2d[:,:]
        maxi = maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)])

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
        #cmds.save_2d([round.(real(out2d[i].full2d),digits=1) for i = 1:length(T)],T,fn)

        ## plot TA (summed 2D spectrum)
        figure()
        ta = [sum(spectra2d[i].full2d,dims=1) for i in 1:length(spectra2d)] 
        ta = vcat(ta...)'
        ta = ta ./ maximum(real(ta))
        plot(ω,ta)
        plot(ω,zeros(size(ω)),"k",linewidth = .5)
        plot(1239.8 ./ TDBC_TA_t0[:,1], -1 * TDBC_TA_t0[:,2]/maximum(-TDBC_TA_t0[:,2]) ,"k--")
        plot(TDBC_agg_abs[:,1],TDBC_agg_abs[:,2] ./ maximum(TDBC_agg_abs[:,2]),"r--")
        plot(-calc_ω,calc_abs./maximum(calc_abs[-calc_ω.>0]),"red")
        bar(ωs_ge,(μs_ge./maximum(μs_ge)).^2,width=0.005)
        #[bar(ωs_ef[i] .- ωs_ge,-μs_ef[:,i].^2 .* (μs_ge./maximum(μs_ge)).^2, width=0.005) for i in 1:length(ωs_ef)]
        for i in 1:length(temp)
            bar(ωtemp[i],-temp[i], width=0.01)
        end
        xlabel("Energy/frequency")
        ylabel("Difference absorption")
        tight_layout()

end
