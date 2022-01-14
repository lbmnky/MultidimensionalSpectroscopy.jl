#= 
This example follows the approach by Luer et al. doi:10.1021/acs.jpclett.6b02704 to simulate 
Frenkel excitons in TDBC J-aggregates.

Essentially it averages many different realisations of a J-aggregate, resulting in an
averaged, smooth lineshape. 

=#

using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles, Random, Combinatorics, Distributions

# make sure to set script directory as pwd()
cd(@__DIR__)

pygui(true)

calc_2d = false

try
    local m = @which dω(1)
    Base.delete_method(m)
catch
end

function dω(stdev)
    return rand_normal(0,stdev)
end

tlist = [0:5.:750;]
corr = zeros(length(tlist))

cmp = create_colormap("bright");

N = 70

#E = [2.41 * dE_n[i] for i in 1:N]
E_mon = 2.41
E_exc = 2.11
E_RWA = E_mon
E_RWA = E_exc
E     = E_mon - E_RWA

b_mon = NLevelBasis(2)
σ⁺ = transition(b_mon,2,1)
σ⁻ = transition(b_mon,1,2)

#rho0 = dm(tensor([nlevelstate(b_mon,1) for i in 1:N]))

B_exc = b_mon^N

H_n(n)  = embed(B_exc,n,E * σ⁺ * σ⁻) # Hamiltonian at site n => Hₙ 

Σ⁺_n(n) = embed(B_exc,n,σ⁺)              # raising operator at site n 
Σ⁻_n(n) = embed(B_exc,n,σ⁻)              # lowering operator at site n

#Σ⁺ = sum(Σ⁺_n.(1:N))
#Σ⁻ = sum(Σ⁻_n.(1:N))


N_deloc = 10
Δ = .299
#J = [-Δ / (2 * cos(pi /(N+1))) for i in 1:N]
J = -Δ / (2 * cos(pi /(N_deloc+1)))
J  = -0.15  # from Luer
J2 = -0.01 # from Luer
#J = 0
#J2 = 0
er = 3
e0 = 8.854e-21 #F/nm
ke = (4 * pi * er * e0)^-1
μd = 1.869e-20 # Cm

R = 1

#J = ke * (μd^2 / R^3) 

#if N > 3
#    H = sum(H_n.(1:N)) + 
#                sum(J[i] * (Σ⁺_n(i) * Σ⁻_n(i+1) + Σ⁻_n(i) * Σ⁺_n(i+1)) for i in 1:N-1) 
#end

b = NLevelBasis(N)
p = transition(b,2,1)
H = dagger(p) * p * 0

rho0 = copy(H)
rho0.data[1,1] = 1

μ    = copy(H)

Ltemp   = copy(H)
L = []

for i in 2:N
    H.data[i,i] = E
    if i > 2
        H.data[i,i-1] = J
        H.data[i-1,i] = J
    end
    if i > 3
        H.data[i,i-2] = J2 # J2 coupling
        H.data[i-2,i] = J2
    end
    μ.data[1,i] = 1
    μ.data[i,1] = 1

    Ltemp.data[1,i] = 1
    append!(L,[copy(Ltemp)])
    Ltemp.data[1,i] = 0

end


#μ    = sum((Σ⁻_n(i) + Σ⁺_n(i)) for i in 1:N)

μ12  = rho0 * μ + μ * rho0
μ23  = copy(μ12)

rho1 = μ12 * rho0 * μ12
μ12  = μ12 / sqrt(tr(rho1))
rho1 = μ12 * rho0 * μ12
rho2 = μ23 * rho1 * μ23
#μ23  = μ23 / sqrt(tr(rho2))
rho2 = μ23 * rho1 * μ23

# TDBC details https://arxiv.org/pdf/1506.01321.pdf
γ01 = 0.00476 # eV ... from 1.15 # ps⁻¹ from photoluminescence
Γ01d = 0.017  # eV from 17 meV 
Γ01 = γ01 / 2 + Γ01d
#L = Γ01 .* [Σ⁻_n(i) * (Σ⁺_n(i) * Σ⁻_n(i)) for i in 1:4]
L = Γ01 .* L


#=
idx = sortperm(real(diag((H).data)))
H.data = H.data[idx,idx]
Σ⁺.data = Σ⁺.data[idx,idx]
Σ⁻.data = Σ⁻.data[idx,idx]
rho0.data = rho0.data[idx,idx]
rho1.data = rho1.data[idx,idx]
rho2.data = rho2.data[idx,idx]
for i in 1:N
    #Σ⁺n[i].data = Σ⁺n[i].data[idx,idx]
    #Σ⁻n[i].data = Σ⁻n[i].data[idx,idx]
    L[i].data = L[i].data[idx,idx]
end

H_si, transf_mat_si, P_si, L_si, rho0_si, rho1_si, rho2_si, μ12_si, μ23_si, Σ⁺_si, Σ⁻_si = 
    cmds.create_subspace([H],"si", L, rho0, rho1, rho2, μ12, μ23, Σ⁺, Σ⁻)

##

H_si = transf_mat_si * H_si[1] * transf_mat_si'
rho0_si = transf_mat_si * rho0_si * transf_mat_si'
rho1_si = transf_mat_si * rho1_si * transf_mat_si'
rho2_si = transf_mat_si * rho2_si * transf_mat_si'
μ12_si = transf_mat_si * μ12_si * transf_mat_si'
μ23_si = transf_mat_si * μ23_si * transf_mat_si'
Σ⁺_si = transf_mat_si * Σ⁺_si * transf_mat_si'
Σ⁻_si = transf_mat_si * Σ⁻_si * transf_mat_si'
=#


##
figure()

D = 1/100
tc = 100
gt = D^2 * tc^2 * (exp.(-tlist./tc) .+ tlist./tc .- 1)
#gt = tc^2 * (exp.(-tlist./tc) .+ tlist./tc .- 1)

F = L

M = 500
for i in 1:M
    println(i)

    global H_in_loop = copy(H)

    # intrachain disorder
    local dE_n = [dω(0.022) for i in 1:N]

    # interchain disorder
    local dE = dω(0.022)
    println(E + dE_n[1])

    #plot(collect(1:N),E .+ dE_n)

    for j in 2:N
        H_in_loop.data[j,j] = (H_in_loop.data[j,j] + dE_n[j-1] ) #+ dE
    end

    # introduce defect sites
    global N_defect = 0
    if rand() > 0.
        N_defect = 3
    else
        N_defect = 0
    end
    #global N_defect = rand(1:10)    
    global N_defect = Int8(round(rand(Normal(5, .01))))
    if N_defect <= 0
        N_defect = 1
    end

    global ind_defect = rand(2:(N-1),N_defect)
    println(ind_defect)
    #if maximum(ind_defect) != N
        for j in ind_defect
            H_in_loop.data[j,j+1] = .6 * H_in_loop.data[j,j+1]
            H_in_loop.data[j+1,j] = .6 * H_in_loop.data[j+1,j]
            H_in_loop.data[j,j-1] = .6 * H_in_loop.data[j,j-1]
            H_in_loop.data[j-1,j] = .6 * H_in_loop.data[j-1,j]
            if j < N-1
                H.data[j,j+2] = .3 * H.data[j,j+2]
                H.data[j+2,j] = .3 * H.data[j+2,j]
            end
        end
    #else
    #end

    println("start corr")


    global corr = ( (timecorrelations.correlation(tlist, dense(rho0), H_in_loop, L, μ12, μ12)) + corr ) ./ 2

    #plot(tlist,corr)
    if i != M
        if i == 1
            plot([1, N], -1 .* [J, J],label="J")
        end
        scatter(collect(1:N) .+ rand(),dE_n .+ E_RWA,s=15,alpha=1/M*10,color="gray")
        scatter(ind_defect .+ .5,[.6 * -J for i in 1:length(ind_defect)], s=15,alpha=.01,color="red")
    elseif i == M
        scatter(collect(1:N) .+ rand(),dE_n .+ E_RWA,s=15,alpha=.25,color="black",label="site energies")
        scatter(ind_defect .+ .5,[.6 * -J for i in 1:length(ind_defect)], s=15,alpha=.25,color="red",label="defects")
    end



    if calc_2d
        if i < 3
            local T = [0]
            global spectra2d = Array{out2d}(undef, length(T))
            spectra2d[1] = make2Dspectra(tlist,[rho0, rho0],[H_in_loop, H_in_loop],[F, F],[μ12, μ12],[μ12, μ12],T[1],
                                                "lindblad";debug=false,use_sub=true,
                                                    t2coh=false,zp=10)
            spectra2d[1].ω = spectra2d[1].ω .+ E_RWA

            if i == 1
                global corr2d = spectra2d[1].corr
            else
                corr2d = (spectra2d[1].corr + corr2d) ./ 2
            end
        end
    end
    
end
ylim((0,3))
legend()

#corr /= M
corr =corr .* exp.(-gt)

freqRWA = exp.(-1im * E_RWA .* tlist) .* exp.(1im * E_RWA .* vcat(-reverse(tlist)[1:end-1],tlist))'

if calc_2d
    #corr2d = corr2d .* exp.(-gt) .* vcat(reverse(exp.(-gt))[1:end-1],exp.(-gt))'
    #spec2d  = corr2spec(corr2d .* freqRWA,10)

    figure()
    pcolor(real(spec2d))
end

zp = 11
corr    = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

ω_exc, spec_exc = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true);

##

figure();
subplot(211)
plot(tnew,corr)
plot(tlist,exp.(-gt))

subplot(212)
plot(-1 .* (ω_exc.-E_RWA),spec_exc)

dataTDBC_mon = readdlm("TDBCmon_abs.dat", ',')
dataTDBC_aggr = readdlm("TDBCaggr_abs.dat", ',')

plot(dataTDBC_mon[:,1],dataTDBC_mon[:,2]/maximum(dataTDBC_mon[:,2]),"k:")
plot(dataTDBC_aggr[:,1],dataTDBC_aggr[:,2]/maximum(dataTDBC_aggr[:,2]),"r:")

