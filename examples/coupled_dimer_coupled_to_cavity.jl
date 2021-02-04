#!/usr/bin/julia
using PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles

# make sure to set script directory as pwd()
cd(@__DIR__)

# include my custom cmds module
if Sys.iswindows()
    include("..\\cmds.jl")
    fn = "01_Output\\"
else
    include("../cmds.jl")
    fn = "01_Output/"
end

# to use PyPlot in GUI mode under Linux
pygui(true)
# run once with calc_2d = false to initialize functions
calc_2d = true

cmp = cmds.create_colormap("bright");

function draw_dipole(x,y,α,d)
    quiver(x- d/2*cos(α), y - d/2*sin(α), d*cos(α), d*sin(α), angles="xy",
            scale_units="xy",scale=1,color="r")
    scatter(x,y)
end

function angle2rad(α)
    A = α / 360 * 2 * π
    return A
end

## geometry of coupled dimer
# here, I use (0,0) as origin for point dipole 1 and an angle of the dipole of 0, such
# that it lies along the x-axis and has a length d. The point dipole is located at
# a distance r under the angle β, and has an dipole angle α and a length d. d
# controls the coupling strength.
#
# position of dipole 1 and angle wrt to x-axis (random)
x₁, y₁, α₁, d₁ = 0, 0, angle2rad(0), 0.25
# r: distance, β: slip angle
r , β          = 1.155, angle2rad(23.147495) # slip angle for TDBC 0.404 rad https://arxiv.org/pdf/1202.5712.pdf (theory)
# position of dipole 2 (calculated)
x₂, y₂         = x₁ + r*cos(β), y₁ + r*sin(β)
# angle and moment of dipole 2
α₂, d₂         = angle2rad(0), 0.25

## set up energies and couplings
# homo dimer    : E1 == E2 and J1 = J2
# hetero dimer  : E1 != E2
E₁ = 2.396
E₂ = 2.396

## draw dimer
#
figure(figsize=(10,6));
ax1 = subplot(331);
grid(b=1,which="both")
plot([x₁,x₂], [y₁,y₂], linestyle="--", color="k")
draw_dipole(x₁,y₁,α₁,d₁); draw_dipole(x₂,y₂,α₂,d₂)
xlim([-0.5, x₂+.5]), ylim([-0.5, y₂+0.5])
title("Dimer geometry")
xlabel("X"); ylabel("Y")

## determine x and y-components of tdm 1 and 2
d₁xy = [cos(α₁)*d₁, sin(α₁)*d₁] # x and y component of point dipole 1
d₂xy = [cos(α₂)*d₂, sin(α₂)*d₂] # x and y component of point dipole 2

## coupling strength will depend on geometry, first cos is out-of-plane rotation
# WHAT IS THE FACTOR (here 0.4) in κ ???
κ = 1 * (cos(0) - 3 * cos(β) * cos(β))
J = κ * (d₁*d₂) / norm([x₁,y₁]-[x₂,y₂])^3

# the following from Cho's book (p.170)
R12 = [x₁,y₁]-[x₂,y₂]
J = dot(d₁xy,d₂xy) / norm(R12)^3 -
        3 * dot(d₁xy,R12) * dot(R12,d₂xy) / norm(R12)^5

## make Hilbert-space
b_mon = NLevelBasis(2)  # Hilbert-space of monomer
b_dim = b_mon ⊗ b_mon   # Hilbert-space of dimer

## create transition operators in monomer basis
j12 = transition(b_mon,1,2) # == g_e, or a / σ₋ (annihilation)
j21 = transition(b_mon,2,1) # == e_g, or a† /σ⁺ (creation)
# g_e / e_g notation from Constantin script

## convert transition operators into dimer basis
mu1hat = j21 ⊗ one(b_mon) + j12 ⊗ one(b_mon)
mu2hat = one(b_mon) ⊗ j21 + one(b_mon) ⊗ j12

## add geometry and magnitude of tdm and calculate total tdm in x and y
μx = d₁xy[1] * mu1hat + d₂xy[1] * mu2hat
μy = d₁xy[2] * mu1hat + d₂xy[2] * mu2hat

## visualize combined tdm (offset in y a bit for clarity)
quiver((x₁+x₂)/2-(d₁xy[1]+d₂xy[1])/2, (y₁+y₂)/2-(d₁xy[2]+d₂xy[2])/2+0.05,
        d₁xy[1]+d₂xy[1], d₁xy[2]+d₂xy[2], angles="xy", scale_units="xy",
        scale=1, color="b", alpha=0.75)

quiver((x₁+x₂)/2-(d₁xy[1]-d₂xy[1])/2, (y₁+y₂)/2-(d₁xy[2]-d₂xy[2])/2+0.05,
        d₁xy[1]-d₂xy[1], d₁xy[2]-d₂xy[2], angles="xy", scale_units="xy",
        scale=1, color="r")

scatter((x₁+x₂)/2,(y₁+y₂)/2)

## define Hamiltonian of the dimer
# H₁ and H₂ are the Hamiltonians of the individual dipoles
H₁ = E₁ * j21 * j12
H₂ = E₂ * j21 * j12
# Hₓ describes the exchange between dipoles 1 and 2
Hₓ = J * (j12 ⊗ j21 + j21 ⊗ j12)
# use the following with TDMs between ground- and biexciton state, this modifies
# the energies of the ground- and doubly excited state slightly
#Hₓ = J * ((j21+j12) ⊗ (j21+j12)) #???????????
# H is the total Hamiltonian of the system
H_dim = one(b_mon) ⊗ H₁ + H₂ ⊗ one(b_mon) + Hₓ


## ADD CAVITY

wcav = 2.2 # cavity frequency
g  = 0.15   # coupling strength
κ = 0.5         # cavity dissipation rate
N = 2               # number of cavity fock states
b_fock = FockBasis(N)
a = destroy(b_fock) ⊗ one(b_dim)
at = sparse(dagger(a))

J12a = j12 ⊗ one(j12)
J12b = one(j12) ⊗ j12
J21a = j21 ⊗ one(j21)
J21b = one(j21) ⊗ j21

sm = one(b_fock) ⊗ (j12 ⊗ j21)
sp = sparse(dagger(sm))

A = at * (one(b_fock) ⊗ J12a) + at * (one(b_fock) ⊗ J12b)
B = a  * (one(b_fock) ⊗ J21a) + a  * (one(b_fock) ⊗ J21b)

H_cav = wcav * at * a

H = H_cav + one(b_fock) ⊗ H_dim + g * (A + B)


## get eigenstates and eigenvalues of H
#  using Julia LinearAlgebra functions instead of QuantumOptics requires to use H.data
states_dim   = eigvecs(dense(H_dim).data)
energies_dim = eigvals(dense(H_dim).data)

states_cav   = eigvecs(dense(H_cav).data)
energies_cav = eigvals(dense(H_cav).data)

states   = eigvecs(dense(H).data)
energies = eigvals(dense(H).data)

#states = eigenstates(dense(H))

## perform bases transformation of μ
#  for this I need to express the matrices as dense, data matrices
#μx =  states * μx * states
#μy =  states' * dense(μy) * states

## plot site basis and exciton basis energy diagrams
ax2 = subplot(132)
hlines([0, E₁], -1.25, -0.75)
hlines([0, E₂],-0.25, 0.25)
hlines(energies_dim[:], 0.75, 1.25)
hlines(energies[:], 1.75, 2.25)
hlines(energies_cav[:], 2.75, 3.25)
xticks([-1, 0, 1, 2, 3], ["monomer 1", "monomer 2", "dimer", "dimer in cavity", "cavity"],
        rotation=45, ha="right")
title("Energy levels")

## make initial density matrix
Psi0 = fockstate(b_fock,0) ⊗ nlevelstate(b_mon,1) ⊗ nlevelstate(b_mon,1)
rho0 = dm(Psi0)
println("rho0")
display(rho0)

## time list for evaluation of master equation
tlist = [0:0.2:150;]

## make collapse operator
ccc = diagonaloperator(b_mon,[1, 1])
#ccc = one(b_mon)
#c_ops = [sqrt(0.25)*one(b_mon)⊗j12, sqrt(0.25)*j12⊗one(b_mon)]
L = [one(b_fock)⊗ccc⊗j12, one(b_fock)⊗j12⊗ccc, a]

Γ_L = [sqrt(0.005), sqrt(0.005),               sqrt(0.02)]
Γ_a = [sqrt(0.2), sqrt(0.2), sqrt(0.2), sqrt(0.2)]

L = Γ_L .* L
## make transition dipole operator FOR DIMER SYSTEM !!!!!
# clearly now μ has a direction dependence, but if sample is isotropic  take:
μ   = (μx + μy)
# and the tdm with the ground state:
#μ12 = one(b_fock) ⊗ (μ .* (j12*j21⊗j21 + j12*j21⊗j12 + j21⊗(j12*j21) + j12⊗(j12*j21)))
#display(dense(μ12))
# ... and with the doubly excited state:
#μ23 = one(b_fock) ⊗ (μ .* (j21*j12⊗j21 + j21*j12⊗j12 + j21⊗(j21*j12) + j12⊗(j21*j12)))
#display(dense(μ23)) 
CC = ccc ⊗ ccc
## make transition dipole operator FOR CAVITY !!!!!
μ12 = at * (one(b_fock) ⊗ (J12a * J21a * CC)) + a * (one(b_fock) ⊗ (J12a * J21a * CC))
μ23 = at * (one(b_fock) ⊗ (J21a * J12a * CC)) + a * (one(b_fock) ⊗ (J21a * J12a * CC)) +
      at * (one(b_fock) ⊗ (J21b * J12b * CC)) + a * (one(b_fock) ⊗ (J21b * J12b* CC))
μ = μ12 + μ23

rho1 = μ12 * rho0 * μ12
rho1 = μ12 * rho0 * μ12

rho2 = μ23 * rho1 * μ23
μ23  = μ23 / sqrt(tr(rho1))
rho2 = μ23 * rho1 * μ23
μ12  = μ12 / sqrt(tr(rho1))

# calculate magnitude of tdm with singly excited states e1 and e2
# HOT TO ORDER THEM ???? ... depending on J>0 (H) or J<0 (J) ? Something more
# fundamental possible ?
if J < 0
        # tdm magnitude <0,0|μ|0,1>
        Ja = states[1,:]' * dense(μ).data * states[2,:]
        # tdm magnitude <0,0|μ|1,0>
        Jb = states[1,:]' * dense(μ).data * states[3,:]
elseif J > 0
        # everything reversed
        Ja = states[1,:]' * dense(μ).data * states[3,:]
        Jb = states[1,:]' * dense(μ).data * states[2,:]
end

# plot
ax3 = subplot(334)
bar(0, real(Ja))
bar(1, real(Jb))
xticks([0, 1], ["μ_g-e1", "μ_g-e2"])
xlim([-1, 2]); ylim([0, 1.5 * maximum(real([Ja, Jb]))])
ylabel("tdm strength")

## calculate and plot expectation values of population
#=
tout, rhot = timeevolution.master(tlist, rho0, H, F)
nc = real(expect(j12 ⊗ j21, rhot))
na = real(expect(j21 ⊗ j12, rhot))
figure(figsize=(6,3));
plot(tout,nc)
plot(tout,na)
=#
#rho0.data[1,1] = 0.9
#rho0.data[2,2] = .0
#rho0.data[3,3] = .1

## calculated 1st order correlation function and associated spectrum
# STILL UNSURE ABOUT FREQUENCY AXIS/CONVENTION
corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)

## zeropad corr and extrapolate tlist
zp = 0
corr = cmds.zeropad(corr,zp)
tnew, ~ = cmds.interpt(tlist,zp)
ω,spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)

## add the results to previous figure
ax4 = subplot(336)
plot(tnew,real(corr),"gray")
title("System correlation function")
xlabel("time"); ylabel("Corr. func. ⟨ ... | ... ⟩")

ax5 = subplot(339)
plot(-[E₁, E₁], [0, .5])
plot(-[E₂, E₂], [0, .5])
plot(ω,spec,"gray")
plot(-[energies[2], energies[2]], [0, .5],"k--")
plot(-[energies[3], energies[3]], [0, .5],"k--")
title("System abs. spectrum")
xlabel("ω"); ylabel("Absorption")


ωC = 1
S = 0.1
# check out: physics.stackexchange.com/questions/474313/ohmic-spectral-density
function spectral_density(ω)
    if ω == 0.0
        return 0
    else
        return #complex(ω)^S / ωC^(S-1) * exp(-ω/ωC);
    end
end

γ1 = 10
function noise_power(ω)
    if ω == 0.0
        return 0
    else
        return 1 / π * γ1 * ωC ./ (ω.^2 .+ ωC^2)
    end
end

λ = 2.2; γ = 22.1;
#check: home.uchicago.edu/~tokmakoff/TDQMS/Notes/7._Fluctuations_3-09.pdf
function debye(ω)
    if ω == 0.0
        return 0
    else
        return 2 * λ * ω * γ ./ (ω.^2 .+ γ^2)
    end
end

### HOW TO GET RIGHT OPERATORS ??????????
### ????????????????????????????????????
a_ops = [Γ_a[1] * one(b_fock) ⊗ one(b_mon) ⊗ (j12*j21) , noise_power,
         Γ_a[2] * one(b_fock) ⊗ (j12*j21) ⊗ one(b_mon) , noise_power,
         Γ_a[3] * one(b_fock) ⊗ (j21*j12) ⊗ (j21*j12)  , noise_power,
         Γ_a[4] * (at*a)                                 , noise_power]

#a_ops = [*(Γ_a[i], a_ops[2*i-1]) for i in 1:length(Γ_a)]

R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops], J=L)
# different elements in R ... https://www.pnas.org/content/pnas/108/52/20908.full.pdf
# changing R[1,1] let's 2D signal disappear after a certain time T
#R.data[1,1] = -0.01
# however ρ[1,1] goes to zero ... shouldn't it go back to one ???


tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0,R,H)
corr = expect(μ,rhot)

zp = 11
corr = cmds.zeropad(corr,zp)
tnew, ~ = cmds.interpt(tout,zp)

#subplot(336)
ax4.plot(tnew,real(corr),"g",linewidth=1,label="RF")

ω,spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)
ax5.plot(ω,spec)

#subplot(236)
#ax5.plot(ω,spec,"g",linewidth=1,label="RF")

tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0*μ12,R,H)

es1 = expect(one(b_fock) ⊗ (j21*j12)  ⊗ one(b_mon), rhot)
es2 = expect(one(b_fock) ⊗ one(b_mon) ⊗ (j21*j12),  rhot)
cav = expect(create(b_fock) * destroy(b_fock) ⊗ one(b_mon) ⊗ one(b_mon), rhot)
gs  = expect(rho0, rhot)

subplot(333)
plot(tout,real(es1),tout,real(es2),tout,cav,tout,gs)
legend(["ES1", "ES2", "cav", "GS"],loc="right")
tight_layout() 

"""
figure(figsize=(6,4))
ww = [0:.1:100;]
plot(ww,debye(ww));
"""

method = "redfield"
if method == "redfield"
    F       = R
    use_sub = false
elseif method == "lindblad"
    F       = L
    use_sub = true
end

if calc_2d

        ## calculate (complex) 3rd order corr function (with T=0)
        zp = 11 # zeropad up to 2^zp

        ## calculate 2D spectra at
        T = [0,20] #fs

        out2d = Array{cmds.out2d}(undef, length(T))

        Threads.@threads for i = 1:length(T)
            out2d[i] = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                                method;debug=true,use_sub=use_sub,zp=zp);
        end

        ## crop 2D data and increase dw
        out2d = [cmds.crop2d(out2d[i],1.3;w_max=3.5,step=1) for i = 1:length(T)]

        ## plot 2D spectra for each(?) T
        # what to plot
        rep="absorptive"
        scal="lin"

        ## make  subplot layout
        nplots = length(T);
        ncols = Int32(ceil(sqrt(nplots)));
        nrows = Int32(ceil(nplots / ncols));

        # determine maximum value in dataset out2d[:].full2d[:,:]
        maxi = maximum([maximum(real(out2d[i].full2d)) for i in 1:length(out2d)])

        # plot 2D spectra
        fig, ax = subplots(nrows,ncols,sharex=true,sharey=true,figsize=(ncols*3.2,nrows*3),squeeze=false)
        fig.suptitle(rep * " 2D spectrum (" * scal * ". scaling)")
        k = 0
        for i = 1:ncols
            for j = 1:nrows
                global k += 1
                if k > nplots
                    continue
                end
                sca(ax[k])
                ax[k].set_aspect = "equal"
                cmds.plot2d(out2d[k].ω,round.(out2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi)
                ax[k].set_title("2D spectrum at $(T[k]) fs")
            end
        end
        tight_layout()
        subplots_adjust(top=0.88)
        

        ω = out2d[1].ω
        ## plot additional things, like energy levels of states
        plot([E₁, E₁], [ω[1], ω[end]],"k--",linewidth=1,alpha=0.25)
        plot([E₂, E₂], [ω[1], ω[end]],"k--",linewidth=1,alpha=0.25)
        plot([ω[1], ω[end]],[E₁, E₁],"k--",linewidth=1,alpha=0.25)
        plot([ω[1], ω[end]],[E₂, E₂],"k--",linewidth=1,alpha=0.25)
        plot([energies[2:3], energies[2:3]], [ω[1], ω[end]],"g--",linewidth=1,alpha=0.25)
        plot([ω[1], ω[end]],[energies[2:3], energies[2:3]],"g--",linewidth=1,alpha=0.25)

        ## Save data
        cmds.save_2d([round.(real(out2d[i].full2d),digits=1) for i = 1:length(T)],T,fn)

        ## plot TA (summed 2D spectrum)
        figure()
        ta = [sum(real(out2d[i].full2d),dims=1) for i in 1:length(out2d)] ./ length(out2d)
        plot(ω,vcat(ta...)')
        plot(ω,zeros(size(ω)),linestyle = "dashed")
        xlabel("Energy/frequency")
        ylabel("Difference absorption")
        tight_layout()

end
