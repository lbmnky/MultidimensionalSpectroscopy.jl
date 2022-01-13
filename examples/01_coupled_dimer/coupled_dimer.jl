using MultidimensionalSpectroscopy, PyPlot, QuantumOptics, LinearAlgebra, FFTW, Colors,
        Printf, DelimitedFiles

# make sure to set script directory as pwd()
cd(@__DIR__)

# to use PyPlot in GUI mode under Linux
pygui(true)

# run once with "calc_2d = false" to initialize functions
calc_2d = false

cmp = create_colormap("bright");
##

## define functions
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
# here, the point x₁, y₁ = (0,0) is the origin for point dipole 1. The dipole has an angle α₁
# where for 0 the dipole lies along the x-axis and for 90 along the y-axis. It has a length d₁.
# The point dipole number 2 is located at a distance r under the angle β, and has its own dipole
# angle α₂ and a length d₂.

# position of dipole 1 and angle wrt to x-axis (input)
x₁, y₁, α₁, d₁ = 0, 0, angle2rad(0), 0.25
# r: distance, β: slip angle
r , β          = 100.155, angle2rad(23.147495)
#r , β          = 19.3155, angle2rad(0)
# position of dipole 2 (calculated from above)
x₂, y₂         = x₁ + r*cos(β), y₁ + r*sin(β)
# angle and moment of dipole 2 (input)
α₂, d₂         = angle2rad(0), 0.25

## set up energies and couplings
# homo dimer    : E1 = E2
# hetero dimer  : E1 ≠ E2
E₁ = 2.396    # energy eV -> t = hbar /eV  , 0.66 fs
E₂ = 2.796
#E₂ = 2.396
E_RWA = E₁
E_RWA = 0
E₁ = E₁ - E_RWA
E₂ = E₂ - E_RWA

## draw dimer
figure(figsize=(10,4));
ax1 = subplot(231);
grid(b=1,which="both")
plot([x₁,x₂], [y₁,y₂], linestyle="--", color="k")
draw_dipole(x₁,y₁,α₁,d₁); draw_dipole(x₂,y₂,α₂,d₂)
xlim([-0.5, x₂+.5]), ylim([-0.5, y₂+0.5])
title("Dimer geometry")
xlabel("X"); ylabel("Y")

## determine x and y-components of tdm 1 and 2
d₁xy = [cos(α₁)*d₁, sin(α₁)*d₁]         # x and y component of point dipole 1
d₂xy = [cos(α₂)*d₂, sin(α₂)*d₂]         # x and y component of point dipole 2

## determine coupling strength (will depend on geometry)
# first cos(0) in κ = ... corresponds to out-of-plane (computer screen) rotation
κ = 1 * (cos(0) - 3 * cos(β) * cos(β))
J = κ * (d₁*d₂) / norm([x₁,y₁]-[x₂,y₂])^3

## alternative way to calculate coupling constant J
# from Cho's book ("Two-Dimensional Optical Spectroscopy", 2009, p.170)
R12 = [x₁,y₁]-[x₂,y₂]
J = dot(d₁xy,d₂xy) / norm(R12)^3 -
        3 * dot(d₁xy,R12) * dot(R12,d₂xy) / norm(R12)^5

## make Hilbert-space
b_mon = NLevelBasis(2)          # Hilbert-space of monomer
b_dim = b_mon ⊗ b_mon          # Hilbert-space of dimer

## create transition operators in monomer basis
j12 = transition(b_mon,1,2)     # == g_e, or a / σ₋ (annihilation)
j21 = transition(b_mon,2,1)     # == e_g, or a† /σ⁺ (creation)
                                # g_e / e_g notation from Constantin's Matlab script

## convert transition operators into dimer basis
mu1hat = (j12 + j21) ⊗ one(b_mon)
mu2hat = one(b_mon) ⊗ (j12 + j21)

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

H₁ = E₁ * j21 * j12                               # dipole 1
H₂ = E₂ * j21 * j12                               # dipole 2

Hₓ = J * (j12 ⊗ j21 + j21 ⊗ j12);              # Hₓ describes the exchange between dipoles 1 and 2
                                                 # use the following with TDMs between ground- and biexciton state, this modifies
                                                 # the energies of the ground- and doubly excited state slightly

H = one(b_mon) ⊗ H₁ + H₂ ⊗ one(b_mon) + Hₓ      # H is the total Hamiltonian of the system

## get eigenstates and eigenvalues of H
#  using Julia LinearAlgebra functions instead of QuantumOptics requires to use H.data
states   = eigvecs(dense(H).data)
energies = eigvals(dense(H).data)

#TODO: do I need this ?
## perform bases transformation of μ
#  for this I need to express the matrices as dense, data matrices
#μx =  states * μx * states
#μy =  states' * dense(μy) * states

## plot site basis and exciton basis energy diagrams
ax2 = subplot(132)
hlines([0, E₁], -1.25, -0.75)
hlines(energies[:],-0.25, 0.25)
hlines([0, E₂], 0.75, 1.25)
xticks([-1, 0, 1], ["monomer 1", "dimer", "monomer 2"])
title("Mon. and dimer energy levels")

## make initial density matrix
Psi0 = nlevelstate(b_mon,1) ⊗ nlevelstate(b_mon,1)
rho0 = dm(Psi0)

## time list for evaluation of master equation
tlist = [0:0.7:400;]

## make collapse operator #CHECK #TODO: understand better (why ccc?)
ccc = diagonaloperator(b_mon,[1, 0])
L1 = one(b_mon) ⊗ j12 #
L1 = ccc⊗j12 # ... ccc don't make sense
L2 = j12 ⊗ one(b_mon) #
L2=j12⊗ccc

Γ = [sqrt(0.025), sqrt(0.025)]
#L = [sqrt(0.25)*one(b_mon)⊗j12, sqrt(0.25)*j12⊗one(b_mon)]
L = Γ .* [L1, L2]

L = [sqrt(0.055) * (j21*j12 ⊗  (j12*j21)+ (j12*j21) ⊗  (j21*j12)), Γ[2] * L2]
#DEBUG
#display(dense(c_ops[1])); display(dense(c_ops[2]))

## make transition dipole operator
# clearly now μ has a direction dependence, but if sample is isotropic take:
μ   = (μx + μy)     #CHECK
# TDM with the ground state:
μ12 = μ .* (j12*j21⊗j21 + j12*j21⊗j12 + j21⊗(j12*j21) + j12⊗(j12*j21))
#DEBUG display(dense(μ12))
# ... and with the doubly excited state:
μ23 = μ .* (j21*j12⊗j21 + j21*j12⊗j12 + j21⊗(j21*j12) + j12⊗(j21*j12))
# DEBUG display(dense(μ23))

## Normalize TDM operators #TODO: Is this the best way to do this ?
rho1 = μ12 * rho0 * μ12
μ12  = μ12 / sqrt(tr(rho1))
μ23  = μ23 / sqrt(tr(rho1))

rho1 = μ12 * rho0 * μ12
rho2 = μ23 * rho1 * μ23
μ23  = μ23 / sqrt(tr(rho2))
rho2 = μ23 * rho1 * μ23

## calculate magnitude of tdm with singly excited states e1 and e2
# HOT TO ORDER THEM ???? ... depending on J>0 (H) or J<0 (J) ? Something more
# fundamental possible ? States not order from diagonalization (?)
if J < 0
        # tdm magnitude <0,0|μ|0,1>
        a = states[1,:]' * dense(μ).data * states[2,:]
        # tdm magnitude <0,0|μ|1,0>
        b = states[1,:]' * dense(μ).data * states[3,:]
elseif J > 0
        # everything reversed
        a = states[1,:]' * dense(μ).data * states[3,:]
        b = states[1,:]' * dense(μ).data * states[2,:]
end

## plot TDM strengths
ax3 = subplot(234)
bar(0, real(a))
bar(1, real(b))
xticks([0, 1], ["μ_g-e1", "μ_g-e2"])
xlim([-1, 2]); ylim([0, 1.5 * maximum(real([a, b]))])
ylabel("tdm strength")

#DELETE (?)
## calculate and plot expectation values of population
#=
tout, rhot = timeevolution.master(tlist, rho0, H, F)
nc = real(expect(j12 ⊗ j21, rhot))
na = real(expect(j21 ⊗ j12, rhot))
figure(figsize=(6,3));
plot(tout,nc)
plot(tout,na)
=#

## calculated 1st order correlation function and associated spectrum
corr = timecorrelations.correlation(tlist, rho0, H, L, μ12, μ12)

## zeropad corr and extrapolate tlist
zp = 10
corr = zeropad(corr,zp)
tnew, ~ = interpt(tlist,zp)

## convert to spectrum
ω,spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)
ω = ω .- E_RWA
## add the results to previous figure
# correlation function
ax4 = subplot(233)
plot(tnew,real(corr),"k")
title("Dimer correlation function")
xlabel("time"); ylabel("Corr. func. ⟨μ(t) μ₀₁ ρ₀⟩")
# spectrum
ax5 = subplot(236)
plot(-[E₁, E₁], [0, .5])
plot(-[E₂, E₂], [0, .5])
plot(ω,spec,"k")
plot(-[energies[2], energies[2]], [0, .5],"k--")
plot(-[energies[3], energies[3]], [0, .5],"k--")
title("Dimer abs. spectrum")
xlabel("ω"); ylabel("Absorption")
tight_layout()

## Redfield
"""
ωC = 10
S = .005
# check out: physics.stackexchange.com/questions/474313/ohmic-spectral-density
function spectral_density(ω)
    if ω == 0.0
        return 0
    else
        return complex(ω).^S / ωC^(S-1) .* exp.(-ω/ωC);
    end
end
"""

"""
ωC = 10
γ1 = 0.1
function noise_power(ω)
    if ω == 0.0
        return 0
    else
        return 1 / π * γ1 * ωC ./ (ω.^2 .+ ωC^2)
    end
end
"""
τc = 1
Δ = .01
temp = 0.3
Λ = 1/τc;
λ = Δ^2 / temp
λ = .2
#CHECK: home.uchicago.edu/~tokmakoff/TDQMS/Notes/7._Fluctuations_3-09.pdf

function debye(ω) #TODO: understand better
    if ω == 0.0
        #return 0
        return λ * temp / Λ
    else
        spec = 2 * λ * ω * Λ ./ ( ω.^2 .+ Λ^2 )
        return ( 1 .+ coth.(ω / 2 / temp) ) .* spec
        #return λ * (2 * Λ * ω ./ (ω.^2 .+ Λ^2))
        #return -1 * (λ * 2 * ω ./ (ω.^2 .+ Λ^2))
        #return -.4 / (ω .+ .02)
    end
end

### HOW TO GET RIGHT OPERATORS ??????????
### ????????????????????????????????????
#a_ops = [sqrt(0.1)*((j12*j21)⊗one(b_mon)+one(b_mon)⊗(j12*j21)),spectral_density]
a_ops = [sqrt(.32)*one(b_mon)⊗(j12*j21), debye, sqrt(.32)*(j12*j21)⊗one(b_mon),debye]
a_ops = [(j12*j21)⊗ j12 + j12⊗(j21*j12), debye, (j21*j12)⊗ j12 + j12⊗(j12*j21), debye]
#L = [sqrt(0.01)*j12⊗one(b_mon),sqrt(0.01)*one(b_mon)⊗j12]
#L = [sqrt(.07)*one(b_mon)⊗(j12*j21), sqrt(.07)*(j12*j21)⊗one(b_mon)]
R, ekets = timeevolution.bloch_redfield_tensor(H, [a_ops]; J=L)
# different elements in R ... https://www.pnas.org/content/pnas/108/52/20908.full.pdf
# changing R[1,1] let's 2D signal disappear after a certain time T
#R.data[1,1] = -0.01
# however ρ[1,1] goes to zero ... shouldn't it go back to one ???

# rate of a_ops affects spectral width
# rate of L (Lindblad-operators) affects t2 evolution

## calculate master equation time evolution
tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0,R,H)
# and the correlation function
corr = expect(μ,rhot)

## zeropad data for smoothing
zp = 11
corr = zeropad(corr,zp)
tnew, ~ = interpt(tout,zp)

## plot correlation function from Redfield
ax4.plot(tnew,real(corr)./maximum(real(corr)),"g",linewidth=1,label="RF")

## calculate spectrum
ω,spec = timecorrelations.correlation2spectrum(tnew, corr; normalize_spec=true)
# and plot
subplot(236)
ax5.plot(ω,spec,"g",linewidth=1,label="RF")

## calculated populations
tout, rhot = timeevolution.master_bloch_redfield(tlist,μ12*rho0*μ12,R,H);


es1 = expect((j21*j12)⊗one(b_mon),rhot)
es2 = expect(one(b_mon)⊗(j21*j12),rhot)
# and plot
ax4.plot(tout,real(es1),tout,real(es2))

#= #DEBUG bath spectral function
figure(figsize=(6,4))
ww = [0:.1:100;]
plot(ww,debye(ww));
=#


# select method for 2D
method = "redfield"
if method == "redfield"
    F       = R
    use_sub = false
elseif method == "lindblad"
    F       = L
    use_sub = false
end

if calc_2d

        zp = 11 # zeropad up to 2^zp

        ## calculate 2D spectra at
        T = [0, 50] #fs

        spectra2d = Array{out2d}(undef, length(T))

        Threads.@threads for i = 1:length(T)
        #for i = 1:length(T)
            spectra2d[i] = make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                                method;debug=true,use_sub=use_sub,
                                                    t2coh="kin",zp=zp);    #IDEA: can check 2D electronic beating
        end

        ## crop 2D data and increase dw
        spectra2d = [crop2d(spectra2d[i],.5 .- E_RWA;w_max=5 .-E_RWA,step=1) for i = 1:length(T)]

        ## plot 2D spectra for each(?) T
        # what to plot
        rep="absorptive"
        scal="asinh"

        ## make  subplot layout
        nplots = length(T);
        ncols = Int32(ceil(sqrt(nplots)));
        nrows = Int32(ceil(nplots / ncols));

        # determine maximum value in dataset out2d[:].full2d[:,:]
        maxi = maximum([maximum(real(spectra2d[i].full2d)) for i in 1:length(spectra2d)])
        ω = spectra2d[1].ω .+ E_RWA
end

if calc_2d

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
                #plot2d(out2d[k].ω,round.(out2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi/2)
                plot2d(ω,round.(spectra2d[k].full2d,digits=1);repr=rep,scaling=scal,norm=maxi/2)
                title("2D spectrum at $(T[k]) fs")
            end
        end
        tight_layout()
        subplots_adjust(top=0.9)

        ## plot additional things, like energy levels of states
        plot([E₁, E₁], [ω[1], ω[end]],"k--",linewidth=1,alpha=0.25)
        plot([E₂, E₂], [ω[1], ω[end]],"k--",linewidth=1,alpha=0.25)
        plot([ω[1], ω[end]],[E₁, E₁],"k--",linewidth=1,alpha=0.25)
        plot([ω[1], ω[end]],[E₂, E₂],"k--",linewidth=1,alpha=0.25)
        plot([energies[2:3], energies[2:3]], [ω[1], ω[end]],"g--",linewidth=1,alpha=0.25)
        plot([ω[1], ω[end]],[energies[2:3], energies[2:3]],"g--",linewidth=1,alpha=0.25)

        ## Save data
        #save_2d(out2d,T,fn)

        ## plot TA (summed 2D spectrum)
        figure()
        ta = [sum(real(spectra2d[i].full2d),dims=1) for i in 1:length(spectra2d)] ./ length(spectra2d)
        plot(ω,vcat(ta...)')
        plot(ω,zeros(size(ω)),linestyle = "dashed")
        xlabel("Energy/frequency")
        ylabel("Difference absorption")
        tight_layout()

end

#
#### use slider for flipping through 2D spectra
#
#using Blink, Interactive

#=
mp = @manipulate for i in slider(1:length(spectra2d))
          clf();cmds.plot2d(spectra2d[i].ω,spectra2d[i].full2d)
          end
w = Window(); body!(w, mp);
=#
