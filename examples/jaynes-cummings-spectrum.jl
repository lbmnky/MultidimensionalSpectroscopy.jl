#!/usr/bin/julia
# using PyPlot
# using QuantumOptics

# make sure to set script directory as pwd()
cd(@__DIR__)

# include my custom cmds module
if Sys.iswindows()
    include("..\\cmds.jl")
    fn = "01_Output"
else
    include("../cmds.jl")
    fn = "01_Output"
end


pygui(true)

calc_2d = true

cmp = cmds.create_colormap("bright");


wc = 2.0  * 2 * pi  # cavity frequency
wa = 2.0  * 2 * pi  # atom frequency
g  = 0.25 * 2 * pi   # coupling strength
κ = 0.5         # cavity dissipation rate
γ = 0.0005         # atom dissipation rate
N = 1               # number of cavity fock states
n_th_a = 0.0        # avg number of thermal bath excitation
#use_rwa = True

tlist = [0:0.1:50;]
b_fock = FockBasis(N)
b_atom = NLevelBasis(3)
Psi0 = fockstate(b_fock,1) ⊗ nlevelstate(b_atom,1)    # start with an excited cavity
println(Psi0)

a = destroy(b_fock) ⊗ one(b_atom)
at = dagger(a)
println(a)
sm = one(b_fock) ⊗ transition(b_atom,1,2)
sp = dagger(sm)
println(sm)

H = wc * at * a + wa * sp * sm + g * (at * sm + a * sp)

println(H)

#J = [sqrt(κ)*sm,sqrt(γ)*a]
# the same ?
J = [sm,a]
rates = [sqrt(κ),sqrt(γ)]

tout, rhot = timeevolution.master(tlist,Psi0,H,J;rates=rates)
nc = real(expect(at * a, rhot))
na = real(expect(sp * sm, rhot))

figure(figsize=(6,3))
plot(tout,nc)
plot(tout,na)

#ω, spec = timecorrelations.correlation2spectrum(tlist, nc; normalize_spec=true)
#figure(figsize=(6,3))
#plot(ω, spec)

println(rhot[1])
println(rhot[end])
println(J)
rho0 = dm(Psi0)

corr = timecorrelations.correlation(tlist, rho0, H, J, at, a; rates=rates)
ω, spec = timecorrelations.correlation2spectrum(tlist, corr; normalize_spec=true)
plot(tlist,corr)

figure(figsize=(6,3))
plot(ω, spec)

display(rhot[end])

Psi0 = fockstate(b_fock,0) ⊗ nlevelstate(b_atom,1)    # start with system in gs
rho0 = dm(Psi0)


μ12 = at * sm*sp + a * sm*sp
μ23 = at * sp*sm + a * sp*sm

if calc_2d
        ## calculate (complex) 3rd order corr function (with T=0)
        zp = 10 # zeropad up to 2^zp

        F = J

        ## calculate 2d spectra at
        T = [0., 10., 20., 30]
        #out2d = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T,"lindblad";debug=false,zp=zp);
        out2d = Array{cmds.out2d}(undef, length(T))
        for i = 1:length(T)
            out2d[i] = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],
                                                "lindblad";debug=false,zp=zp);
        end

        ## crop 2D data and increase dw
            out2d = [cmds.crop2d(out2d[i],3.5;w_max=20,step=2) for i = 1:length(T)]

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

show();
