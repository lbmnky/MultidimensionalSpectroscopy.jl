module MultidimensionalSpectroscopy

export create_colormap, zeropad, interpt, make2Dspectra, correlations,
        view_dm_evo, save_2d, load_2d, plot2d, crop2d, round2d, tri, absorptionSpectrum, plot2d_comps,
         pretty_show_mat, vib_analysis, create_subspace, rand_normal, plot_levels, logger, out2d

using QuantumOptics, FFTW, LinearAlgebra, PyPlot, Colors, DelimitedFiles, Printf, Random, JLD2, Interact, Blink

"""
    plot2d(ω, data; repr = "absorptive", norm_spec=false, scaling="lin")

Plots the 2D data.

# Arguments
ω         : frequency / energy axis
data      : data structure of 2D output
repr      : pick a representation ("absorptive", "dispersive", "absolute")
norm_spec : normalize to abs.(maximum(data)), where data is the spectrum to
            be plotted
scaling   : pick z-scaling ("lin", "tan", "cube", "asinh")
norm      : if norm == 0 normalize each 2D spectrum individually, else normalize
            all spectra to norm = <val>

"""
function plot2d(ω, data; repr = "absorptive", norm_spec=false, scaling="lin", norm=0)

    ## select 2D spectrum to plot
    if repr in ["absorptive", "abst", "real"]
        data = real(data)
        cmp = create_colormap()
    elseif repr in ["dispersive", "disp", "imaginary", "imag", "i"]
        data = imag(data)
        cmp = create_colormap()
    elseif repr in ["absolute", "absl"]
        #data = sqrt.(real(data).^2 + imag(data).^2) #BUG doing this calculation results in Inf values when data{ComplexF16}
        data = abs.(data)
        cmp = create_colormap()
    elseif repr in ["phase", "ph"]
        data_cntl = real(data) ./ maximum(abs.(data))
        data = atan.(imag(data) ./ real(data)) 
        cmp = "rainbow"
    end

    ## normalize 2D spectrum to plot
    if norm_spec
        data = data/maximum(abs.(data));
    end

    sensitivity = "++"

    ## make levels for contourplot with different scaling
    lvls = [-1.025:0.05:1.025;]
    if sensitivity == "+++"
        lvls = [-1,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,-.08,-.06,-.04,-.02,-.01,-0.001, 0,
                .001,.01,.02,.04,.06,.08,.1,.2,.3,.4,.5,.6,.7,.8,.9,1];
    elseif sensitivity == "++"
         lvls = [-1,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,-.08,-.06,-.04,-.02, 0,
                .02,.04,.06,.08,.1,.2,.3,.4,.5,.6,.7,.8,.9,1];
    elseif sensitivity == "+"
        lvls = [-1,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,-.075,-.05, 0,
                .05,.075,.1,.2,.3,.4,.5,.6,.7,.8,.9,1];
    else
    end

    if scaling == "lin"
        #
    elseif scaling == "tan"
        lvls = tan.(lvls); lvls = lvls ./ maximum(lvls);
    elseif scaling == "cube"
        lvls = lvls.^3;
    elseif scaling == "asinh"
        lvls = asinh.(0.1*lvls); lvls = lvls ./ maximum(lvls);
    end


    # rescale lvls to data (don't rescale data, takes more time)
    m = -1*(maximum(abs.(data))); M = (maximum(abs.(data)))
    if m == 0 && M == 0
        lvls = 0.001 * lvls;
    elseif m == Inf || M == Inf
        println("\nERROR: z-value is Inf!\n")
    else
        lvls = lvls * maximum(abs.([m, M]));
        lvls_ticks = [-1:0.2:1;] * maximum(abs.([m, M]))
    end

    # normalize to global maximum with evolution time T scan
    if norm == 0
    else
        lvls = lvls ./ maximum(lvls) * norm;
    end

    
    ## plot data as filled contour
    conv = "ω3ω1"
    if conv == "ω1ω3"
        cs  = contourf(ω,ω,data,lvls,cmap=create_colormap()); colorbar();
        if repr == "phase"
            cs2 = contour(ω,ω,data_cntl,lvls,colors="white",linewidths=(.5,))
        else
            cs2 = contour(cs,levels=cs.levels[:],colors="gray",linewidths=(.5,))
        end
        xlabel("detection (ω₃)"); ylabel("excitation (ω₁)"); #clim([m, M]);
        
    elseif conv == "ω3ω1"
        cs  = contourf(ω,ω,transpose(data),lvls,cmap=cmp);
        if repr == "phase"
            cs2 = contour(ω,ω,transpose(data_cntl),lvls,colors="white",linewidths=(.5,))
        else
            cs2 = contour(cs,levels=cs.levels[:],colors="gray",linewidths=(.5,))
        end
        xlabel("excitation (ω₁)"); ylabel("detection (ω₃)");
    end
    cbar = colorbar(cs)
    cbar.add_lines(cs2)
    plot([ω[1], ω[end]], [ω[1], ω[end]],"k--");

    ## plot positive and negative contourlines with different linestyles
    #contour(ω,ω,temp.*(temp.>0),lvls[1:1:end],cmap="Greys",
    #        linewidths=0.5,linestyles="solid");
    #contour(ω,ω,temp.*(temp.<0),lvls[1:1:end],cmap="Greys",
    #        linewidths=0.5,linestyles="dashed");

end

function plot_timeTrace(dat2d,T,ω,w1,w3)
    figure(figsize=(8.5,3.5))
    subplot(121)
    plot2d(ω,dat2d[1])
    for i in 1:length(w1)
        plot(w1[i],w3[i],"o",fillstyle="none",mew=3,ms=10,marker="s")
    end
    subplot(122)
    idx1 = [argmin(abs.(ω .- i)) for i in w1]
    idx3 = [argmin(abs.(ω .- i)) for i in w3]
    for j in 1:length(idx1)
        I = [dat2d[i][idx1[j],idx3[j]] for i in 1:length(T)]
        plot(T,I,"o-",label=string(w1[j])*";"*string(w3[j]))
    end
    xlabel("Time"); ylabel("Intensity at w1;w3")
    legend()
    tight_layout()
    #sel = slider(1:length(dat2d))
    #w = Window()
    #Interact.@on 
    #body!(w, sel);
end

"""
    save_2d(spec2d,T,fn_base)

Saves 2D spectra to specified directory. Can be used in a loop to store all
data.

# Arguments
- spec2d  : 2D array of data to be saved
- T       : array of T times
- fn_base : directory to store data in

# Example
```julia
    cmds.save_2d([out2d[i].full2d for i = 1:length(T)],T,fn)
```
"""
function save_2d(spec2d,T,fn_base)
    if isdir(fn_base)
        """ nothing """
    else
        mkdir(fn_base)
    end
    writedlm(fn_base * "T_steps.dat"    , T, ',')
    writedlm(fn_base * "energy_axis.dat", spec2d[1].ω, ',')
    for i = 1:length(T)
        fn = fn_base * @sprintf("2Dspec_%03i.dat",i)
        writedlm(fn, round.(real.(spec2d[i].full2d), digits=2), ',')
        fn = fn_base * @sprintf("2Dspec_%03i_gsb.dat",i)
        writedlm(fn, round.(real.(spec2d[i].gsb), digits=2), ',')
        fn = fn_base * @sprintf("2Dspec_%03i_se.dat",i)
        writedlm(fn, round.(real.(spec2d[i].se), digits=2), ',')
        fn = fn_base * @sprintf("2Dspec_%03i_esa.dat",i)
        writedlm(fn, round.(real.(spec2d[i].esa), digits=2), ',')
        println("\nFile saved as: " * fn)
    end
    println("\nFiles saved!")
end

function save_2d(spec2d,T,fn_base,type)
    if isdir(fn_base)
        """ nothing """
    else
        mkdir(fn_base)
    end
    #TODO: correctly save struct to .h5
    fn = fn_base * "2Dspec.h5"
    h5open(fn, "w") do file
        write(file, "timesteps", T)
        for i = 1:length(T)
            write(file, @sprintf("2Dspec_%03i",i), spec2d[i])  # alternatively, say "@write file A"
        end
    end

    println("\nFiles saved!")
end

"""
    load 2D spectra
"""
function load_2d(fn_base;type="full")
    #
    T = readdlm(fn_base * "T_steps.dat")
    ω = readdlm(fn_base * "energy_axis.dat")
    #
    dat2d = []
    for i in 1:length(T)
        if type == "full"
            push!(dat2d,readdlm(fn_base * @sprintf("2Dspec_%03i.dat",i),','));
        end
    end
    return T, ω, dat2d;
end

"""
    view_dm_evo(rhot,dt)

Quick function to view evolution of density matrix elements as 2D plot.

# Arguments
rhot : time evolution of density matrix
dt   : time steps between subplots

"""
function view_dm_evo(rhot,dt)
    figure()
    for i = 1:12
        subplot(3,4,i,aspect="equal")
        j = (i-1) * dt + 1
        pcolormesh(real(rhot[j].data),cmap="seismic"); clim([-1,1]);
        xticks([1:1:9;]); yticks([1:1:9;])
        grid("on",which="major")
        #title("tlist # $j")
    end
    tight_layout()
end

## return a random sample from a normal (Gaussian) distribution
# https://www.johndcook.com/blog/services-2/
function rand_normal(mean, stdev)
    if stdev <= 0.0
        error("standard deviation must be positive")
    end
    u1 = rand()
    u2 = rand()
    r = sqrt( -2.0*log(u1) )
    theta = 2.0*pi*u2
    mean + stdev*r*sin(theta)
end

"""
    cmp = create_colormap(scheme)

Create a "bright" or "dark" custom colormap to display 2D spectra.
"""
function create_colormap(scheme="bright")

    ## define colors from most negative (MMM) ... to ... most positive (PPP)
    MMM = RGB(0.718 , 0     , 0.718 );
    MM  = RGB(0.516 , 0.516 , 0.991 );
    M   = RGB(0     , 0.559 , 0.559 );
    P   = RGB(0     , 0.592 , 0     );
    PP  = RGB(0.827 , 0.827 , 0     );
    PPP = RGB(0.947 , 0.057 , 0.057 );

    ## zero value is either "bright" or "dark"
    if scheme == "bright"
        Z = RGB(.93,.93,.93);
    elseif scheme == "dark"
        Z = RGB(.23,.23,.23);
    else
        error("cmds -- no colormap scheme selected")
    end

    return cmp  = ColorMap("cmps",[MMM,MM,M,Z,P,PP,PPP]);
end


"""
    dat = zeropad(dat, N)

Add zeros to end end of correlation data (after it has decayed to zero) in order
to increase frequency/energy resolution in the spectral domain (i.e. after Fourier
transform).

# Arguments
* 'dat' : correalation data (1D or 2D)
* 'N'   : data will be zeropadded up to 2^N, if N = 0 return input data

"""
function zeropad(dat, N)

    if 0 < N < log2(size(dat,1))
        # if N smaller than log2(size of data), increase N by 1
        N = ceil(log2(size(dat,1)))
    elseif N == 0
        # if N == zero, return original data
        println("\nNo zero padding.\n")
        return dat;
    end

    # Sanity check that N is not too large
    if N > 2^12
        error("Large N detected. Check if N isn't too large.")
    end

    if N > 0
        # switch between 1D and 2D array
        if size(dat,2) == 1
            out = complex(zeros(2^N));
            println("\nZeropadding 1D-data to $(2^N) values\n")
        elseif size(dat,2) > 1
            out = complex(zeros(2^N,2^N));
            println("\nZeropadding 2D-data to $(2^N) x $(2^N) values\n")
        end

        out[1:size(dat,1),1:size(dat,2)] = dat;

        return out;
    end
end


"""
    tlist = interpt(tlist, N)

Interpolate frequency/energy vector after zeropadding of data using "zeropad".

# Arguments
* 'tlist'   :
* 'N'       :

# Output
* 'tlist'   : extrapolated time vector
* 'ω'       : interpolated frequency/energy vector
"""
function interpt(tlist, N)

    # if N smaller than log2(size of data), increase N by 1
    if 0 < N < log2(length(tlist))
        N = ceil(log2(size(dat,1)))
    elseif N == 0
        tlist = tlist;
    end

    # Sanity check that N is not too large
    if N > 2^12
        error("Large N detected. Check if N isn't too large.")
    end

    if N > 0
        # extrapolate tlist / get new tlist[end]
        dt = abs(tlist[end]-tlist[1]) / (length(tlist)-1)
        tlist = [0:(2^N)-1;] * dt .+ tlist[1];
    end

    # create ω axis corresponding to zeropadded data
    dω = 1/tlist[end]
    ω = ([1:length(tlist);] .- length(tlist)/2 ) * dω * 2 * π;
    # return zeropadded data, corresponding ω and tlist

    return tlist, ω;
end

"""
    out2d

Structure that stores the output of "make2Dspectra"

# Elements:
* 'full2d'      : complex 2D spectrum, take real() or imag() to get absorptive
                  or refractive 2D spectrum
* 'full2d_r'    : rephasing complex 2D spectrum
* 'full2d_nr'   : non-rephasing complex 2D spectrum
* 'gsb'         : complex 2D spectrum with only GSB pathways
* 'se'          : complex 2D spectrum with only SE pathways
* 'esa'         : complex 2D spectrum with only ESA pathways
"""
mutable struct out2d{T}
    ω
    full2d::T
    full2d_r::T
    full2d_nr::T
    gsb::T
    gsb_r::T
    gsb_nr::T
    se::T
    se_r::T
    se_nr::T
    esa::T
    esa_r::T
    esa_nr::T
    corr
end

"""
    create_subspace(H, manifold, F, dat...)

Function accepts variable number of arguments (dat) in addition to H and returns tuple 
of said arguments in Eigen basis (exciton basis/polariton basis/...).

# Arguments
* 'H'        : Hamiltonian
* 'manifold' : up to which excitation manifold/sector should the subspace extent, "si", "bi", "bi_lowest", "bi_polariton"
* 'F'        : operator list, e.g. list of Lindblad dissipation operators
* 'dat'      : either a single operator, or a tuple of operators, e.g.:

H, ..., dat_a                = create_subspace(H, ..., dat_a )

H, ..., dat_a, dat_b         = create_subspace(H, ..., dat_a, dat_b)

H, ..., dat_a, dat_b, dat_c  = create_subspace(H, ..., dat_a, dat_b, dat_c)

...

# Output
* 'H_sub'       : subspace Hamiltonian
* 'transf_op'   : operator to convert between Eigen and site basis
* 'P'           : projection operator
* 'F_sub'       : operator list, e.g. Lindblad dissipation operators, in subspace
* 'dat_sub'     : single or tuple of operators in subspace

`transf_op` can be used to convert between Eigen basis and site basis, e.g.:

`H_site    = transf_op     * H_eigen     * transf_op'`
"""
function create_subspace(H, manifold, dat...)
    Eivecs = eigvecs(dense(H[1]).data)
    L = length(H[1].basis_l.bases)
    # size of subspace 
    if manifold == "bi"
        n_sub = 1 + L + binomial(L,2)
    elseif manifold == "si"
        n_sub = 1 + L + 1
    end

    # eigenvectors of subspace
    Eivecs_sub = [Ket(H[1].basis_l,Eivecs[:,i]) for i in 1:n_sub]
    # subspace basis
    b_sub = SubspaceBasis(H[1].basis_l,Eivecs_sub)
    # projection operators
    P = sparse(projector(b_sub, (H[1].basis_l)))
    Pt = dagger(P)
    println("dim(superspace): ", *(H[1].basis_l.shape...))
    # output new dimension
    println("dim(subspace): ", length(b_sub))
    H_sub   = [P * H[i] * Pt for i in 1:length(H)]
    dat_sub = [P * dat[i] * Pt for i in 1:length(dat)]
    println(dat)
#    F_sub    = [P * F[i] * Pt for i in 1:length(F)]
    return H_sub, Eivecs, P, dat_sub...;
end

function create_subspace(H, manifold, F, dat...)
    Eivecs = eigvecs(dense(H[1]).data)
    L = length(H[1].basis_l.bases)
    # size of subspace 
    if manifold == "bi"
        #n_sub = 1 + L + binomial(L,2)
        n_sub = collect(1:1+L+binomial(L,2)+4)   #TODO: +1 for 2nd fock state... w/o +1 only true for all TLSs (cavity and matter)
    elseif manifold == "bi_lowest"
        #n_sub = 1 + L + 1                       # just use a single bi-excitonic state as an approximation
        n_sub = collect(1:1+L+1)
    elseif manifold == "bi_polariton"
        n_sub = append!(collect(1:1+L+1),1+L+binomial(L,2)+1)   # just use bi-excitonic states that are coupled to field ...
    elseif manifold == "si"
        #n_sub = 1 + L
        n_sub = collect(1:1+L)
    end

    # eigenvectors of subspace
    Eivecs_sub = [Ket(H[1].basis_l,Eivecs[:,i]) for i in n_sub]
    # subspace basis
    b_sub = SubspaceBasis(H[1].basis_l,Eivecs_sub)
    # projection operators
    P = sparse(projector(b_sub, (H[1].basis_l)))
    Pt = dagger(P)
    println("dim(superspace): ", *(H[1].basis_l.shape...))
    # output new dimension
    println("dim(subspace): ", length(b_sub))
    H_sub   = [P * H[i] * Pt for i in 1:length(H)]
    dat_sub = [P * dat[i] * Pt for i in 1:length(dat)]
    F_sub   = [P * F[i] * Pt for i in 1:length(F)]

    transf_op = Operator(H[1].basis_l, H[1].basis_r, Eivecs)
    transf_op = P * transf_op * Pt
    #transf_op = Eivecs
    return H_sub, transf_op, P, F_sub, dat_sub...;
end



"""
    out2d = make2dspectra(tlist, rho0, H, F, μ12, μ23, method; debug=false, zp=0)


make2Dspectra.jl calculates 2D spectra (GSB, SE, and ESA) using 2D correlation
functions by evaluating the Lindblad master equation.

# Arguments
* 'tlist' : time steps at which to evaluate correlation functions τ and t
* 'rho0'  : ground state density matrix
* 'H'     : Hamiltonian of the system
* 'F'     : collapse operators (i.e. relaxation pathways during τ and t) OR
            relaxation tensor R for Redfield
* 'μ12'   : transition dipole operator between ground- and excited states,
            responsible for GSB and SE
* 'μ23'   : transition dipole operator between excited- and higher excited states,
            responsible for ESA
* 'method': select either Lindblad or Redfield master equation solver
* 'debug' : if true, plot response functions and spectra of GSB, SE and ESA and
            show intermediate density matrices
* 'zp'    : zero pad data up to 2^zp in both dimensions. Improves data quality
            after fft, whenever timeevolution has already decayed to 0. Remember
            to also scale tlist and ω via function "interpt.jl". Default zp=0,
            data is not zeropadded.

Hidden is a lineshape function (work in progress).
"""
function make2Dspectra(tlist, rho0, H, F, μ12, μ23, T, method; debug=false, use_sub=true, zp=0, t2coh="kin")

    if use_sub 
    H_si    = [[pop!(H)]]
    H_bi    = [[pop!(H)]]
    rho0_si = pop!(rho0)
    rho0_bi = pop!(rho0)
    F_si    = pop!(F)
    F_bi    = pop!(F)
    μ12_si  = pop!(μ12)
    μ12_bi  = pop!(μ12)
    μ23_si  = pop!(μ23)
    μ23_bi  = pop!(μ23)
    else
        H_si = [[H]]
        H_bi = [[H]]
        rho0_si = rho0
        rho0_bi = rho0
        F_si = F
        F_bi = F
        μ12_si = μ12
        μ12_bi = μ12
        μ23_si = μ23
        μ23_bi = μ23
    end

    println(rho0)

    println("############# ########################### #############");
    println("############# calculate R and NR pathways #############");
    println("############# ########################### #############\n");

    #TODO: implement T-dependent H
    N = 1
    # if !isempty(methods(H))
    #     H = [H(T) for i in 1:N]
    #     println.(H) 
    # elseif length(H) == 2
    #     println("L2")
    #     #println(H)
    #     N = length(H[1])
    #     H = [[H[1][i], H[2][i]] for i in 1:N]
    # else
    #     N = 1
    #     H = [[H]]
    # end

    # keep around for now
    if use_sub
        println("Using system subspace!")
        #H_bi, μ12_bi, μ23_bi, rho0_bi, F_bi =  create_subspace(H, μ12, μ23, rho0, F, "bi")
    else
    #    H_bi, μ12_bi, μ23_bi, rho0_bi, F_bi = H, μ12, μ23, rho0, F
    end

    # TODO: can use a 2D lineshape function to cause inhomogeneous broadening, but need to decrease rate
    # of Lindblad operators beforehand ... very narrow transition, then broadened by lineshape function
    # inhomogeneously
    D = .0001
    τc = 20
    tc = 20 #+ T
    #τc = D / 0.000125
    #tc = D / 0.000125
    #gτ = D^2 * τc^2 * (exp.(-tlist./τc) .+ tlist./τc .- 1)
    gτ_(τ) = D^2 * τc^2 * (exp.(-τ./τc) .+ τ./τc .- 1)
    #gt = D^2 * tc^2 * (exp.(-tlist./tc) .+ tlist./tc .- 1)
    gt_(t) = D^2 * tc^2 * (exp.(-t./tc) .+ t./tc .- 1)
    test = true

    # excited state absorption
    corr_func_NR_esa = sum(correlations(tlist,rho0_bi,H_bi[i],F_bi,μ12_bi,μ23_bi,T,
                            "NR_esa",method,false; t2coh=t2coh) for i in 1:N) ./ N
    corr_func_R_esa  = sum(correlations(tlist,rho0_bi,H_bi[i],F_bi,μ12_bi,μ23_bi,T,
                            "R_esa", method,false; t2coh=t2coh) for i in 1:N) ./ N
    #rhots_NR_esa = rhots_NR_esa - conj(rhots_NR_esa)
    #rhots_R_esa  = rhots_R_esa  - conj(rhots_R_esa)
    if test
        corr_func_NR_esa = corr_func_NR_esa .* (exp.(-gτ_(-tlist)) * exp.(-gt_(-tlist))')
        corr_func_R_esa  = corr_func_R_esa  .* (exp.(-gτ_( tlist)) * exp.(-gt_( tlist))')
    end

    # excited state absorption ... from GS #TODO!: Can't use new Ham(T) here, needs to be same as for GSB
    corr_func_NR_esax = sum(correlations(tlist,rho0_bi,H_bi[i],F_bi,μ12_bi,μ23_bi,T,
                                "NR_esax",method,false; t2coh=t2coh) for i in 1:N) ./ N
    corr_func_R_esax  = sum(correlations(tlist,rho0_bi,H_bi[i],F_bi,μ12_bi,μ23_bi,T,
                                "R_esax", method,false; t2coh=t2coh) for i in 1:N) ./ N
    if test
        corr_func_NR_esax = corr_func_NR_esax .* (exp.(-gτ_(-tlist)) * exp.(-gt_(-tlist))')
        corr_func_R_esax  = corr_func_R_esax  .* (exp.(-gτ_( tlist)) * exp.(-gt_( tlist))')
    end

    # keep this around for now
    if use_sub
       # H_si, μ12_si, μ23_si, rho0_si, F_si =  create_subspace(H, μ12, μ23, rho0, F, "si")
    else
    #    H_si, μ12_si, μ23_si, rho0_si, F_si = H, μ12, μ23, rho0, F
    end

    # ground state bleach
    corr_func_NR_gsb = sum(correlations(tlist,rho0_si,H_si[i],F_si,μ12_si,μ23_si,T,
                            "NR_gsb",method,false; t2coh=t2coh) for i in 1:N) ./ N
    corr_func_R_gsb  = sum(correlations(tlist,rho0_si,H_si[i],F_si,μ12_si,μ23_si,T,
                            "R_gsb", method,false; t2coh=t2coh) for i in 1:N) ./ N
    #rhots_NR_gsb = rhots_NR_gsb - conj(rhots_NR_gsb)
    #rhots_R_gsb  = rhots_R_gsb  - conj(rhots_R_gsb)
    if test
        corr_func_NR_gsb = corr_func_NR_gsb .* (exp.(-gτ_(-tlist)) * exp.(-gt_(-tlist))')
        corr_func_R_gsb  = corr_func_R_gsb  .* (exp.(-gτ_( tlist)) * exp.(-gt_( tlist))')
    end

    # stimulated emission
    corr_func_NR_se = sum(correlations(tlist,rho0_si,H_si[i],F_si,μ12_si,μ23_si,T,
                            "NR_se",method,false; t2coh=t2coh) for i in 1:N) ./ N
    corr_func_R_se  = sum(correlations(tlist,rho0_si,H_si[i],F_si,μ12_si,μ23_si,T,
                            "R_se", method,false; t2coh=t2coh) for i in 1:N) ./ N
    #rhots_NR_se = rhots_NR_se - conj(rhots_NR_se)
    #rhots_R_se  = rhots_R_se  - conj(rhots_R_se)
    if test
        corr_func_NR_se = corr_func_NR_se .* (exp.(-gτ_(-tlist)) * exp.(-gt_(-tlist))')
        corr_func_R_se  = corr_func_R_se  .* (exp.(-gτ_( tlist)) * exp.(-gt_( tlist))')
    end

    use_alt = true # use alternative way of calculating signal

    if use_alt
    
        corr_gsb = hcat(reverse(corr_func_R_gsb,dims=2), corr_func_NR_gsb[:,2:end])
        corr_se  = hcat(reverse(corr_func_R_se,dims=2),  corr_func_NR_se[:,2:end])
        corr_esa = hcat(reverse(corr_func_R_esa,dims=2), corr_func_NR_esa[:,2:end])
        corr_esax = hcat(reverse(corr_func_R_esax,dims=2),corr_func_NR_esax[:,2:end])

        corr = corr_gsb + corr_se - corr_esa - corr_esax

        spec2dd     = corr2spec(corr, zp)
        spec2dd_gsb = corr2spec(corr_gsb, zp)
        spec2dd_se  = corr2spec(corr_se, zp)
        spec2dd_esa = -corr2spec(corr_esa, zp)
        spec2dd_esax = -corr2spec(corr_esax, zp)

        #spec2d_NR_gsb = []
        #spec2d_R_gsb  = []
        #spec2d_NR_se  = []
        #spec2d_R_se   = []
        #spec2d_NR_esa = []
        #spec2d_R_esa  = []
        #spec2d_NR_esax = []
        #spec2d_R_esax  = []

    end
    
        # divide first value by .5 (see Hamm and Zanni)
        corr_func_NR_gsb[1,:] = corr_func_NR_gsb[1,:] / 2
        corr_func_NR_gsb[:,1] = corr_func_NR_gsb[:,1] / 2
        corr_func_R_gsb[1,:]  = corr_func_NR_gsb[1,:] / 2
        corr_func_R_gsb[:,1]  = corr_func_NR_gsb[:,1] / 2

        corr_func_NR_se[1,:] = corr_func_NR_se[1,:] / 2
        corr_func_NR_se[:,1] = corr_func_NR_se[:,1] / 2
        corr_func_R_se[1,:]  = corr_func_NR_se[1,:] / 2
        corr_func_R_se[:,1]  = corr_func_NR_se[:,1] / 2

        corr_func_NR_esa[1,:] = corr_func_NR_esa[1,:] / 2
        corr_func_NR_esa[:,1] = corr_func_NR_esa[:,1] / 2
        corr_func_R_esa[1,:]  = corr_func_NR_esa[1,:] / 2
        corr_func_R_esa[:,1]  = corr_func_NR_esa[:,1] / 2

        corr_func_NR_esax[1,:] = corr_func_NR_esax[1,:] / 2
        corr_func_NR_esax[:,1] = corr_func_NR_esax[:,1] / 2
        corr_func_R_esax[1,:]  = corr_func_NR_esax[1,:] / 2
        corr_func_R_esax[:,1]  = corr_func_NR_esax[:,1] / 2
        
        # zeropad data prior to fft to increase resolution
        corr_func_NR_gsb = zeropad(corr_func_NR_gsb,zp)
        corr_func_R_gsb  = zeropad(corr_func_R_gsb ,zp)
        
        corr_func_NR_se  = zeropad(corr_func_NR_se,zp)
        corr_func_R_se   = zeropad(corr_func_R_se ,zp)

        corr_func_NR_esa = zeropad(corr_func_NR_esa,zp)
        corr_func_R_esa  = zeropad(corr_func_R_esa ,zp)
        
        corr_func_NR_esax = zeropad(corr_func_NR_esax,zp)
        corr_func_R_esax  = zeropad(corr_func_R_esax ,zp)
        
        #######################################################
        ######### plot 2nd order correlation function #########
        #######################################################

        println("done\n")

        println("############# ########################### #############");
        println("############ Fourier transform of pathways ############");
        println("############# ########################### #############\n");

        spec2d_NR_gsb = fftshift(fft((corr_func_NR_gsb)))
        spec2d_R_gsb  = fftshift(fft((corr_func_R_gsb)))

        spec2d_NR_se  = fftshift(fft((corr_func_NR_se)))
        spec2d_R_se   = fftshift(fft((corr_func_R_se)))

        spec2d_NR_esa = -fftshift(fft((corr_func_NR_esa)))
        spec2d_R_esa  = -fftshift(fft((corr_func_R_esa)))

        spec2d_NR_esax = -fftshift(fft((corr_func_NR_esax)))
        spec2d_R_esax  = -fftshift(fft((corr_func_R_esax)))

        println("done\n")

        println("############# ########################### #############");
        println("########## Construct absorptive 2D spectrum ###########");
        println("############# ########################### #############\n");

        spec2d_R_gsb  = circshift(reverse(spec2d_R_gsb ,dims=1),(1,0))
        spec2d_R_se   = circshift(reverse(spec2d_R_se  ,dims=1),(1,0))
        spec2d_R_esa  = circshift(reverse(spec2d_R_esa ,dims=1),(1,0))
        spec2d_R_esax = circshift(reverse(spec2d_R_esax ,dims=1),(1,0))

    

    # calculate the absorptive spectrum of GSB, ESA, and SE
    spec2d_gsb  = spec2d_NR_gsb  + spec2d_R_gsb
    spec2d_se   = spec2d_NR_se   + spec2d_R_se
    spec2d_esa  = spec2d_NR_esa  + spec2d_R_esa
    spec2d_esax = spec2d_NR_esax + spec2d_R_esax

    # sum rephasing and non-rephasing spectra individually
    spec2d_r =  spec2d_R_gsb +
                spec2d_R_se  +
                spec2d_R_esa# +
                #spec2d_R_esax

    spec2d_nr = spec2d_NR_gsb +
                spec2d_NR_se  +
                spec2d_NR_esa #+
#                spec2d_NR_esax;

    # Subtract esax from gsb. Additional channel to simulate ground state recovery
    spec2d_gsb = spec2d_gsb + spec2d_esax

    # calculate the total 2D spectrum (complex)
    spec2d = spec2d_gsb +
             spec2d_se  +
             spec2d_esa

    println("done\n")


    #BUG #TODO: So far I need to change tlist and then crop ω when using alternative way of calculating 2D signals
    #tlist = [tlist; tlist[2:end] .+ tlist[end]]
    tlist, ω = interpt(tlist,zp)
    ω = ω[end÷2+1:end]

    # could return corr for looking at it ... set [] to save space
    #corr = []

    #FACT: spec2dd behaves correctly (for coupled dimer) when plotting absorptive, absolute, dispersive and phase // spec2d now as well ! ... 
    #TODO: absl is stretched too much in ω3!
    out = out2d{Array{ComplexF32,2}}(ω, spec2dd, spec2d_r, spec2d_nr, spec2d_gsb, spec2d_R_gsb,
                                     spec2d_NR_gsb, spec2d_se, spec2d_R_se, spec2d_NR_se, spec2d_esa,
                                     spec2d_R_esa, spec2d_NR_esa, corr)

    #out = crop2d(out,1;w_max=10,step=1) 
    out = round2d(out,2)

    return out
end

"""
    function crop2d()

Use to crop 2D spectra to smaller size.
"""
function crop2d(data_struct,w_min;w_max=100,step=1)

    min_i = argmin(abs.(data_struct.ω .- w_min))
    max_i = argmin(abs.(data_struct.ω .- w_max))

    data_struct.ω         = data_struct.ω[min_i:step:max_i]
    data_struct.full2d    = data_struct.full2d[min_i:step:max_i,min_i:step:max_i]
    data_struct.full2d_r  = data_struct.full2d_r[min_i:step:max_i,min_i:step:max_i]
    data_struct.full2d_nr = data_struct.full2d_nr[min_i:step:max_i,min_i:step:max_i]
    data_struct.gsb       = data_struct.gsb[min_i:step:max_i,min_i:step:max_i]
    data_struct.gsb_r     = data_struct.gsb_r[min_i:step:max_i,min_i:step:max_i]
    data_struct.gsb_nr    = data_struct.gsb_nr[min_i:step:max_i,min_i:step:max_i]
    data_struct.se        = data_struct.se[min_i:step:max_i,min_i:step:max_i]
    data_struct.se_r      = data_struct.se_r[min_i:step:max_i,min_i:step:max_i]
    data_struct.se_nr     = data_struct.se_nr[min_i:step:max_i,min_i:step:max_i]
    data_struct.esa       = data_struct.esa[min_i:step:max_i,min_i:step:max_i]
    data_struct.esa_r     = data_struct.esa_r[min_i:step:max_i,min_i:step:max_i]
    data_struct.esa_nr    = data_struct.esa_nr[min_i:step:max_i,min_i:step:max_i]

    return data_struct

end

"""
    function round2d()

Use to round 2D spectra to smaller size.
"""
function round2d(data_struct,digits=2)

    data_struct.ω         = round.(data_struct.ω,digits=6)
    data_struct.full2d    = round.(data_struct.full2d,digits=digits)
    data_struct.full2d_r  = round.(data_struct.full2d_r,digits=digits)
    data_struct.full2d_nr = round.(data_struct.full2d_nr,digits=digits)
    data_struct.gsb       = round.(data_struct.gsb,digits=digits)
    data_struct.gsb_r     = round.(data_struct.gsb_r,digits=digits)
    data_struct.gsb_nr    = round.(data_struct.gsb_nr,digits=digits)
    data_struct.se        = round.(data_struct.se,digits=digits)
    data_struct.se_r      = round.(data_struct.se_r,digits=digits)
    data_struct.se_nr     = round.(data_struct.se_nr,digits=digits)
    data_struct.esa       = round.(data_struct.esa,digits=digits)
    data_struct.esa_r     = round.(data_struct.esa_r,digits=digits)
    data_struct.esa_nr    = round.(data_struct.esa_nr,digits=digits)

    return data_struct

end

#CHECK I am considering R + conj(R) ... seems that 2D spectra turn out distorted
"""
    corr = correlations(tlist, rho0, H, F, μa, μb, pathway, method, debug; t2)

Calculate the response functions necessary for 2D spectroscopy.

# Arguments
* 'tlist'   : time steps at which to evaluate correlation functions τ and t
* 'rho0'    : ground state density matrix
* 'H'       : Hamiltonian of the system
* 'F'       : collapse operators (i.e. relaxation pathways during τ and t) OR
              relaxation tensor R for Redfield
* 'μa'      : transition dipole operator between ground- and excited states,
              responsible for GSB and SE
* 'μb'      : transition dipole operator between excited- and higher excited
              states, responsible for ESA
* 'pathway' : select response pathway 'R', or 'NR' _ 'gsb', 'se', or 'esa'
              'NR_gsb', 'R_esa', etc.
* 'method'  : select either Lindblad or Redfield master equation solver
* 'debug'   : if true, plot response functions and spectra of GSB, SE and ESA and
              show intermediate density matrices
* 't2'      : behaviour during t2 time "kin", "vib" or "full"

# Output
* 'corr'    : 2D correlation function

# Additional considerations:

* GS        : ground state
* ESs       : excited state or states
* FS        : final state
* μa and μb : TDMs, coupling states, such that GS <-μa-> ESs <-μb-> FS
* F         : is either list of collapse operators (c_ops; Lindblad) or
              a relaxation tensor (R; Redfield)


# Pathways calculated (w/o cc)
rephasing:
```
GSB               { μ * [(rho0 * μ)  * μ]}  + cc
            copy: ( μ * ((rho0 * μ)  * μ))  
SE                {[μ *  (rho0 * μ)] * μ}   + cc  
            copy: ((μ *  (rho0 * μ)) * μ)  
ESA               { μ *  [μ * (rho0  * μ)]} + cc  
            copy: ( μ *  (μ * (rho0  * μ)))  
```

non-rephasing:  
```
GSB               {μ * [μ * (μ * rho0)]} + cc 
            copy: (μ * (μ * (μ * rho0)))  
SE                {[(μ * rho0) * μ] * μ} + cc  
            copy: (((μ * rho0) * μ) * μ)  
ESA               {μ * [(μ * rho0) * μ]} + cc  
            copy: (μ * ((μ * rho0) * μ))  
```
"""
function correlations(tlist, rho0, H, F, μ_ge, μ_ef, T, pathway, method, debug; t2coh=false)

    #=
    #IDEA either create subspace in function that generates Hamiltonian, or create in t_ev
    try
        if pathway[end] == 'a' || pathway[end] == 'x'
            H, μa, μb, rho0, F =  create_subspace(H, μa, μb, rho0, F, "bi")
        elseif pathway[end] == 'b' || pathway[end] == 'e'
            H, μa, μb, rho0, F =  create_subspace(H, μa, μb, rho0, F, "si")
        else 
            println("error creating subspace")
        end
    catch
        println("error creating subspace")
    end
    =#


    ## use full transition dipole operator matrix
    #μ = μa + μb #DELETE
    ## or not ...
    #μ_ge = μa
    #μ_ef = μb

    #TODO!: T-dependent H not working with Redfield
    if length(H) == 1
        H   = H[1]
        H_T = H
    elseif length(H) == 2
        H_T = H[2]
        H   = H[1]
    else
        println("error T-dep H")
    end

    ## switch between Lindblad and Redfield master equation
    t_ev(A,B,C,D) =
        if method == "lindblad"
            #if isempty(methods(H))
                timeevolution.master(A,B,C,D)
            #elseif !isempty(methods(H))
            #    timeevolution.master_dynamic(A,B,C) #TODO: implement T-dep Hamiltonian, need to figure out how to create subspace when H is passed as a function
            #end
        elseif method == "redfield"
            timeevolution.master_bloch_redfield(A,B,D,C)
        else
            error("cmds -- No valid method selected")
        end

    ## add functionality to vary T
    T = [0.:1.:T;]     

    # initialize output matrix
    corr    = complex(zeros(length(tlist),length(tlist)))
    corr_cc = complex(zeros(length(tlist),length(tlist)))

    # use right-hand excitation for rephasing and left-hand excitation for
    # non-rephasing to ensure that emission is from the left of the double-
    # sided vertical Feynman diagram
    if pathway[1] == 'R'        # first interaction from right if rephasing
        #rho1    = rho0 * μ_ge
        rho1    = rho0 * tri(μ_ge,"U")
        #rho1_cc = μ_ge * rho0
        #rho1_cc = tri(μ_ge,"L") * rho0
    elseif pathway[1] == 'N'    # first interaction from left if non-rephasing
        #rho1    = μ_ge * rho0
        rho1    = tri(μ_ge,"L") * rho0
        #rho1_cc = rho0 * μ_ge
        #rho1_cc = rho0 * tri(μ_ge,"U")
    else
        error("cmds -- no pathway selected")
    end

    # define transition dipole operator to evaluate expectation value
    if pathway[end] == 'a' # if esa pathway use μb
        #μexp = -μb         # minus sign follows from Feynman diagrams
        μexp    = tri(μ_ef,"U") # ... removed -1
        #μexp_cc = -μ_ef
    elseif pathway[end] == 'e'
        μexp    = tri(μ_ge,"U")
        #μexp_cc = tri(μ_ge,"L")
    elseif pathway[end] == 'x'
        μ_ef    =  μ_ge
        μexp    =  μ_ge
        #μexp_cc = -μ_ge
        pathway = chop(pathway)
    else
        μexp    =  μ_ge
        #μexp_cc =  μ_ge
    end

    # calculate evolution during t1(τ)
    τ = tlist; tout, rho1_τ    = t_ev(τ,rho1,H,F)
    #τ = tlist; tout, rho1_τ_cc = t_ev(τ,rho1_cc,H,F)

    # add progress bar #TODO: make compatible with multiple Threads
    
    R = length(collect(1:20:length(rho1_τ)))
    Tpr = T[end]
    print(Threads.threadid())
    print(string(rpad(pathway,6), " at T = $Tpr fs. Progress: [" * repeat(" ", R-1)) * "]")
    print(repeat("\b", R))

    # TODO: damping function τ, implement somehow in redfield ??
    """
    τc = 8
    tc = 20
    gτ = (exp.(-τ./τc) .+ τ./τc .- 1)
    gt = (exp.(-τ./tc) .+ τ./tc .- 1)
    rho1_τ = rho1_τ .* exp.(-gτ)
    rho1_τ_cc = rho1_τ_cc .* exp.(-gt)
    """

    # now, use every rho(τ) as starting point to calculate second coherence
    # time (t) dependence. Note that first interaction with μ puts density
    # matrix back into a population without any time delay to the next
    # interaction (T=0)!
    for i=1:length(rho1_τ)

        if pathway == "R_gsb"
            #------------------------#
            # {μ * [(rho0 * μ) * μ]}
            # copy: (μ * ((rho0 * μ) * μ))
            #------------------------#
            #rho2    =   rho1_τ[i]     * μ_ge
            rho2    =   rho1_τ[i]    * tri(μ_ge,"L")
            #rho2_cc =   μ_ge         * rho1_τ_cc[i]
    #         rho2_cc =   tri(μ_ge,"U")* rho1_τ_cc[i]
            rho2a = rho2
            # eliminate on-diagonal or off-diagonal elements
            if t2coh == "kin"
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif t2coh == "vib"
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            elseif t2coh == "full"
                #
            end
            rho2b = rho2
            # time evolution during t2(T)-time
            if T[end] != 0
                tout, rho2_T = t_ev(T,rho2,H,F)
                rho2 = rho2_T[end]
    #            tout, rho2_T_cc = t_ev(T,rho2_cc,H,F)
    #            rho2_cc = rho2_T_cc[end]
            end
            # third field interaction
            #rho3   = μ_ge * rho2
            rho3   = tri(μ_ge,"L")    * rho2
            #rho3_cc   = μ_ge * rho2_cc
    #        rho3_cc   = tri(μ_ge,"U") * rho2_cc

        elseif pathway == "NR_gsb"
            #-------------------------#
            # {μ * [μ * (μ * rho0)]}
            # copy: (μ * (μ * (μ * rho0)))
            #-------------------------#
            #rho2 = μ_ge * rho1_τ[i]
            rho2 = tri(μ_ge,"U") * rho1_τ[i]
            #rho2_cc = rho1_τ_cc[i] * μ_ge
    #        rho2_cc = rho1_τ_cc[i] * tri(μ_ge,"L")

            if t2coh == "kin"
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif t2coh == "vib"
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            elseif t2coh == "full"
                #
            end

            if T[end] != 0
                tout, rho2_T = t_ev(T,rho2,H,F)
                rho2 = rho2_T[end]
    #            tout, rho2_T_cc = t_ev(T,rho2_cc,H,F)
    #            rho2_cc = rho2_T_cc[end]
            end

            #rho3   = μ_ge * rho2
            rho3   = tri(μ_ge,"L") * rho2
            #rho3_cc = rho2_cc * μ_ge
    #        rho3_cc = rho2_cc * tri(μ_ge,"U")

        elseif pathway == "R_se"
            #------------------------------#
            # {[μ * (rho0 * μ)] * μ} + cc
            # copy: ((μ * (rho0 * μ)) * μ)
            #------------------------------#
            #rho2 = μ_ge * rho1_τ[i]
            rho2 = tri(μ_ge,"L") * rho1_τ[i]
            #rho2_cc = rho1_τ_cc[i] * μ_ge
    #        rho2_cc = rho1_τ_cc[i] * tri(μ_ge,"U")

            if t2coh == "kin"
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif t2coh == "vib"
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            elseif t2coh == "full"
                #
            end

            if T[end] != 0
                tout, rho2_T = t_ev(T,rho2,H,F)
                rho2 = rho2_T[end]
    #            tout, rho2_T_cc = t_ev(T,rho2_cc,H,F)
    #            rho2_cc = rho2_T_cc[end]
            end

            #rho3   =       rho2 * μ_ge
            rho3   =       rho2 * tri(μ_ge,"L")
            #rho3_cc = μ_ge * rho2
    #        rho3_cc = tri(μ_ge,"L") * rho2

        elseif pathway == "NR_se"
            #------------------------------#
            # {[(μ * rho0) * μ] * μ} + cc
            # copy: (((μ * rho0) * μ) * μ)
            #------------------------------#
            #rho2 =      rho1_τ[i]   * μ_ge
            rho2 =      rho1_τ[i]   * tri(μ_ge,"U")
            #rho2_cc = μ_ge * rho1_τ_cc[i]
    #        rho2_cc = tri(μ_ge,"L") * rho1_τ_cc[i]

            if t2coh == "kin"
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif t2coh == "vib"
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            elseif t2coh == "full"
                #
            end

            if T[end] != 0
                tout, rho2_T = t_ev(T,rho2,H,F)
                rho2 = rho2_T[end]
    #            tout, rho2_T_cc = t_ev(T,rho2_cc,H,F)
    #            rho2_cc = rho2_T_cc[end]
            end

            #rho2.data[1,1] = 0

            #rho3   =       rho2 * μ_ge
            rho3   =       rho2 * tri(μ_ge,"L")
            #rho3_cc = μ_ge * rho2
    #        rho3_cc = tri(μ_ge,"U") * rho2

        elseif pathway == "R_esa"
            #μ_ef = μ_ge# + μ_ef # after relaxation, reabsorption is part of esa (!?)
            #------------------------------#
            # {μ * [μ * (rho0 * μ)]} + cc
            # copy: (μ * (μ * (rho0 * μ)))
            #------------------------------#
            #rho2 = μ_ge * rho1_τ[i]
            rho2 = tri(μ_ge,"L")    * rho1_τ[i]
            #rho2_cc = rho1_τ_cc[i]  * μ_ge
    #        rho2_cc = rho1_τ_cc[i]  * tri(μ_ge,"U")

            if t2coh == "kin"
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif t2coh == "vib"
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            elseif t2coh == "full"
            end

            if T[end] != 0
                tout, rho2_T = t_ev(T,rho2,H,F)
                rho2 = rho2_T[end]
    #            tout, rho2_T_cc = t_ev(T,rho2_cc,H,F)
    #            rho2_cc = rho2_T_cc[end]
            end

            #rho3   = μ_ef * rho2
            rho3    = tri(μ_ef,"L") * rho2
            #rho3_cc = rho2_cc       * μ_ef
    #        rho3_cc = rho2_cc       * tri(μ_ef,"U")

        elseif pathway == "NR_esa"
            #μ_ef = μ_ge #+ μ_ef
            #-----------------------------#
            # {μ * [(μ * rho0) * μ]} + cc
            # copy: (μ * ((μ * rho0) * μ))
            #-----------------------------#
            #rho2 =      rho1_τ[i]   * μ_ge
            rho2 =      rho1_τ[i]   * tri(μ_ge,"U")
            #rho2 =      μ_ge * rho1_τ[i]

            #rho2_cc = μ_ge * rho1_τ_cc[i]
    #        rho2_cc = tri(μ_ge,"L") * rho1_τ_cc[i]
            #rho2_cc = rho1_τ_cc[i] * μ_ge

            if t2coh == "kin"
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif t2coh == "vib"
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            elseif t2coh == "full"
            end

            if T[end] != 0
                tout, rho2_T = t_ev(T,rho2,H,F)
                rho2 = rho2_T[end]
    #            tout, rho2_T_cc = t_ev(T,rho2_cc,H,F)
    #            rho2_cc = rho2_T_cc[end]
            end

            #rho3    = μ_ef * rho2
            rho3    = tri(μ_ef,"L") * rho2
            #rho3_cc = rho2_cc * μ_ef
    #        rho3_cc = rho2_cc * tri(μ_ef,"U")
            
        end

        # calculate time evolution of ρ during detection time t
        t = tlist; tout, rho3_t     = t_ev(t,rho3,H_T,F)

    #    t = tlist; tout, rho3_t_cc  = t_ev(t,rho3_cc,H,F)

        # calc. exp. v. for GSB and SE (μexp = μ12) and for ESA (μexp = μ23)
        corr[i,:]    = (expect(μexp,conj(rho3_t))) # .* exp.(-gt) #TODO: lineshape f
    #    corr_cc[i,:] = real(expect(μexp_cc,conj(rho3_t_cc)))

        #DELETE
        #corr = corr + corr_cc
        #TODO: damping function t
        #D = .4
        #τc = 10
        #gt = (exp.(-t./τc) .+ t./τc .- 1)
        #corr[i,:] = corr[i,:] .* exp.(-gt)

        # ALTERNATIVELY take the trace: Tr(Aρ)
        #for j=1:length(t)
        #    corr_cc[i, j] = tr(real(rhot_temp_cc[j].data * μexp.data))
        #end

        # progress bar update
        if mod(i,20) == 0
            print("=")
        end

        # print a few combinations to check in debug mode
        if debug
            if i in [1]
                println("\n")
                println("###################################################")
                println("###################### DEBUG ######################")
                println("###################################################\n")
                println("showing density matrices from during the evolution.\n")
                println("Which path ?")
                println(pathway)
                println(" ")
                println("μ_ge")
                println(dense(μ_ge))
                println(" ")
                println("μ_ef")
                println(dense(μ_ef))
                println(" ")
                println("rho0")
                println(rho0)
                println(" ")
                println("rho1")
                println(rho1)
                println(" ")
                if T[end] != 0
                    println(" ")
                    println("rho2_T[1] (before T evolution)")
                    println(rho2_T[1])
                end
                println("rho2a")
                println(rho2a)
                println("rho2b")
                println(rho2b)
                println("rho2")
                println(rho2)
                println(" ")
                println("rho3")
                println(rho3)
                println(" ")
                println("μexp")
                println(μexp)
                println(" ")
                print("μexp * rho3")
                println(μexp * rho3)
                println(" ")
                println("\n")
            end
        end
    end
    print("\n")
    return corr
end

"""
    ω, corr, spec

Quick function to calculate absorption spectrum.
"""
function absorptionSpectrum(t, rho, H, F, μ; zp=11)

    corr = timecorrelations.correlation(t, rho, H, F, μ, μ)
    corr_intp = zeropad(corr,zp)
    t_intp, ω_intp = interpt(t,zp)
    ω, spec = timecorrelations.correlation2spectrum(t, corr; normalize_spec=true)
    return ω, corr, spec
end

#FACT: tri is required for properly making ... σ⁺ ρ ... (tf does this mean?) 
#FACT: deactivate tri for displaced harmonic oscillator
"""
    out = tri(dat, uplo)

function used internally by cmds.jl to select upper triagular (uplo="U") or lower
triangular (uplo="L") matrix.
"""
function tri(dat,uplo)
    out = copy(dat)
    if uplo == "L"
        out.data = tril(out.data)
    elseif uplo == "U"
        out.data = triu(out.data)
    end
    return out;
end

"""
    plot2d_comps(data)

Plots GSB, SE and ESA component of 2D spectrum.

#Arguments
* data  : data struct for single time step T (e.g. out2d[1])

"""
function plot2d_comps(data)

    fig, ax = subplots(1,4,sharex=true,sharey=true,figsize=(15,3))

    sca(ax[1])
    plot2d(data.ω, data.full2d)
    PyPlot.title("absorptive")

    sca(ax[2])
    plot2d(data.ω, data.gsb)
    PyPlot.title("GSB")

    sca(ax[3])
    plot2d(data.ω, data.se)
    PyPlot.title("SE")

    sca(ax[4])
    plot2d(data.ω, data.esa)
    PyPlot.title("ESA")

    tight_layout()
end

#TODO: is this useful ? 
function pretty_show_mat(data)
    fig, ax = plt.subplots()
    im = ax.imshow(data)
    for i in 1:size(data,1)
        for j in 1:size(data,2)
            ax.text(i-1,j-1,round(data[i,j],digits=2),ha="center",va="center")
        end
    end
end

"""
    function vib_analysis(data)

Used to quickly plot GSB, SE and ESA rephasing (R) and non-rephasing (NR) spectra. Useful for checking
for beating patterns, but could also be used to look at individual components of "normal" 2D dataset.

# Arguments
* 'data'   : just the 2D as (e.g. out2d[1])
"""
function vib_analysis(data)
     fig, ax = subplots(3,2,sharex=true,sharey=true,figsize=(2*3.2,3*3))
     sca(ax[1,1])
     plot2d(data.ω,abs.(data.gsb_r))
     title("GSB R")
     sca(ax[1,2])
     plot2d(data.ω,abs.(data.gsb_nr))
     title("GSB NR")
     sca(ax[2,1])
     plot2d(data.ω,abs.(data.se_r))
     title("SE R")
     sca(ax[2,2])
     plot2d(data.ω,abs.(data.se_nr))
     title("SE NR")
     sca(ax[3,1])
     plot2d(data.ω,abs.(data.esa_r))
     title("ESA R")
     sca(ax[3,2])
     plot2d(data.ω,abs.(data.esa_nr))
     title("ESA NR")
     tight_layout()
end

"""
    plot_levels(H,center_x;col,ls)

Plots the energy levels of H at position center_x within an axis. Create a figure() first, 
to plot the energy levels of different systems, subspaces, monomers, etc. together. E.g.

```julia
figure(); title("Energy level diagram")
cmds.plot_levels(H,0)
cmds.plot_levels(Hcav,-1)
cmds.plot_levels(Hexc,1)
cmds.plot_levels(Hmon,2)
xticks([-1, 0, 1, 2], ["cavity", "full", "matter", "monomer"])
```

"""
function plot_levels(H, center_x; col="k", ls="solid")
    # get energies and eigenvectors of full system
    E      = eigvals(dense(H).data)
    vecs   = eigvecs(dense(H).data)
    eivecs = [Ket(H.basis_l,vecs[:,i]) for i in 1:length(E)]
    # plot eigenenergies
    xn = 0
    lx = center_x - .45
    rx = center_x + .45
    for i in 1:length(E)
        if i > 1 #&& E[i] - E[i-1] < 0.1
            hlines(E[i],lx,rx,col,ls)
            xn += .02
        else
            hlines(E[i],lx,rx,col,ls)
            xn = 0
        end
    end

end

function logger(what, fn)
    open(fn,"a") do io
        # do stuff with the open file
        println(io,what)
    end
end

function logger(what,dat, fn)
    open(fn,"a") do io
        # do stuff with the open file
        println(io,what,dat)
    end
end

"""
    corr2spec(corr, zp)

Auxiliary function to convert correlation map to 2D spectrum.
"""
function corr2spec(corr, zp)
        corr = vcat(corr[end:-1:1,end:-1:1],corr[2:end,:])
        corr = vcat(zeros(size(corr,2))',corr)
        corr = hcat(corr,zeros(size(corr,1)))

        N = 2^zp
        temp = complex(zeros(N,N));

        temp[(N-size(corr,1))÷2+1:end-(N-size(corr,1))÷2,(N-size(corr,1))÷2+1:end-(N-size(corr,1))÷2] = corr
        corr = temp

        corr = circshift(corr,(0,1))

        spec = fft(fftshift(corr))
        spec = fftshift(spec)
        #spec = circshift(spec,(-100,-100))
        spec = spec[1:end÷2,1:end÷2]
        spec = reverse(spec,dims=1)
        spec = reverse(spec,dims=2)
        return spec
    end

end # module end
