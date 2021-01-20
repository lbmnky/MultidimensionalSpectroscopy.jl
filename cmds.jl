module cmds

export create_colormap, zeropad, interpt, make2Dspectra, correlations,
        view_dm_evo, save_2d, plot2d, crop2d, tri

using QuantumOptics, FFTW, LinearAlgebra, PyPlot, Colors, DelimitedFiles, Printf



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
    if repr == "absorptive"
        data = real(data)
    elseif repr == "dispersive"
        data = imag(data)
    elseif repr == "absolute"
        data = abs.(data)
    end

    ## normalize 2D spectrum to plot
    if norm_spec
        data = data/maximum(abs.(data));
    end

    ## make levels for contourplot with different scaling
    lvls = [-1.025:0.05:1.025;]
    if scaling == "lin"
        #
    elseif scaling == "tan"
        lvls = tan.(lvls); lvls = lvls ./ maximum(lvls);
    elseif scaling == "cube"
        lvls = lvls.^3;
    elseif scaling == "asinh"
        lvls = asinh.(3*lvls); lvls = lvls ./ maximum(lvls);
    end

    ## rescale lvls to data (don't rescale data, takes more time)
    m = -1*(maximum(abs.(real(data)))); M = (maximum(abs.(real(data))))
    if m != 0 && M != 0
        lvls = lvls * maximum(abs.([m, M]));
        lvls_ticks = [-1:0.2:1;] * maximum(abs.([m, M]))
    else
        lvls = 0.001 * lvls;
    end

    # normalize to global maximum with evolution time T scan
    if norm == 0
    else
        lvls = lvls ./ maximum(lvls) * norm;
    end

    ## plot data as filled contour
    conv = "ω3ω1"
    if conv == "ω1ω3"
        contourf(ω,ω,data,lvls,cmap=create_colormap());
        xlabel("ω₃"); ylabel("ω₁"); #clim([m, M]);
    elseif conv == "ω3ω1"
        contourf(ω,ω,transpose(data),lvls,cmap=create_colormap());
        xlabel("ω₁"); ylabel("ω₃");
    end
    plot([ω[1], ω[end]], [ω[1], ω[end]],"k--");

    ## plot positive and negative contourlines with different linestyles
    #contour(ω,ω,temp.*(temp.>0),lvls[1:1:end],cmap="Greys",
    #        linewidths=0.5,linestyles="solid");
    #contour(ω,ω,temp.*(temp.<0),lvls[1:1:end],cmap="Greys",
    #        linewidths=0.5,linestyles="dashed");

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
    writedlm(fn_base * "T_steps.dat", T, ',')
    for i = 1:length(T)
        fn = fn_base * @sprintf("2Dspec_%03i.dat",i)
        writedlm(fn, spec2d[i], ',')
        println("\nFile saved as: " * fn)
    end
    println("\nFiles saved!")
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
        title("tlist # $j")
    end
    tight_layout()
end

"""
    cmp = create_colormap(scheme)

Create a "bright" or "dark" custom colormap to display 2D spectra.
"""
function create_colormap(scheme="bright")

    ## define colors from most negative (MMM) ... to ... most positive (PPP)
    MMM = RGB(0.718, 0, 0.718);
    MM = RGB(0.516, 0.516, 0.991);
    M = RGB(0, 0.559, 0.559);
    P = RGB(0, 0.592, 0);
    PP = RGB(0.827, 0.827, 0);
    PPP = RGB(0.947, 0.057, 0.057);

    ## zero value is either "bright" or "dark"
    if scheme == "bright"
        Z = RGB(.93,.93,.93);
    elseif scheme == "dark"
        Z = RGB(.03,.03,.03);
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
        return dat
    end

    # Sanity check that N is not too large
    if N > 2^12
        error("Large N detected. Check if N isn't too large.")
    end

    if N > 0
        # switch between 1D and 2D array
        if size(dat,2) == 1
            out = complex(zeros(2^N))
            println("\nZeropadding 1D-data to $(2^N) values\n")
        elseif size(dat,2) > 1
            out = complex(zeros(2^N,2^N))
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
        tlist = tlist
    end

    # Sanity check that N is not too large
    if N > 2^12
        error("Large N detected. Check if N isn't too large.")
    end

    if N > 0
        # extrapolate tlist / get new tlist[end]
        dt = abs(tlist[end]-tlist[1]) / (length(tlist)-1)
        tlist = [0:(2^N)-1;] * dt .+ tlist[1]
    end

    # create ω axis corresponding to zeropadded data
    dω = 1/tlist[end]
    ω = ([1:length(tlist);] .- length(tlist)/2 ) * dω * 2 * π
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
mutable struct out2d
    ω
    full2d
    full2d_r
    full2d_nr
    gsb
    gsb_r
    gsb_nr
    se
    se_r
    se_nr
    esa
    esa_r
    esa_nr
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
"""
function make2Dspectra(tlist, rho0, H, F, μ12, μ23, T, method; debug=false, use_sub=true, zp=0)

    println("############# ########################### #############");
    println("############# calculate R and NR pathways #############");
    println("############# ########################### #############\n");

    # TODO make subspace for GSB and SE pathways with only single excited manifold -- DOne

    Eivecs = eigvecs(dense(H).data)

    L = length(H.basis_l.bases)
    n_sub = 1 + L + binomial(L,2)
    Eivecs_sub = [Ket(H.basis_l,Eivecs[:,i]) for i in 1:n_sub]
    b_sub = SubspaceBasis(H.basis_l,Eivecs_sub)

    P = sparse(projector(b_sub, (H.basis_l)))
    Pt = dagger(P)

    if use_sub
        println("dim(subspace): ", length(b_sub))
        H_bi    = P * H * Pt
        μ12_bi  = P * μ12 * Pt
        μ23_bi  = P * μ23 * Pt
        rho0_bi = P * rho0 * Pt
        F_bi    = [P * F[i] * Pt for i in 1:length(F)]
    else
        H_bi    = H
        μ12_bi  = μ12
        μ23_bi  = μ23
        rho0_bi = rho0
        F_bi    = F
    end

    # excited state absorption
    corr_func_NR_esa = correlations(tlist,rho0_bi,H_bi,F_bi,μ12_bi,μ23_bi,T,"NR_esa",method,false)
    corr_func_R_esa  = correlations(tlist,rho0_bi,H_bi,F_bi,μ12_bi,μ23_bi,T,"R_esa", method,false)
    #rhots_NR_esa = rhots_NR_esa - conj(rhots_NR_esa)
    #rhots_R_esa  = rhots_R_esa  - conj(rhots_R_esa)

    # excited state absorption ... from GS
    corr_func_NR_esax = correlations(tlist,rho0_bi,H_bi,F_bi,μ12_bi,μ23_bi,T,"NR_esax",method,false)
    corr_func_R_esax  = correlations(tlist,rho0_bi,H_bi,F_bi,μ12_bi,μ23_bi,T,"R_esax", method,false)

    n_sub       = 1 + L
    Eivecs_sub  = [Ket(H.basis_l,Eivecs[:,i]) for i in 1:n_sub]
    b_sub       = SubspaceBasis(H.basis_l,Eivecs_sub)

    P  = sparse(projector(b_sub, (H.basis_l)))
    Pt = dagger(P)

    if use_sub
        println("dim(subspace): ", length(b_sub))
        H_si    = P * H * Pt
        μ12_si  = P * μ12 * Pt
        μ23_si  = P * μ23 * Pt
        rho0_si = P * rho0 * Pt
        F_si    = [P * F[i] * Pt for i in 1:length(F)]
    else
        H_si    = H
        μ12_si  = μ12
        μ23_si  = μ23
        rho0_si = rho0
        F_si    = F
    end

    # ground state bleach
    corr_func_NR_gsb = correlations(tlist,rho0_si,H_si,F_si,μ12_si,μ23_si,T,"NR_gsb",method,false)
    corr_func_R_gsb  = correlations(tlist,rho0_si,H_si,F_si,μ12_si,μ23_si,T,"R_gsb", method,true)
    #rhots_NR_gsb = rhots_NR_gsb - conj(rhots_NR_gsb)
    #rhots_R_gsb  = rhots_R_gsb  - conj(rhots_R_gsb)

    # stimulated emission
    corr_func_NR_se = correlations(tlist,rho0_si,H_si,F_si,μ12_si,μ23_si,T,"NR_se",method,false)
    corr_func_R_se  = correlations(tlist,rho0_si,H_si,F_si,μ12_si,μ23_si,T,"R_se", method,false)
    #rhots_NR_se = rhots_NR_se - conj(rhots_NR_se)
    #rhots_R_se  = rhots_R_se  - conj(rhots_R_se)


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

    spec2d_NR_esa = fftshift(fft((corr_func_NR_esa)))
    spec2d_R_esa  = fftshift(fft((corr_func_R_esa)))

    spec2d_NR_esax = fftshift(fft((corr_func_NR_esax)))
    spec2d_R_esax  = fftshift(fft((corr_func_R_esax)))

    println("done\n")

    println("############# ########################### #############");
    println("########## Construct absorptive 2D spectrum ###########");
    println("############# ########################### #############\n");

    spec2d_R_gsb = circshift(reverse(spec2d_R_gsb ,dims=2),(0,1))
    spec2d_R_se  = circshift(reverse(spec2d_R_se ,dims=2),(0,1))
    spec2d_R_esa = circshift(reverse(spec2d_R_esa ,dims=2),(0,1))
    spec2d_R_esax = circshift(reverse(spec2d_R_esax ,dims=2),(0,1))

    # calculate the absorptive spectrum of GSB, ESA, and SE
    spec2d_gsb = spec2d_NR_gsb + spec2d_R_gsb
    spec2d_se  = spec2d_NR_se  + spec2d_R_se
    spec2d_esa = spec2d_NR_esa + spec2d_R_esa
    spec2d_esax = spec2d_NR_esax + spec2d_R_esax

    # sum rephasing and non-rephasing spectra individually
    spec2d_r = circshift(reverse(spec2d_R_gsb ,dims=2),(0,1)) +
               circshift(reverse(spec2d_R_se  ,dims=2),(0,1)) +
               circshift(reverse(spec2d_R_esa ,dims=2),(0,1))# +
#               circshift(reverse(spec2d_R_esax,dims=2),(0,1));

    spec2d_nr = spec2d_NR_gsb +
                spec2d_NR_se  +
                spec2d_NR_esa #+
#                spec2d_NR_esax;

    # calculate the total 2D spectrum (complex)
    spec2d = spec2d_gsb +
             spec2d_se  +
             spec2d_esa +
             spec2d_esax;

    println("done\n")

    # cut spec2d
    tlist, ω = interpt(tlist,zp)

    out = out2d(ω, spec2d, spec2d_r, spec2d_nr, spec2d_gsb, spec2d_R_gsb,
                spec2d_NR_gsb, spec2d_se, spec2d_R_se, spec2d_NR_se, spec2d_esa,
                spec2d_R_esa, spec2d_NR_esa)

    if debug == true
        #=
        figure(figsize=(8,6))
        subplot(431); pcolormesh(real(rhots_R_gsb)); title("GSB"); ylabel("Rephasing")
        subplot(434); pcolormesh(real(rhots_NR_gsb)); ylabel("Non-rephasing")
        subplot(437); pcolormesh(spec2d_R_gsb); ylabel("Rephasing")
        subplot(4,3,10); pcolormesh(spec2d_NR_gsb); ylabel("Non-rephasing")

        subplot(432); pcolormesh(real(rhots_R_se)); title("SE")
        subplot(435); pcolormesh(real(rhots_NR_se));
        subplot(438); pcolormesh(spec2d_R_se);
        subplot(4,3,11); pcolormesh(spec2d_NR_se);

        subplot(433); pcolormesh(real(rhots_R_esa)); title("ESA")
        subplot(436); pcolormesh(real(rhots_NR_esa));
        subplot(439); pcolormesh(spec2d_R_esa);
        subplot(4,3,12); pcolormesh(spec2d_NR_esa);
        tight_layout()
        =#
        """
        figure(figsize=(10,6))
        L = size(spec2d,1);
        subplot(221); pcolormesh(real(spec2d_gsb)); title("GSB abs."); colorbar()
        clim([-100, 100])
        plot([1:L;],[1:L;])
        subplot(222); pcolormesh(real(spec2d_se)); title("SE abs."); colorbar()
        clim([-100, 100])
        plot([1:L;],[1:L;])
        subplot(223); pcolormesh(real(spec2d_esa)); title("ESA abs."); colorbar()
        clim([-100, 100])
        tight_layout()
        plot([1:L;],[1:L;])
        subplot(224); pcolormesh(real(spec2d_esax)); title("ESAx abs."); colorbar()
        clim([-100, 100])
        tight_layout()
        plot([1:L;],[1:L;])
        """
        return out

    else

        return out #crop2d(out,0)

    end
end

"""
need help
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

### I am considering R + conj(R) ... seems that 2D spectra turn out distorted
"""
    corr = correlations(tlist, rho0, H, F, μa, μb, pathway, method, debug)

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

# Output
* 'corr'    : 2D correlation function

# Additional considerations:

* GS        : ground state
* ESs       : excited state or states
* FS        : final state
* μa and μb : TDMs, coupling states, such that GS <-μa-> ESs <-μb-> FS
* F         : is either list of collapse operators (c_ops; Lindblad) or
            of relaxation tensor (R, Redfield)
"""
function correlations(tlist, rho0, H, F, μa, μb, T, pathway, method, debug)

    ## use full transition dipole operator matrix
    μ = μa + μb
    ## or not ...
    μ_ge = μa
    μ_ef = μb

    ## switch between Lindblad and Redfield master equation
    t_ev(A,B,C,D) =
        if method == "lindblad"
            timeevolution.master(A,B,C,D)
        elseif method == "redfield"
            timeevolution.master_bloch_redfield(A,B,D,C)
        else
            error("cmds -- No valid method selected")
        end

    ## add functionality to vary T
    print(string(pathway, " at T = $T fs ..."))
    T = [0.:1.:T;]     # cAN I START FROM -1.0 to include T = 0 ?

    # SELECTING ONLY DIAGONAL ELEMENTS AT THIS POINT ELIMINATES PATHWAYS
    # THAT OSCILLATE DURING t2
    XX = false   # get rid of off-diagonal elements during t2
    YY = true  # get rid of on-diagonal elements during t2

    # initialize output matrix
    corr    = zeros(length(tlist),length(tlist))
    corr_cc = zeros(length(tlist),length(tlist))


    # use right-hand excitation for rephasing and left-hand excitation for
    # non-rephasing to ensure that emission is from the left of the double-
    # sided vertical Feynman diagram
    if pathway[1] == 'R'        # first interaction from right if rephasing
        #rho1    = rho0 * μ_ge
        rho1    = rho0 * tri(μ_ge,"U")
        #rho1_cc = μ_ge * rho0
        rho1_cc = tri(μ_ge,"L") * rho0
    elseif pathway[1] == 'N'    # first interaction from left if non-rephasing
        #rho1    = μ_ge * rho0
        rho1    = tri(μ_ge,"L") * rho0
        #rho1_cc = rho0 * μ_ge
        rho1_cc = rho0 * tri(μ_ge,"U")
    else
        error("cmds -- no pathway selected")
    end

    # define transition dipole operator to evaluate expectation value
    if pathway[end] == 'a' # if esa pathway use μb
        #μexp = -μb         # minus sign follows from Feynman diagrams
        μexp    = -tri(μ_ef,"U")
        μexp_cc = -μ_ef
    elseif pathway[end] == 'e'
        μexp    = tri(μ_ge,"U")
        μexp_cc = tri(μ_ge,"L")
    elseif pathway[end] == 'x'
        μ_ef    =  μ_ge
        μexp    = -μ_ge
        μexp_cc = -μ_ge
        pathway = chop(pathway)
    else
        μexp    =  μ_ge
        μexp_cc =  μ_ge
    end

    # calculate evolution during t1(τ)
    τ = tlist; tout, rho1_τ    = t_ev(τ,rho1,H,F)
    τ = tlist; tout, rho1_τ_cc = t_ev(τ,rho1_cc,H,F)

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
            # second field interaction A(τ)
            #rho2    =   rho1_τ[i]     * μ_ge
            rho2    =   rho1_τ[i]    * tri(μ_ge,"L")
            #rho2_cc =   μ_ge         * rho1_τ_cc[i]
    #         rho2_cc =   tri(μ_ge,"U")* rho1_τ_cc[i]
            rho2a = rho2
            # eliminate on-diagonal or off-diagonal elements
            if XX
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif YY
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            else
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

            if XX
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif YY
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            else
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

            if XX
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif YY
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            else
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

            if XX
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif YY
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            else
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

            if XX
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif YY
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            else
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

            if XX
                rho2.data = tril(triu(rho2.data))
    #            rho2_cc.data = tril(triu(rho2_cc.data))
            elseif YY
                rho2.data = tril(rho2.data,-1) + triu(rho2.data,1)
    #            rho2_cc.data = tril(rho2_cc.data,-1) + triu(rho2_cc.data,1)
            else
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
        t = tlist; tout, rho3_t     = t_ev(t,rho3,H,F)

    #    t = tlist; tout, rho3_t_cc  = t_ev(t,rho3_cc,H,F)



        # calc. exp. v. for GSB and SE (μexp = μ12) and for ESA (μexp = μ23)
        corr[i,:]    = real(expect(μexp,conj(rho3_t)))
    #    corr_cc[i,:] = real(expect(μexp_cc,conj(rho3_t_cc)))

        #corr = corr + corr_cc

        # ALTERNATIVELY take the trace: Tr(Aρ)
        #for j=1:length(t)
        #    corr_cc[i, j] = tr(real(rhot_temp_cc[j].data * μexp.data))
        #end

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
    print(" done!\n")
    return corr
end

"""
    out = tri(dat, uplo)

function used internally by cmds.jl to select upper triagular (uplo="U") or lower
triangular (uplo="L") matrix.
"""
function tri(dat,uplo)
    out = copy(dat)
    #if uplo == "L"
    #    out.data = tril(out.data)
    #elseif uplo == "U"
    #    out.data = triu(out.data)
    #end
    return out;
end

end # module
