using MultidimensionalSpectroscopy, QuantumOptics, LinearAlgebra, FFTW, Colors,
    Printf, DelimitedFiles, Plots

# make sure to set script directory as pwd()
cd(@__DIR__)

# to use PyPlot in GUI mode under Linux
#pygui(true)

# run once with "calc_2d = false" to initialize functions
calc_2d = false

cmp = create_colormap("bright");
##

## define functions
function draw_dipole(x, y, α, d)
    quiver(x - d / 2 * cos(α), y - d / 2 * sin(α), d * cos(α), d * sin(α), angles="xy",
        scale_units="xy", scale=1, color="r")
    scatter(x, y)
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
r, β = 1.155, angle2rad(23.147495)
#r , β          = 19.3155, angle2rad(0)
# position of dipole 2 (calculated from above)
x₂, y₂ = x₁ + r * cos(β), y₁ + r * sin(β)
# angle and moment of dipole 2 (input)
α₂, d₂ = angle2rad(0), 0.25

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
#figure(figsize = (10, 4));
#ax1 = subplot(231);
#grid("on")



f = scatter([x₁, x₂], [y₁, y₂], legend=false)#, linestyle = "--", color = "k")

f = plot(f,f,layout=(2,2))
display(f)

