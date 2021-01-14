# CMDS

Simulate coherent multidimensional spectroscopy signals from quantum mechanical models.

## Introduction

The code relies primarily on [qojulia/QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) and a tutorial can be found [here](https://docs.qojulia.org/). Further examples and functionalities are found the Python project [QuTiP](http://qutip.org/).

The module [CMDS.jl](/cmds.jl) contains the necessary functions to calculate 2D spectra and will be described below. [examples/](/examples) shows example scenarios.

## Installation

CMDS requires the Julia language and [qojulia/QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl), which can be installed from via the standard sources:

[Julia](https://docs.julialang.org/en/v1/manual/getting-started/)

[QoJulia](https://docs.qojulia.org/installation/)

## CMDS.jl - Functions

Type ``?cmds.<function>`` into the REPL to access the documentation for a certain function.

### Available functions:

- create_colormap

- zeropad

- interpt

- make2Dspectra

- correlations

- view_dm_evo

- save_2d

- plot2d

- crop2d

### How to use:

initialize output array

`` out2d = Array{cmds.out2d}(undef, length(T)) ``

where T is a vector containing population time steps.

Next call `` cmds.make2dspectra `` in a for loop

```julia
for i = 1:length(T)
    out2d[i] = cmds.make2Dspectra(tlist,rho0,H,F,μ12,μ23,T[i],"lindblad";debug=false,zp=zp)
end
```

with tlist, rho0, H, F, ...

Using multithreading, several population time steps can be evaluated simultaneously (make sure to disable all output plots within cmds.make2Dspectra(), as these might crash the execution):
```julia
Threads.@threads for i = 1:length(T)
    out2d[i] = cmds.make2Dspectra(...)
end
```
## Examples

The following [/examples](/examples) are available.

### coupled_dimer.jl

The properties (angles, coupling strength, etc.) of a coupled dimer system are calculated and QuantumOptics.jl is used to calculate the correlation function and linear absorption spectrum.

![coupledDimer](/example_images/coupledDimer.png)

CMDS.jl uses QuantumOptics.jl to calculate the response functions in a four-wave mixing experiment and calculates the expected 2D spectrum.

![coupledDimer 2D spectrum](/example_images/coupledDimer2D.png)

### displaced_harmonic_oscillator_model.jl

[text]

![displacedHarmonicOscillator](/example_images/displHarmOsc.png)

![displacedHarmonicOscillator 2D spectrum](/example_images/displHarmOsc2D.png)

### FCF_morse-potential.jl

[text]

![FCF Morse Potential](/example_images/FCfactorsMorsePot.png)

### Jaynes-Cummings model

The coupling between a quantized optical field and a two-level system is described by the Jaynes-Cummings Hamiltonian

H = ωᵣ a† a + ωₐ σ₊ σ₋ + ( a† σ₋ + a σ₊ )

```julia
H = wc * at * a + wa * sp * sm + g * (at * sm + a * sp)
```

blabla

![Jaynes-Cummings](/example_images/JaynesCummingsSpectrum.png)

2D spectrum of the Jaynes-Cummings model at different delays of the population time T.

![Jaynes-Cummings 2D spectrum](/example_images/JaynesCummingsSpectrum2D.png)


### Ensemble of two-level systems with disorder

...

### Tavis-Cummings

In order to go beyond the Jaynes-Cummings model ...
