# CMDS

v0.2

Simulate coherent multidimensional spectroscopy signals using from quantum mechanical models.

## Introduction

The code relies primarily on [qojulia/QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl). The module [CMDS.jl](https://github.com/lbmnky/CMDS/blob/master/cmds.jl) contains all the necessary functions. [examples/](https://github.com/lbmnky/CMDS/tree/master/examples) shows example scenarios.

## Installation

CMDS requires the Julia language and [qojulia/QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl), which can be installed from via the standard sources:

Julia: [Installation](https://docs.julialang.org/en/v1/manual/getting-started/)

QuantumOptics.jl:
[Installation](https://docs.qojulia.org/installation/)


## examples/

The following examples are available.

### coupled_dimer.jl

The properties (angles, coupling strength, etc.) of a coupled dimer system are calculated and QuantumOptics.jl is used to calculate the correlation function and linear absorption spectrum.

![coupledDimer](/example_images/coupledDimer.png)

CMDS.jl uses QuantumOptics.jl to calculate the response functions in a four-wave mixing experiment and calculates the expected 2D spectrum.

![coupledDimer 2D spectrum](/example_images/coupledDimer2D.png)


### displaced_harmonic_oscillator_model.jl

[text]

### FCF_morse-potential.jl

[text]

### tavis-cummings

[text]
