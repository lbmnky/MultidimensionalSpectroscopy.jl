# Examples "systems coupled to cavities"

## Jaynes-Cummings model

The coupling between a quantized optical field and a two-level system is described by the Jaynes-Cummings Hamiltonian

H = ω<sub>r</subR-NR> a<sup>†</sup> a + ω<sub>a</sub> σ<sub>+</sub> σ<sub>-</sub> + ( a<sup>†</sup> σ<sub>-</sub> + a σ<sub>+</sub> )

or for you to copy:

```julia
H = wc * at * a + wa * sp * sm + g * (at * sm + a * sp)
```

Here, ω<sub>r</sub> is the energy/frequency/... of the cavity mode, a<sup>†</sup>(a) is the  ... and σ<sub>+</sub>(σ<sub>-</sub>)  the ... . The calculated linear absorption spectrum of the system looks pointy:

![Jaynes-Cummings](ex_images/JaynesCummingsSpectrum.png)

2D spectrum of the Jaynes-Cummings model at different delays of the population time T.

![Jaynes-Cummings 2D spectrum](ex_images/JaynesCummingsSpectrum2D.png)
