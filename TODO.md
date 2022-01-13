# To do !

- Comments, comments, comments ...
- generic way to construct Hamiltonians *
- understand many-body formalism

- symbolic matrix algebra to follow the "flow"

- ESA twice as strong as GSB/SE ? They need to cancel out, whenever there is no shift between them  ... (true ?) #FACT: True depending on model, neglect Lindblad terms between 2nd and 1st exc sector, since their presence seem to lead to broadening of ESA peaks in w3

- can convolute 2D signal with laser by #CHECK Does this only affect spectrum ? #TODO: how to convolute in time ? 
```julia
E = exp.(-(out2d[1].ω.-4.25).^2/(2*(1/2.355)^2))
```
 and
 ```julia
2Dconv = out2d[1].full2d .* (E.^2*(E.^2)')
```

- Add TA generation from 2D ... done by integrating over full 2d spectrum, #TODO: filter by excitation





How to construct full Hamiltonian
http://hitoshi.berkeley.edu/221A/tensorproduct.pdf
