#!/usr/bin/julia


rho0 = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]
mua = [0 1 1 0; 1 0 0 0; 1 0 0 0; 0 0 0 0]
mub = [0 0 0 0; 0 0 0 1; 0 0 0 1; 0 1 1 0]

println("rho0")
display(rho0)
print("mu g-e")
display(mua)
print("mu e-f")
display(mub)

println("GSB")
println("R")
println("step0")
display(rho0)
println("step1")
display(rho0 * mua)
println("step2")
display((rho0 * mua) * mua)
println("step3")
display(mua * ((rho0 * mua) * mua))
