module main
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
include("wigner.jl")
include("build.jl")
using .wigner
using .build



Nfock=50         # Size of the Fock space
a=10              # Hamiltonian par
theta=0.1

#  ------ Caulculating Open dynamics----------------
aop    = build.anhilation(Nfock)
adop   = transpose(conj(aop))
xop    = (1/2^(1/2))*(aop+adop)
pop    = im*(1/2^(1/2))*(adop-aop)
HH=(1/2)*pop^2+a*xop^2-xop^4
HHs=(1/2)*exp(-2*im*theta)*pop^2+a*exp(2*im*theta)*xop^2-exp(2*im*theta)*xop^4
ev=eigvals(HH)
evs=eigvals(HHs)
println("eigvals  : ",ev)
println("eigvals_s: ",evs)
# --------------------------------------



end


