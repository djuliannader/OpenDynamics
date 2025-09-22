module test_oscar
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("build.jl")
include("opendynamics.jl")
include("wigner.jl")
using .build
using .opendynamics
using .wigner


Nfock= 50         # Size of the Fock space
N    = 400
L    = 20      
Delta = -1.0
ai = 2
af = -5
b  = 1 
tm = 0.5        # Maximum time for survival probability
tint = 0.1       # time intervals for survival probability


#  ------ Caulculating Open dynamics----------------
nintervals = floor(Int, tm/tint)
timep=[i*tint for i in 1:nintervals]
HH0 = build.HamiltonianQuartic(Nfock,ai,b)
Hvecs=eigvecs(HH0)
psi0 = [Hvecs[i,1] for i in 1:(Nfock+1)]
HH = build.HamiltonianQuartic(Nfock,af,b)
psit = exp(-im*tm*HH)*psi0
wf1=wigner.focktowf(psit,L,N)
ww=wigner.wignerf(wf1,L,N,"../output/test_wignerquartic.dat")
# --------------------------------------





# ------  Printing results ----------------
# println("Size of the Fock space: ",Nfock)
#println("Survivaval amplitude: ",sp)
# -----------------------------------------


end


