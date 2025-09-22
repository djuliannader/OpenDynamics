module test_oscar
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("build.jl")
include("opendynamics.jl")
using .build
using .opendynamics


Nfock=50         # Size of the Fock space
Delta = -1.0
K = 0.2
epi = 1.5 
epf = 0.4        # 
tm = 10.0        # Maximum time for survival probability
tint = 0.1       # time intervals for survival probability


#  ------ Caulculating Open dynamics----------------
nintervals = floor(Int, tm/tint)
timep=[i*tint for i in 1:nintervals]
HH0 = build.HamiltonianKerr(Nfock,Delta,epi,K)
Hvecs=eigvecs(HH0)
psi0 = [Hvecs[i,1] for i in 1:(Nfock+1)]
HH = build.HamiltonianKerr(Nfock,Delta,epf,K)
sp = opendynamics.survivalamplitude(Nfock,HH,psi0,tm,tint)
# --------------------------------------





# ------  Printing results ----------------
# println("Size of the Fock space: ",Nfock)
#println("Survivaval amplitude: ",sp)
# -----------------------------------------


end


