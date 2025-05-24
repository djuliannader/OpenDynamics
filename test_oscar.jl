module test_oscar
push!(LOAD_PATH, pwd())
using LinearAlgebra
import build
import opendynamics


Nfock=50         # Size of the Fock space
xav=0.0          # Average position of the initial coherent state
pav=0.0          # Average momentum of the initial coherent state
w0=1.0           # Parameter of the Oscar
E0=0.0           # Parameter of the Oscar
xi=0.0           # Parameter of the Oscar
xf=2.0           # Parameter of the Oscar
tm = 40.0        # Maximum time for survival probability
tint = 0.1       # time intervals for survival probability


#  ------ Caulculating Open dynamics----------------
nintervals = floor(Int, tm/tint)
timep=[i*tint for i in 1:nintervals]
HH0 = build.HamiltonianOsc(Nfock,w0,xi,E0)
Hvecs=eigvecs(HH0)
psi0 = [Hvecs[i,1] for i in 1:(Nfock+1)]
HH = build.HamiltonianOsc(Nfock,w0,xf,E0)
sp = opendynamics.survivalamplitude(Nfock,HH,psi0,tm,tint)
# --------------------------------------





# ------  Printing results ----------------
 println("Size of the Fock space: ",Nfock)
println("Survivaval amplitude: ",sp)
# -----------------------------------------


end


