module main
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
import wigner
import build
import opendynamics


Nfock=20         # Size of the Fock space
N=400            # Finite differences
L=20             # Size of the phase space
xav=-2.0         # Average position of the initial coherent state
pav=2.0          # Average momentum of the initial coherent state
Delta=1.0        # Parameter of the Kerr Hamiltonian
epsilon=0.0      # Parameter of the Kerr Hamiltonian
K=1.0            # Parameter of the Kerr Hamiltonian
tm = 1.0         # Maximum time for survival probability
nshots=6         # Number of shoots for the Wigner function
tint = 0.1       # Time interval for the Wigner function shoots


#  ------ Caulculating ----------------
timep=[(i-1)*tint for i in 1:nshots]
outputlist=["wignerfunction"*string(i-1)*"_out.dat" for i in 1:nshots]
HH = build.HamiltonianKerr(Nfock,Delta,epsilon,K)
rho0 = build.initialrho(Nfock,xav,pav)
#rho0 = build.initialrhomix(Nfock,xav,pav)
a = build.creation(Nfock)
ad = transpose(conj(a))
jumpop = [a,ad]
jumppar = [0.05,0.05]
sp = opendynamics.survivalp(Nfock,HH,rho0,jumppar,jumpop,tm)
wopen = opendynamics.wigneropen_t(Nfock,HH,rho0,timep,L,N,jumppar,jumpop,outputlist)
# --------------------------------------





# ------  Printing results ----------------
println("Survivaval probability: ",sp)
println("-------------   Go to file wignerresults.dat to see the following   ----------------")
println("-------------        - 1st column: time                             ----------------")
println("-------------        - 2nd column: Volume of the Wigner function    ----------------")
println("-------------        - 3rd column: Wigner negativities              ----------------")
open("wignerresults.dat","w") do io
  for i in 1:length(timep)
    println(io,timep[i]," ",trunc(wopen[1][i],digits=6)," ",trunc(wopen[2][i],digits=6))
  end
end
println("-- Go to files _out.dat to see the shoots of the dynamics of the Wigner function --")
# -----------------------------------------


end


