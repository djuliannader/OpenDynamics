module main
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
import wigner
import build
import opendynamics


Nfock=20                # Size of the Fock space
N=400                   # Finite differences
L=20                    # Size of the phase space
xav=-2.0                # Average position of the initial coherent state
pav=2.0                 # Average momentum of the initial coherent state
Delta=-2.0               # Parameter of the Kerr Hamiltonian
epsilon=0.0             # Parameter of the Kerr Hamiltonian
K=1.0                   # Parameter of the Kerr Hamiltonian
tm = 10.0               # Maximum time for survival probability
nshots=9                # Number of shoots for the Wigner function
tint = 1.0              # Time interval for the Wigner function shoots
jumppar = [0.05,0.0]     # Jump parameters


#  ------ Caulculating Open dynamics----------------
#timep=[(i-1)*tint for i in 1:nshots]
timep=[0 pi/2 pi 3*pi/2 2*pi 5*pi/2 3*pi 7*pi/2 4*pi]
outputlist=["output/wignerfunction"*string(i-1)*"_out.dat" for i in 1:nshots]
#outputlist=["wignerfunction0_out.dat" for i in 1:nshots]
HH = build.HamiltonianKerr(Nfock,Delta,epsilon,K)
rho0 = build.initialrho(Nfock,xav,pav)
psi0 = build.initialpsi(Nfock,xav,pav)
#rho0 = build.initialrhomix(Nfock,xav,pav)
a    = build.anhilation(Nfock)
ad   = transpose(conj(a))
jumpop = [a,ad]
sp = opendynamics.survivalp(Nfock,HH,rho0,jumppar,jumpop,tm)
#sa = opendynamics.survivalamplitude(Nfock,HH,psi0,tm,0.01)
wopen = opendynamics.wigneropen_t(Nfock,HH,rho0,timep,L,N,jumppar,jumpop,outputlist)
# --------------------------------------





# ------  Printing results ----------------
 println("Size of the Fock space: ",Nfock)
 println("Number of subintervals for the finite difference: ",N)
println("-------------   Go to file output/wignerresults.dat to see the following   ----------------")
println("-------------        - 1st column: time                             ----------------")
println("-------------        - 2nd column: Volume of the Wigner function    ----------------")
println("-------------        - 3rd column: Wigner negativities              ----------------")
open("output/wignerresults.dat","w") do io
  for i in 1:length(timep)
    println(io,timep[i]," ",trunc(wopen[1][i],digits=6)," ",trunc(wopen[2][i],digits=6))
    println("Wigner function at time: ", timep[i], " obtained")
  end
end
println("-- Go to files output/_out.dat to see the shoots of the dynamics of the Wigner function --")
# -----------------------------------------


end


