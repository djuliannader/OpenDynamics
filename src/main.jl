module main
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
include("modules/wigner.jl")
include("modules/build.jl")
include("modules/opendynamics.jl")
using .wigner
using .build
using .opendynamics



Nfock=50         # Size of the Fock space
N=400            # Finite differences
L=30             # Size of the phase space
xav=-2.0         # Average position of the initial coherent state
pav=2.0          # Average momentum of the initial coherent state
Delta=-1.0       # Parameter of the Kerr Hamiltonian
epsilon0=1.5     # Quench parameter of the Kerr Hamiltonian  (initial)
epsilonf=0.4     # Quench parameter of the Kerr Hamiltonian (final)
K=0.5            # Parameter of the Kerr Hamiltonian
nis = 1           # initial state (NFock-1)
tm = 7            # Maximum time for survival probability
tint = 0.01       # Time interval for the SP and LE
nshots = 2        # Number of shoots for the Wigner function
tint2 = 0.05      # Time interval for the Wigner function shoots
dper = 1.0        # period of the drive
dpar = 0.0        # amplitude of the drive
acc = 1e-16       # Accuraccy goal for the differential equations


jumppar = [0.0,0.0]      # Jump parameters


#  ------ Caulculating Open dynamics----------------
#timep=[(i-1)*tint2 for i in 1:nshots]
timep=[0,4]
outputlist=["output/wignerfunction"*string(i-1)*"_out.dat" for i in 1:nshots]
#timep=[]
#for i in 1:(nshots)
# append!(timep,(i-1)*tint2+0.5)
#end
#outputlist=["output/wignerfunction_out.dat" for i in 1:nshots]
HH = build.HamiltonianKerr(Nfock,Delta,epsilonf,K)
#diagH0  = [0.0+0*im for i in 1:(Nfock+1)]
#HH = Diagonal(diagH0)
#rho0 = build.initialrhofock(Nfock,xav,pav,nis)
#rho0 = build.initialrhocat(Nfock,al)
rho0 = build.initialrhoquench(Nfock,epsilon0,Delta,K)
a    = build.anhilation(Nfock)
ad   = transpose(a)
dop   = ad*a
nop   = ad*a
jumpop = [a,ad]
Nexp = tr(rho0*nop)
sp = opendynamics.survivalp(Nfock,HH,rho0,jumppar,jumpop,dpar,dop,dper,tm,tint,acc)
wopen = opendynamics.wigneropen_t(Nfock,HH,rho0,timep,L,N,K,jumppar,jumpop,dpar,dop,dper,acc,outputlist)
# --------------------------------------


# ------  Printing results ----------------
 println("Size of the Fock space: ",Nfock)
 println("Expectation value <n> of the initial state: ",real(Nexp))
 println("Number of subintervals for the finite difference: ",N)
println("-------------   Go to file output/wignerresults.dat to see the following   ----------------")
println("-------------        - 1st column: time                             ----------------")
println("-------------        - 2nd column: Volume of the Wigner function    ----------------")
println("-------------        - 3rd column: Wigner negativities              ----------------")
println("-- Go to files output/_out.dat to see the shoots of the dynamics of the Wigner function --")
# -----------------------------------------


end


