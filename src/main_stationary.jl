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


Nfock=30         # Size of the Fock space
N=500            # Finite differences
L=20             # Size of the phase space
xav=-2.0       # Average position of the initial coherent state
pav=2.0          # Average momentum of the initial coherent state
Delta=-2.0       # Parameter of the Kerr Hamiltonian
epsilon=0.0      # Parameter of the Kerr Hamiltonian
K=1.0            # Parameter of the Kerr Hamiltonian
tm = 0.1         # Maximum time for survival probability
t0 = 0.0         # First shot of Wigner Fucntions
nshots=1        # Number of shoots for the Wigner function
tint = 0.1       # Time interval for the Wigner function shoots
kk=1   		 # Fock state to be displaced
jumppar = [0.05,0.05]      # Jump parameters


#----------Calculating GS Wigner function--------------
timep = collect(t0:tint:t0 + (nshots-1)*tint)
   a = build.anhilation(Nfock)
   ad = transpose(conj(a))
   jumpop = [a,ad]
   jumppar = [0.0,0.0]
for i in 3.0:0.1:3.3
   outputlist=["output/gsDelta=$(i).dat" for j in 1:nshots]
   #output="gsDelta=$(i)"
   HH=build.HamiltonianKerr(Nfock,i,epsilon,K)
   rho0 = build.initialrhoGS(HH)
   wopen = opendynamics.wigneropen_t(Nfock,HH,rho0,timep,L,N,jumppar,jumpop,outputlist)
end

# ------  Printing results ----------------
 println("Size of the Fock space: ",Nfock)
 println("Number of subintervals for the finite difference: ",N)
println("-------------   Go to file output/wignerresults.dat to see the following   ----------------")
println("-------------        - 1st column: time                             ----------------")
println("-------------        - 2nd column: Volume of the Wigner function    ----------------")
println("-------------        - 3rd column: Wigner negativities              ----------------")
println("-- Go to files output/_out.dat to see the shoots of the dynamics of the Wigner function --")
# -----------------------------------------


end


