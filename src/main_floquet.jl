module main
push!(LOAD_PATH, pwd())
using LinearAlgebra
include("modules/troterization.jl")
include("modules/build.jl")
using .troterization
using .build


Nfock=50         # Size of the Fock space
N=400            # Finite differences
L=20             # Size of the phase space
Delta=-2.0        # Parameter of the Kerr Hamiltonian
epsilon=8.0      # Parameter of the Kerr Hamiltonian
K=1.0            # Parameter of the Kerr Hamiltonian
nu = 27.3135     # frequency of the drive 
xi = 0.25         # strength of the drive
Nf = 1000         # Number of intervals for the troterization
kk = 1            # state of interest
name = "wigner_floquet0.dat"  # name of the data for wigner funtion


#  ------ Caulculating Floquet state----------------
tt = troterization.floquetstate(K,epsilon,Delta,Nfock,nu,xi,Nf,L,N,name,kk)
# --------------------------------------





# ------  Printing results ----------------
println("Results for the ",kk," Floquet state")
println("Wigner volume of the ",kk," floquet state : ",round(tt[1][1], digits=8))
println("Negativity of the ",kk," floquet state    : ",round(tt[1][2], digits=8))
println("Expectation value <Ek|H0|Ek> for k= ",kk,"  eigenstate       : ",round(tt[2], digits=8))
println("Expectation value <Fk|H0|Fk> for k= ",kk,"  floquet state    : ",round(real(tt[3]), digits=8))
# -----------------------------------------


# ------  Printing results ----------------
open("resonances3.dat","w") do io
nulist = [(120.0+0.1*i) for i in 1:300]
for i in 1:length(nulist)
tt = troterization.floquetstate(K,epsilon,Delta,Nfock,nulist[i],xi,Nf,L,N,name,kk)
println(io,nulist[i],"   ",round(tt[1][2], digits=8))
println(nulist[i],"   ",round(tt[1][2], digits=8))
end
end
# -----------------------------------------





#HH = build.HamiltonianKerr(Nfock,Delta,epsilon,K)
#egg = eigvals(HH)
#println("flag: ",HH)


end


