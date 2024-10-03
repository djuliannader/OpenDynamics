module main
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
import wigner
import build
import opendynamics


Nfock=22
N=400
L=22
xav=2.0
pav=2.0
Delta=0.0
epsilon=0.0
K=0.0
tm = 8.0
timep=3.0
HH = build.Hamiltonian(Nfock,Delta,epsilon,K)
rho0 = build.initialrho(Nfock,xav,pav)
#rho0 = build.initialrhomix(Nfock,xav,pav)
a = build.creation(Nfock)
ad = transpose(conj(a))
jumpop = [a,ad]
jumppar = [0.1,0.1]
sp = opendynamics.survivalp(Nfock,HH,rho0,jumppar,jumpop,tm)

println(sp)

wopen = opendynamics.wigneropen_t(Nfock,HH,rho0,timep,L,N,jumppar,jumpop,"Wopentest.dat")

println(wopen)

end


