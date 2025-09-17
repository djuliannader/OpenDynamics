module troterization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export floquetstate
include("build.jl")
include("wigner.jl")
using .build
using .wigner


function troter(H0,VV,b,nu,nn,Nmax)
 pi=acos(-1)
 T=2*pi/nu
# b=nu*b
 dt=T/nn
#--- building the time independent Hamiltonian ----
 # diagonal matrix elements
 H1=H0
 U0=exp(-im*H1*dt/2)
#--- builiding the time dependent Hamiltonian
 ddiag=[1 for i in 1:(Nmax+1)]
 Ut=b*VV
 U = Diagonal(ddiag)
#---- Troterization ----------------- 
 for i in 1:nn
  facu1=(im/nu)*(-sin(i*nu*dt)+sin((i-1)*nu*dt))
  U1=exp(facu1*Ut)
  UM=U0*U1*U0
  U=UM*U
 end
 #println(H1)
 #println(U)
 return U
end


function floquetstate(K,ep,delta,Nmax,nu,xi,Nf,L,N,name,kk)
  H0   = build.HamiltonianKerr(Nmax,delta,ep,K)
  aop  = build.creation(Nmax)
  adop = transpose(aop)
  #VV   = (adop+aop)/2^(1/2)
  VV   = adop*aop
  floquet = troter(H0,VV,xi,nu,Nf,Nmax)
  evs = eigvecs(floquet)
  evls = eigvals(H0)
  floquetord =  Floquetstatesordering(floquet,H0,Nmax)
  listvec=[evs[i,floquetord[kk]] for i in 1:(Nmax+1)]
  wf = wigner.focktowf(listvec,L,N)
  wwig = wigner.wignerf(wf,L,N,name)
  expvalf = transpose(listvec)*H0*listvec
  return [wwig,evls[kk],expvalf]
end

function Floquetstatesordering(U,H0,Nmax)
   eigvecsf=eigvecs(U)
   evfvec = zeros(0)
   for i in 1:(Nmax+1)
     eigvecf=[eigvecsf[j,i] for j in 1:(Nmax+1)]
     eigvecf_t=Array{Complex{Float64}}(undef,1,length(eigvecf)) 
     for k in 1:length(eigvecf)
       eigvecf_t[1,k]=conj(eigvecf[k])
     end
     evf=eigvecf_t*H0*eigvecf
     append!(evfvec, real(evf) )
   end
   evford=sortperm(evfvec)
   #println(evford)   #  descomentar para ver como se ordenan las quasienergias
   return evford
   end
   
end


