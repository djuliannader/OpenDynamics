module opendynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
import wigner
import build


function survivalp(Nmax,Ham,rho0,jpar,jop,tmax)
  times=(0.0,tmax)
  f(u,p,t) = -im*(Ham*u-u*Ham) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) ) 
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob)
 tint=0.05
 nt = floor(Int, tmax/tint)
 surprob = []
 open("survivalprobability.dat","w") do io
 for i in 1:(nt+1)
   rhot = sol(tint*(i-1))
   fidinst = (tr((rhot^(1/2)*rho0*rhot^(1/2))^(1/2)))^2
   println(io,tint*(i-1)," ", round(real(fidinst),digits=16))
 end
 end
  return "done"
end


function wigneropen_t(Nmax,Ham,rho0,time,L,N,jpar,jop,string)
  times=(0.0,time)
  f(u,p,t) = -im*(Ham*u-u*Ham) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) ) 
  prob = ODEProblem(f,rho0,times)
  sol = solve(prob)
  rhot = sol(time)
  evals=eigvals(rhot)
  evecs=eigvecs(rhot)
  eps = 0.001
  iflag = length(evals)
  w=[]
  wstates=[]
  while real(evals[iflag])>eps
    append!(w,real(evals[iflag]))
    comp = [evecs[i,iflag] for i in 1:(Nmax+1)]
    compwf = wigner.focktowf(comp,L,N)
    append!(wstates,[compwf])
    iflag=iflag-1
  end
  ww=wigner.wigner_mix(w,wstates,L,N,string)
  return "done"
end




end


