module opendynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
import wigner
import build


function survivalp(Nmax,Ham,rho0,jpar,jop,tmax)
  times=(0.0,tmax)
  tint=0.05
  f(u,p,t) = -im*(Ham*u-u*Ham) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) )
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=tint)
 #println("--here --")
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
 println("-------------   Go to file survivalprobability.dat to see the survival probability  ----------------")
  return "done"
end


function wigneropen_t(Nmax,Ham,rho0,time,L,N,jpar,jop,string)
  times=(0.0,time[length(time)])
  tint=0.05
  f(u,p,t) = -im*(Ham*u-u*Ham) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) ) 
  prob = ODEProblem(f,rho0,times)
  sol = solve(prob,Tsit5(),alg_hints = [:stiff],dt=tint)
  wlist=[]
  wnlist=[]
  for k in 1:length(time)
  rhot = sol(time[k])
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
  ww=wigner.wigner_mix(w,wstates,L,N,string[k])
  append!(wlist,ww[1])
  append!(wnlist,ww[2])
  end
  return [wlist,wnlist]
end

function expectation2modes(op,Ham,rho0,jpar,jop,tmax,tmsp)
  times=(0.0,tmax)
  f(u,p,t) = -im*(Ham*u-u*Ham) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) ) 
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob)
 tint=0.05
 nt = floor(Int, tmax/tint)
 surprob = []
 open("expectation.dat","w") do io
 for i in 1:(nt+1)
   rhot = sol(tint*(i-1))
   expinst = tr(rhot*op)
   println(io,tint*(i-1)," ", round(real(expinst),digits=16))
 end
 end
   rhot = sol(tmsp)
   expinst = tr(rhot*op)
   println("Nb at the time period: ",expinst)
   rhot = sol(tmax)
   expinst = tr(rhot*op)
   println("Nb asymptotic value: ",expinst)
  return "done"
end



end


