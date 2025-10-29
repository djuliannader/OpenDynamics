module opendynamics
push!(LOAD_PATH, pwd())
using LinearAlgebra
using DifferentialEquations
include("wigner.jl")
include("fisher.jl")
include("build.jl")
using .wigner
using .build
using .fisher


function survivalp(Nmax,Ham,rho0,jpar,jop,dop,dpar,dper,Aop,tmax,tint,acc)
  times=(0.0,tmax)
  f(u,p,t) = -im*((Ham+dpar*dop*cos(2*pi*t/dper))*u-u*(Ham+dpar*dop*cos(2*pi*t/dper))) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) )
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff],dt=tint)
    #println("--here --")
 
 nt = floor(Int, tmax/tint)
 surprob = []
 open("output/survivalprobability.dat","w") do io
 open("output/loschmidt_echo.dat","w") do io2
 open("output/fisher.dat","w") do io3
 for i in 1:(nt+1)
   rhot = sol(tint*(i-1))
   QFI = fisher.FisherInf(rhot, Aop)
   sp = (tr((rhot^(1/2)*rho0*rhot^(1/2))^(1/2)))^2
   fidinst = (tr((rhot^(1/2)*rho0*rhot^(1/2))^(1/2)))
   mat = rho0*rhot
   ls = [mat[i,i] for i in 1:(Nmax+1)]
   lecho   = sum(ls)
   println(io,tint*(i-1)," ", round(real(sp),digits=8))
   println(io2,tint*(i-1)," ", round(real(lecho),digits=8)," ",round(imag(lecho),digits=8)," ",round(real(fidinst),digits=8)," ",round(imag(fidinst),digits=8))
   println(io3,tint*(i-1)," ", round(real(QFI),digits=8))
end
 end
 end
 end
 println("-------------   Go to file output/survivalprobability.dat to see the survival probability  ----------------")
 println("-------------   Go to file output/loschmidt_echo.dat to see the Loschmidt echo amplitude   ----------------")
 println("-------------   Go to file output/fisher.dat to see the Fisher Information   ----------------")
  return "done"
end





function wigneropen_t(Nmax,Ham,rho0,time,L,N,K,jpar,jop,dpar,dop,dper,acc,string)
  times=(0.0,time[length(time)])
  tint=0.01
  f(u,p,t) = -im*((Ham+dpar*dop*cos(2*pi*t/dper))*u-u*(Ham+dpar*dop*cos(2*pi*t/dper))) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) ) 
  prob = ODEProblem(f,rho0,times)
  sol = solve(prob,abstol = acc,Tsit5(),alg_hints = [:stiff],dt=tint)
  wlist=[]
  wnlist=[]
  open("output/wignerresults.dat","w") do io
  for k in 1:length(time)
  println("Calculating Wigner function at time:",time[k])
  rhot = sol(time[k])
  #rhot = exp(-im*Ham*time[k])*rho0*exp(im*Ham*time[k])
  evals=eigvals(rhot)
  evecs=eigvecs(rhot)
  eps = 0.01
  iflag = length(evals)
  w=[]
  wstates=[]
  while real(evals[iflag])>eps
    append!(w,real(evals[iflag]))
    println("component: ",real(evals[iflag]))
    comp = [evecs[i,iflag] for i in 1:(Nmax+1)]
    compwf = wigner.focktowf(comp,L,N)
    append!(wstates,[compwf])
    iflag=iflag-1
  end
  ww=wigner.wigner_mix(w,wstates,L,N,K,string[k])
  println(io,time[k]," ",trunc(ww[1],digits=6)," ", trunc(ww[2],digits=6))
  append!(wlist,ww[1])
  append!(wnlist,ww[2])
  end
  end
  return "done"
end



function expectation2modes(op,Ham,rho0,jpar,jop,tmax,tmsp)
  times=(0.0,tmax)
  f(u,p,t) = -im*(Ham*u-u*Ham) + jpar[1]*(jop[1]*u*transpose(conj(jop[1])) - (1/2)*(transpose(conj(jop[1]))*jop[1]*u + u*transpose(conj(jop[1]))*jop[1]) ) +  jpar[2]*(jop[2]*u*transpose(conj(jop[2])) - (1/2)*(transpose(conj(jop[2]))*jop[2]*u + u*transpose(conj(jop[2]))*jop[2]) ) 
 prob = ODEProblem(f,rho0,times)
 sol = solve(prob)
 tint=0.05
 nt = floor(Int, tmax/tint)
 surprob = []
 open("output/expectation.dat","w") do io
 for i in 1:(nt+1)
   rhot = sol(tint*(i-1))
   expinst = tr(rhot*op)
   println(io,tint*(i-1)," ", round(real(expinst),sigdigits=8))
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

function survivalamplitude(N,Ham,psi0,tm,tint)
    psi0a = transpose(conj(psi0))
    time = [i for i in 0:tint:tm]
    open("../output/survivalamplitudetest.dat","w") do io
    for ti in time
        psit = exp(-im*ti*Ham)*psi0
        sa = psi0a*psit
        println(io,ti," ", round(real(sa[1]),digits=8)," ",round(imag(sa[1]),digits=8))
   end
   end
end

end
