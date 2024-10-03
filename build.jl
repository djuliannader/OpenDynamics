module build
push!(LOAD_PATH, pwd())
using LinearAlgebra


function Hamiltonian(Nmax,delta,epsilon,K)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab)
  id = Diagonal(diagid)
  adop = transpose(aop)
  Ham = adop*aop+(1/2)*id
  #Ham = -delta*adop*aop-epsilon*(adop^2+aop^2)-K*adop^2*aop^2 
  return Ham
end

function creation(Nmax)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab) 
  return aop
end

function coherentstate(Nmax,xc,pc)
  al=(1/2^(1/2))*(xc+im*pc)
  cs = [ exp(-abs2(al)/2)*al^n/sqrt(factorial(big(n))) for n in 0:Nmax]
  return cs
end

function initialrho(Nmax,xc,pc)
   cs = coherentstate(Nmax,xc,pc)
   cstp = transpose(conj(cs))
   rhoin = cs*cstp
   return rhoin
end

function initialrhomix(Nmax,xc,pc)
   cs1 = coherentstate(Nmax,xc,pc)
   cstp1 = transpose(conj(cs1))
   cs2 = coherentstate(Nmax,-xc,-pc)
   cstp2 = transpose(conj(cs2))
   rhoin1 = cs1*cstp1
   rhoin2 = cs2*cstp2
   return (1/2)*rhoin1 + (1/2)*rhoin2
end



end


