module build
push!(LOAD_PATH, pwd())
using LinearAlgebra


function HamiltonianKerr(Nmax,delta,ep,K)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab)
  id = Diagonal(diagid)
  adop = transpose(aop)
  xop=(1/2^(1/2))*(aop+adop)
  pop=(-im)*(1/2^(1/2))*(aop-adop)
  apar=10.0
  #Ham = adop*aop+(1/2)*id
  #Ham = -delta*adop*aop - ep*(adop^2+aop^2) + K*adop^2*aop^2
  Ham = pop^2/2 - apar*xop^2 + xop^4
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
  csn = []
  for i in 1:length(cs)
  recs = convert(Float64, real(cs[i]))
  imcs = convert(Float64, imag(cs[i]))
  append!(csn,recs+im*imcs)
  end
  return csn
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

function Hamiltonian2modes(Nmaxa,Nmaxb,lam)
  ADdiag = [0.0 for i in 1:(Nmaxa+1)*(Nmaxb+1)]
  BDdiag = [0.0 for i in 1:(Nmaxa+1)*(Nmaxb+1)]
  AD = diagm(ADdiag)
  BD = diagm(BDdiag)
  for i in 1:Nmaxa
     for j in 1:(Nmaxb+1)
        AD[i+(Nmaxa+1)*(j-1)+1,i+(Nmaxa+1)*(j-1)] = i^(1/2)
     end
  end
  for i in 1:(Nmaxa+1)
     for j in 1:Nmaxb
        BD[i+(Nmaxa+1)*(j),i+(Nmaxa+1)*(j-1)] = j^(1/2)
     end
  end
  A=transpose((AD))
  B=transpose((BD))
  Ham=lam*(A*BD + B*AD)
  return [Ham,A,B]
end

function rho0coherent2modes(Nmaxa,Nmaxb,xav,pav)
   al = (1/2^(1/2))*(xav + im*pav)
   rho0diag = [0.0+0.0*im for i in 1:(Nmaxa+1)*(Nmaxb+1)]
   rho0 = diagm(rho0diag)
   for k in 1:(Nmaxa+1)
      rho0[k,k]=abs(exp(-abs(al)^2)*(al^(2*(k-1)))/factorial(k-1))
   end
   return rho0
end

function rho0thermal2modes(Nmaxa,Nmaxb,be)
   rho0diag = [0.0+0.0*im for i in 1:(Nmaxa+1)*(Nmaxb+1)]
   rho0 = diagm(rho0diag)
   for k in 1:(Nmaxa+1)
      rho0[k,k]=exp(-be*(k-1))*(1 - exp(-be))
   end
   return rho0
end



end


