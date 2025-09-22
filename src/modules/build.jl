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
  #xop=(1/2^(1/2))*(aop+adop)
  #pop=(-im)*(1/2^(1/2))*(aop-adop)
  #apar=10.0
  #Ham = adop*aop+(1/2)*id
  Ham = -delta*adop*aop - ep*(adop^2+aop^2) + K*adop^2*aop^2
  #Ham = pop^2/2 - apar*xop^2 + xop^4
  return Ham
end

function HamiltonianOsc(Nmax,w0,x0,E0)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab)
  id = Diagonal(diagid)
  adop = transpose(aop)
  xop=(1/2^(1/2))*(aop+adop)
  pop=(-im)*(1/2^(1/2))*(aop-adop)
  Ham = pop^2/2 + (w0/2)*(xop-x0*id)^2 + id*E0
  return Ham
end

function HamiltonianQuartic(Nmax,a,b)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab)
  id = Diagonal(diagid)
  adop = transpose(aop)
  xop=(1/2^(1/2))*(aop+adop)
  pop=(-im)*(1/2^(1/2))*(aop-adop)
  Ham = pop^2/2 + a*xop^2 + b*xop^4
  return Ham
end

function anhilation(Nmax)
  diag  = [0.0 for i in 1:(Nmax+1)]
  diagid  = [1.0 for i in 1:(Nmax+1)]
  diaga = [0.0 for i in 1:Nmax]
  diagab = [sqrt(i) for i in 1:Nmax]
  aop = Tridiagonal(diaga, diag, diagab) 
  return aop
end

function coherentstate(Nmax,xc,pc)
  al=(1/2^(1/2))*(xc+im*pc)
  vstate = [0 for i in 1:(Nmax+1)]
  vstate[1] = 1.0
  aop = anhilation(Nmax)
  adop = transpose(aop)
  Dop= exp(al*adop - conj(al)*aop)
  csn = Dop*vstate
  return csn
end

function displacedfock(Nmax,xc,pc,k)
  al=(1/2^(1/2))*(xc+im*pc)
  vstate = [0 for i in 1:(Nmax+1)]
  vstate[k] = 1.0
  aop = anhilation(Nmax)
  adop = transpose(aop)
  Dop= exp(al*adop - conj(al)*aop)
  csn = Dop*vstate
  return csn
end

function generalCoherentStates(Nmax,k,xc,pc)		#k symbolizes the kth displaced state
  al=(1/2^(1/2))*(xc+im*pc) 		     		#alpha
  diagid  = [1.0 for i in 1:(Nmax+1)]
  aop = anhilation(Nmax)
  id = Diagonal(diagid)
  adop = transpose(aop)
  cs= coherentstate(Nmax,xc,pc)				#initialCoherentState
  gcs=(1/sqrt(factorial(big(k)))*(adop-conj(al)*id)^k*cs)
  return gcs
end

function initialrho(Nmax,xc,pc)
   cs = coherentstate(Nmax,xc,pc)
   cstp = transpose(conj(cs))
   rhoin = cs*cstp
   return rhoin
end

function initialrhocat(Nmax,alpha)
   xc = 2^(1/2)*real(alpha)
   pc = 2^(1/2)*imag(alpha)
   cat = (1/2^(1/2))*(coherentstate(Nmax,xc,pc)+coherentstate(Nmax,-xc,-pc))
   cattp = transpose(conj(cat))
   rhoin = cat*cattp
   return rhoin
end

function initialrhofock(Nmax,xc,pc,k)
   ds = displacedfock(Nmax,xc,pc,k)
   dtp = transpose(conj(ds))
   rhoin = ds*dtp
   return rhoin
end


function initialrhoGS(Ham)
  eigval, eigvec=eigen(Ham)
  gsvec=eigvec[:,1]
  gsvecad=transpose(conj(gsvec))
  rhogs = gsvec*gsvecad
  return ComplexF64.(rhogs)
end

function initialpsi(Nmax,xc,pc)
   cs = coherentstate(Nmax,xc,pc)
   return cs
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

function initialrhoquench(Nmax,epsilon0,Delta,K)
  HH = HamiltonianKerr(Nmax,Delta,epsilon0,K)
  eigvecsHH = eigvecs(HH)
  eigvalsHH = eigvals(HH)
  psi0 = [eigvecsHH[i,1]+0.0*im for i in 1:length(eigvalsHH)]
  psi0t = transpose(conj(psi0))
  rhogs = psi0*psi0t
  return rhogs
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


