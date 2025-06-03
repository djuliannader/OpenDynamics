module wigner
push!(LOAD_PATH, pwd())
export wignerf
export focktowf
export wigner_mix
#import body
#import norm
using LinearAlgebra
import build


function wignerf(psi,L,N,name::String)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
	 imax=length(psi)-imin
	 d=2*L/N
	 x=[-L+i*d for i in 1:(N-1)]
	 eps=0.001
	 if (abs(psi[imin])+abs(psi[imax]))>eps
	   println("*Please, consider a larger domain to see the Wigner function*")
	   return "Done"
         end
	 open(name,"w") do io
	 sumw=0.0
	 sumnw=0.0
	 sumnx=0.0
         for i in imin:imax
	   xinst=x[i]
	   for j in imin:imax
	     pinst=x[j]
	     sum=0.0+0*im
	     ie=1
	     for k in imin:imax
	        y=x[k]
	        sum=sum+(exp(2*im*pinst*y))*conj.(psi[i-iint+ie])*psi[i+iint-ie+1]*d
		ie=ie+1
	     end
	     w=sum/(pi)
	     println(io,xinst," ",pinst," ",round(real(w),digits=16))
	     sumw=sumw+d*d*w
	     sumnw=sumnw+d*d*abs(w)
	     sumnx=sumnx+d*d*w*(xinst*xinst)
           end
	 end
	 println("Go to file ",name," to see the wigner function")
	 #println("Volume of the wigner function: ",real(sumw))
	 #println("Volume of the negative region: ",real(sumnw)-1)
	 #println("Expectation value <x^2> : ",real(sumnx))
	 return [real(sumw),real(sumnw)-1]
	 end
	 end


function wignerf_component(psi,L,N)
	 pi=acos(-1)
         imin=trunc(Int64,N/4)
	 iint=imin
	 imax=length(psi)-imin
	 d=2*L/N
	 x=[-L+i*d for i in 1:(N-1)]
	 eps=0.001
	 wigner=[]
	 if (abs(psi[imin])+abs(psi[imax]))>eps
	   println("*Please, consider a larger domain to see the Wigner function*")
	   return "Done"
         end
	 #open(name,"w") do io
	 sumw=0.0
	 sumnw=0.0
	 sumnx=0.0
         for i in imin:imax
	   xinst=x[i]
	   for j in imin:imax
	     pinst=x[j]
	     sum=0.0+0*im
	     ie=1
	     for k in imin:imax
	        y=x[k]
	        sum=sum+(exp(2*im*pinst*y))*conj.(psi[i-iint+ie])*psi[i+iint-ie+1]*d
		ie=ie+1
	     end
	     w=sum/(pi)
	     wc = convert(Float64,real(w))
	     #println(io,xinst," ",pinst," ",)
	     #sumw=sumw+d*d*w
	     #sumnw=sumnw+d*d*abs(w)
	     #sumnx=sumnx+d*d*w*(xinst*xinst)
	   push!(wigner,[xinst,pinst,wc])  
           end
	 end
	 #end
         return wigner
	 end


function wigner_mix(wp,psis,L,N,K,name::String)
     winst=wignerf_component(psis[1],L,N)
     wres = [0.0 for i in 1:length(winst)]
     d=2*L/N
   for i in 1:length(wp)
     winst=wignerf_component(psis[i],L,N)
     wr = [winst[j][3] for j in 1:length(winst)]
     wres = wres  + wp[i]*wr
   end
   sumw = 0.0
   sumabs = 0.0
   open(name,"w") do io
      for i in 1:length(wres)
        println(io,winst[i][1]*K^(1/2)," ",winst[i][2]*K^(1/2)," ",wres[i])
	sumw=sumw+d*d*wres[i]
	sumabs=sumabs+d*d*abs(wres[i])
      end
   end
   res1=convert(Float64,real(sumw))
   res2=convert(Float64,real(sumabs))
   return [res1,res2-res1]
end

function focktowf(lfock,L,N)
  d=2*L/N
  Nfock=length(lfock)-1
  x=[-L+i*d for i in 1:(N-1)]
  wf=[]
  for i in 1:length(x)
    val=0.0 + 0.0*im
    for n in 0:Nfock
        val = val + lfock[n+1]*fockbasis(n,x[i])
    end
    append!(wf,val)
  end
  return wf
end

function hermite(x,n)
  nf = trunc(Int, n/2)
  sum=0.0
  for m in 0:nf
   sum = sum + ((-1)^m/(factorial(big(m))*factorial(big(n-2*m))))*((2*x)^(n-2*m))
  end
  ff = factorial(big(n))*sum
  ffc = convert(ComplexF64, ff)
  return ffc 
end

function fockbasis(n,x)
  fb = ((1/(2^n*factorial(big(n))*pi^(1/2))^(1/2)))*hermite(x,n)*exp(-x^2/2)
  fbc = convert(ComplexF64, fb)
  return fb
end

#N=400
#Nmax=20
#L=20
#smix=12
#wp=[0.1 for i in 1:smix]
#wfm=[]
#for i in 1:smix
#  cs = build.coherentstate(Nmax,0.1*i,0.1*i)
#  rho = cs*transpose(conj(cs))
#  evecs = eigvecs(rho)
#  comp = [evecs[i,Nmax+1] for i in 1:(Nmax+1)]
#  compwf = focktowf(comp,L,N)
#  append!(wfm,[compwf])
#end


#wigner_mix(wp,wfm,L,N,"output/wigmix_test.dat")



end