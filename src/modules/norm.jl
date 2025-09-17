module norm
export normalizing, derivative2

function normalizing(psi,d)
    s=[abs2(i) for i in psi]
	 nfac=sum(d*s)
	 t=[i/(nfac)^(1/2) for i in psi]
	 return t
	 end



function derivative2(psi,h)
	n=trunc(Int,length(psi)/2)+1
	der=(psi[n+1]-2*psi[n]+psi[n-1])/h^2
        return der
	end

end