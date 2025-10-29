module fisher
push!(LOAD_PATH, pwd())
export FisherInf
using LinearAlgebra

function FisherInf(rho, A; epsilon = 0.00001) #rho es la matriz densidad que necesitas, A el operador
   eigval, eigvec = eigen(Hermitian(rho))
   N = length(eigval)
   suma1 = 0
   suma2 = 0
   
  #Terminos diagonales
   for i in 1:N
       psii= eigvec[:,i]
       pis = eigval[i]
       expA = psii' * A * psii
       expA2= psii' * A^2 * psii
       varA = real(expA2 - expA^2)
       suma1 += pis * varA
   end
   
   #Terminos no diagonales
   for m in 1:N
        for n in 1:N
            if m != n
                pm = eigval[m]
                pn = eigval[n]
                denom = pm + pn
                if abs(denom) > epsilon  # evita división por cero
                    psim = eigvec[:, m]
                    psin = eigvec[:, n]
                    Mmn = psim' * A * psin
                    term = (pm * pn / abs(denom)) * abs2(Mmn)
                    suma2 += real(term)  # puede haber residuos complejos pequeños
                end
            end
        end
    end
    QFI = 4*suma1 - 8*suma2
    return QFI	
end

end #This is for closing module