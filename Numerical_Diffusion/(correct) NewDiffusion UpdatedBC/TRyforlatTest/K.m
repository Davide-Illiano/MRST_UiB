function sol = K(x,y,Nmod,KMean,varK,C1,C2,phi)

coeff = sqrt(varK*2/Nmod) ;

ak = coeff*sum( cos( (C1*x + C2*y)*(2*pi) + phi) ) ;

sol = KMean * exp(-varK/2) * exp(ak) ;

end