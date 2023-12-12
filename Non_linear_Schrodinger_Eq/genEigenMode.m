function Psi = genEigenMode(n,l,r,phi)

Psi = sqrt(factorial((n-abs(l))/2)/(pi*factorial((n+abs(l))/2)))*r.^(abs(l)).*laguerreL((n-abs(l))/2,abs(l),r.^2).*exp(-r.^2/2).*exp(1i*l*phi);