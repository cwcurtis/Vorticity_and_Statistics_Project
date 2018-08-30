function [Om,cg,ad,anl] = param_maker(k0,om,sig)

s = sign(k0);

Om = .5*(s*om + sqrt(om^2 + 4*abs(k0)*(1+sig*k0^2)));

cg = (1+3*sig*k0^2)/(2*s*Om-om);

ad = (cg^2 - 3*abs(k0)*sig)/(2*Om - s*om);

av = (cg*k0*s - 2*Om*s)*om^4 + k0*(4*k0^2*s*sig + 2*Om*cg - s)*om^3 ...
     + k0*(16*cg*k0^3*sig - 8*Om*k0^2*sig + 10*cg*k0 - 6*Om)*om^2 ...
     - k0^2*(15*Om*cg*k0^2*s*sig - 16*k0^4*sig^2 - 24*k0^2*sig - 2)*om ...
     + k0^3*(2*cg*k0^4*s*sig^2 + cg*k0^2*s*sig - 15*Om*k0*s*sig + 8*cg*s);

 anl = k0/((2*sign(k0)*Om-om)*(1+cg*om)*(4*Om^2-sign(k0)*(2*k0*(1+4*sig*k0^2)+2*om*Om)))*(sign(k0)*k0^3*(8 + sig*k0^2 + 2*(sig*k0^2)^2)+om*av);
 
%disp(ad)
%disp(anl)