function rhs = dno_nonlinearity(K,eta,Q,om,ep,sig,Kmesh)

KT = 2*K;
% Find the wave numbers to implement the 2/3 de-aliasing throughout
Kc = floor(2*K/3);
Kuc = KT-Kc+1;
Kc = Kc+1;

Dx = 1i*Kmesh;
Dx2 = Dx.^2;
%Dx3 = Dx2.*Dx;
Hop = 1i*sign(Kmesh);

eta(Kc:Kuc)=0;
Q(Kc:Kuc)=0;

etax = real(ifft(Dx.*eta)); 
etaxp2 = etax.*etax;
etaxx = real(ifft(Dx2.*eta)); 
etap = real(ifft(eta));
etap2 = etap.^2;
etap3 = etap2.*etap;

Qp = real(ifft(Q));
Qpn = Qp + om*etap;
Qpn2 = Qpn.^2;

G0 = -real(ifft(Hop.*Q));
G1 = real(ifft( Dx.*( Hop.*fft(etap.*G0) - fft(om/2*etap.^2 + etap.*Qp) )));
G2 = real(ifft( Dx.*Hop.*fft( etap.*G1 ) + Dx2.*( fft( etap2.*G0/2 ) + Hop.*fft(om/3*etap3 + etap2.*Qp/2) )));
%G3 = real(ifft( Dx.*Hop.*fft( etap.*G2 ) + Dx2.*fft(etap2.*G1/2) + Dx3.*( Hop.*fft( -etap3.*G0/6 )+ fft(om/8*etap3.*etap + etap3.*Qp/6) )));
%rhs1 = fft( ep*G1 + ep^2*G2 + ep^3*G3 );
%rhs2 = -om*rhs1 + Dx.*fft( ep/2*( -Qpn2 + G0.^2 + 2*ep*G0.*G1 + ep^2*(2*G0.*G2 + G1.^2) ) ...
%                         + ep^2*( etax.*Qpn.*(G0 + ep*G1) - 3/2*sig*etaxp2.*etaxx ) ...
%                         + ep^3*etaxp2.*( Qpn2 - G0.^2 )/2 );

rhs1 = fft( ep*G1 + ep^2*G2 );
rhs2 = -om*rhs1 + Dx.*fft( ep/2*( -Qpn2 + G0.^2 + 2*ep*G0.*G1 ) ...
                         + ep^2*( etax.*Qpn.*G0 - 3/2*sig*etaxp2.*etaxx ) );
                         

rhs1(Kc:Kuc) = 0;
rhs2(Kc:Kuc) = 0;

rhs = [rhs1;rhs2];
