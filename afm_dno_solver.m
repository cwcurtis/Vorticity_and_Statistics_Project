function etanls = afm_dno_solver(K,k0,ep,Llx,sig,om,Om,tf,dt,uintsac)
    
    uintphys = ifft(uintsac);
    etan = ep*fft(2*real(uintphys));
    %etan(1) = 0;
    Qn = ep*Om*fft(-2*real(uintphys));
    Qn(1) = 0; %Enforce zero average in the surface velocity
    
    KT = 2*K;
    % Find the wave numbers to implement the 2/3 de-aliasing throughout
    
    Kmesh = pi/Llx*[0:K -K+1:-1]';
    Xmesh = (-Llx:Llx/K:Llx-Llx/K)';
        
    nmax = round(tf/dt);
    
    G0 = -1i*sign(Kmesh);
    L1 = -1i*Kmesh.*(1+sig.*Kmesh.^2);
    dti = 3*dt/4;
    
    Linvd = (1 + dti*om*G0 - dti^2*G0.*L1).^(-1);
    Linv11 = (1 + om*dti*G0).*Linvd;
    Linv12 = dti*G0.*Linvd;
    Linv21 = dti*L1.*Linvd;
    
    etanm1 = etan;
    Qnm1 = Qn;
    
    nln = dno_nonlinearity(K,etan,Qn,om,ep,sig,Kmesh);
    
    nlnm1 = nln;
    nlnm2 = nlnm1;
    nlnm3 = nlnm2;
    
    for jj=1:nmax
        
        nln = dno_nonlinearity(K,etan,Qn,om,ep,sig,Kmesh);
        
        nlvecn = 55/24*nln(1:KT) - 59/24*nlnm1(1:KT) + 37/24*nlnm2(1:KT) - 3/8*nlnm3(1:KT);    
        nlvecq = 55/24*nln(KT+1:2*KT) - 59/24*nlnm1(KT+1:2*KT) + 37/24*nlnm2(KT+1:2*KT) - 3/8*nlnm3(KT+1:2*KT);
        
        etatot = etan + etanm1/3 + dt*nlvecn;
        qtot = Qn + Qnm1/3 + dt*nlvecq;
        
        eta1 = Linv11.*etatot;
        eta2 = Linv12.*qtot;
        
        Q1 = Linv21.*etatot;
        Q2 = Linvd.*qtot;
        
        etanp1 = -etanm1/3 + eta1 + eta2;
        Qnp1 = -Qnm1/3 + Q1 + Q2;        
        
        etanm1 = etan;
        etan = etanp1;
        
        Qnm1 = Qn;
        Qn = Qnp1;
        
        nlnm3 = nlnm2;
        nlnm2 = nlnm1;
        nlnm1 = nln;
        
    end        
    
    etanx = real(ifft(1i*Kmesh.*etan));
    etanp = real(ifft(etan));
    theta = exp(1i*k0*Xmesh)*exp(1i*tf*Om);
    
    etanls = fft((etanp - 1i*etanx/k0)./(2*theta));