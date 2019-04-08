function full_simul_tester(K,ep,Llx,sig,om,tf,dt,width,Nens)

    k0 = 10;
    [Om,cg,ad,anl] = param_maker(k0,om,sig);    
    Xmesh = -Llx:Llx/K:Llx-Llx/K;
    KT = 2*K;
    
    Kc = floor(KT/3);
    Kuc = KT - Kc + 1;
    Kc = Kc + 1;    
    
    Kvec = [ 0:K -K+1:-1 ]';
    Kmesh = pi/Llx*Kvec;
    
    rvecac = exp(-(Kmesh-k0).^2/(2*width^2*ep^2));
    uintsac = zeros(KT,Nens);
    
    for jj=1:Nens
        arand = exp(2*pi*1i*rand(KT,1));
        uintsac(:,jj) = KT*sqrt(sqrt(pi)/(Llx*width*ep))*rvecac.*arand;   
        uintsac(Kc:Kuc,jj) = 0;
    end
    
    parfor jj=1:Nens
        act_sol(:,jj) = afm_dno_solver(K,k0,ep,Llx,sig,om,Om,tf,dt,uintsac)
    end
    
    figure(1)
    plot(fftshift(Kmesh),ep^2*fftshift(abs(uintsac(:,1).*conj(uintsac(:,1)))),'k-',fftshift(Kmesh)+k0,fftshift(abs(act_sol(:,1).*conj(act_sol(:,1)))),'k--','LineWidth',2)
    