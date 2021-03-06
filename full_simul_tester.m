function full_simul_tester(K,ep,Llx,sig,om,tf,dt,width)

    k0 = 1;
    [Om,cg,ad,anl] = param_maker(k0,om,sig);    
    Xmesh = -Llx:Llx/K:Llx-Llx/K;
    KT = 2*K;
    
    Kc = floor(KT/3);
    Kuc = KT - Kc + 1;
    Kc = Kc + 1;    
    
    Kvec = [ 0:K -K+1:-1 ]';
    Kmesh = pi/Llx*Kvec;
    
    rvecac = exp(-(Kmesh-k0).^2/(2*width^2*ep^2));
    
    arand = exp(2*pi*1i*rand(KT,1));
    uintsac = KT*sqrt(sqrt(pi)/(Llx*width*ep))*rvecac.*arand;   
    uintsac(Kc:Kuc) = 0;
    act_sol = afm_dno_solver(K,k0,ep,Llx,sig,om,Om,tf,dt,uintsac);
    
    figure(1)
    plot(fftshift(Kmesh),(ep/KT)^2*fftshift(abs(uintsac.*conj(uintsac))),'k-',fftshift(Kmesh)+k0,(ep/KT)^2*fftshift(abs(act_sol.*conj(act_sol))),'k--','LineWidth',2)
    
    figure(2)
    plot(Xmesh,ep*real(ifft(act_sol)),'k-','LineWidth',2)
    