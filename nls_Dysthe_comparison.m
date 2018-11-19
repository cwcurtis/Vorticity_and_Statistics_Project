function nls_Dysthe_comparison(Llx,K,ep,tf,dt,om,sig,width,Nens)
    
    % [-Llx,Llx] is size of simulation
    % K is number of modes in simulation
    % ep is magnitude of perturbation in Dysthe equation
    % tf is length of time to let simulation run
    % dt is time step in RK4 solver for time simulation
    % k0 is wave number of carrier wave 
    % om is magnitude of vorticity 
    % sig is surface tension
    % Nens is number of ensemble members

    Llx = ep*Llx;
    
    disp('Spatial mesh size is')
    disp(Llx/K)
    
    nmax = round(tf/dt); % Step size for time integrator and number of time steps
    k0 = 1;
    disp('Carrier wave number is')
    disp(k0)
    
    [Om,cg,ad,anl] = param_maker(k0,om,sig);    
    
    KT = 2*K;
    
    Kc = floor(KT/3);
    Kuc = KT - Kc + 1;
    Kc = Kc + 1;    
    
    Kvec = [ 0:K -K+1:-1 ]';
    Kmesh = pi/Llx*Kvec;
        
    Xmesh = -Llx:Llx/K:Llx-Llx/K;
    
    uints = zeros(KT,Nens);
    nls_sol = zeros(KT,Nens);
    dysthe_sol = zeros(KT,Nens);
    
    rvec = exp(-(Kmesh).^2/(2*width^2));
    kvec = pi/Llx*(-K+1:K)';
    %deadinds = abs(ifftshift(kvec)) > sqrt(2)*width;
    
    for jj=1:Nens
        arand = 2*pi*1i*rand(KT,1);
        uints(:,jj) = KT*sqrt(sqrt(pi)/(Llx*width))*rvec.*exp(arand);    
        %uints(deadinds,jj) = 0;
        uints(Kc:Kuc,jj) = 0;
        %plot(Xmesh,abs(ifft(uints(:,jj))),'k','LineWidth',2)
        %pause
    end   
    
    aval = sqrt(sqrt(pi)/(Llx*width)*sum(rvec.^2));
    disp('Root mean square NLS amplitude')
    disp(aval)
    disp('Minimal Stable Envelope Width')
    disp(sqrt(2*anl/ad))
    
    parfor jj=1:Nens
        nls_sol(:,jj) = nls_solver_slider(K,Llx,nmax,ad,anl,dt,uints(:,jj));
        dysthe_sol(:,jj) = vor_Dysthe_solver_slider(K,Llx,nmax,ad,anl,cg,k0,Om,om,sig,ep,dt,uints(:,jj));
    end
    
    nls_spec_mean = ep^2*fftshift(abs(mean(nls_sol.*conj(nls_sol),2))/KT^2);
    dysthe_spec_mean = ep^2*fftshift(abs(mean(dysthe_sol.*conj(dysthe_sol),2))/KT^2);
    mean_int = ep^2*fftshift(abs(mean(uints.*conj(uints),2))/KT^2);
    
    
    plot(ep*kvec+k0,mean_int,'k-',ep*kvec+k0,nls_spec_mean,'k--',ep*kvec+k0,dysthe_spec_mean,'k-.','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$k$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\left<\left|\hat{\eta}(k,t_{f})\right|^{2}\right>$','Interpreter','LaTeX','FontSize',30)
    legend({'$Initial$','$NLS$','$Dysthe$'},'Interpreter','LaTeX','FontSize',30)
end

