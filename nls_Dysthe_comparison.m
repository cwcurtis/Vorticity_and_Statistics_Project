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
    Llxf = Llx;
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
    Kmeshf = pi/Llxf*Kvec;
        
    Xmesh = -Llx:Llx/K:Llx-Llx/K;
    
    uints = zeros(KT,Nens);
    uintsac = zeros(KT,Nens);
    nls_sol = zeros(KT,Nens);
    dysthe_sol = zeros(KT,Nens);
    act_sol = zeros(KT,Nens);
    
    rvec = exp(-(Kmesh).^2/(2*width^2));
    rvecac = exp(-(Kmeshf-k0).^2/(2*width^2*ep^2));
    kvec = pi/Llx*(-K+1:K)';
    %deadinds = abs(ifftshift(kvec)) > sqrt(2)*width;
    
    for jj=1:Nens
        arand = exp(2*pi*1i*rand(KT,1));
        uints(:,jj) = KT*sqrt(sqrt(pi)/(Llx*width))*rvec.*arand;   
        uintsac(:,jj) = KT*sqrt(sqrt(pi)/(Llxf*width*ep))*rvecac.*arand;   
        uints(Kc:Kuc,jj) = 0;
        uintsac(Kc:Kuc,jj) = 0;
    end   
    
    aval = sqrt(sqrt(pi)/(Llx*width)*sum(rvec.^2));
    disp('Root mean square NLS amplitude')
    disp(aval)
    disp('Minimal Stable Envelope Width')
    disp(sqrt(2*anl/ad))
    tic
    %Since these work over such different time scales.  
    parfor jj=1:Nens
        nls_sol(:,jj) = nls_solver(K,Llx,nmax,ad,anl,dt,uints(:,jj));
        dysthe_sol(:,jj) = vor_Dysthe_solver(K,Llx,nmax,ad,anl,cg,k0,Om,om,sig,ep,dt,uints(:,jj));        
    end    
    parfor jj=1:Nens
       act_sol(:,jj) = afm_dno_solver(K,k0,ep,Llxf,sig,om,Om,tf/ep^2,dt,uintsac(:,jj)) 
    end
    toc
    
    ekT = (ep/KT)^2;
    nls_spec_mean = ekT*fftshift(abs(mean(nls_sol.*conj(nls_sol),2)));
    dysthe_spec_mean = ekT*fftshift(abs(mean(dysthe_sol.*conj(dysthe_sol),2)));
    mean_int = ekT*fftshift(abs(mean(uints.*conj(uints),2)));
    mean_act = ekT*fftshift(abs(mean(act_sol.*conj(act_sol),2)));
    
    dk = ep*pi/Llx;
    
    nls_spec_mean = nls_spec_mean/(dk/2*(nls_spec_mean(1)+nls_spec_mean(end) + 2*sum(nls_spec_mean(2:end-1))));
    dysthe_spec_mean = dysthe_spec_mean/(dk/2*(dysthe_spec_mean(1)+dysthe_spec_mean(end) + 2*sum(dysthe_spec_mean(2:end-1))));
    mean_act = mean_act/(dk/2*(mean_act(1)+mean_act(end) + 2*sum(mean_act(2:end-1))));
    mean_int = mean_int/(dk/2*(mean_int(1)+mean_int(end) + 2*sum(mean_int(2:end-1))));
    
    plot(ep*kvec+k0,mean_int,'k-',ep*kvec+k0,nls_spec_mean,'k--',ep*kvec+k0,dysthe_spec_mean,'k-.',ep*kvec+k0,mean_act,'k:','LineWidth',2)
    %plot(ep*kvec+k0,mean_int,'k-',ep*kvec+k0,nls_spec_mean,'k--',ep*kvec+k0,mean_act,'k:','LineWidth',2)
    %plot(ep*kvec+k0,mean_int,'k-',ep*kvec+k0,nls_spec_mean,'k--',ep*kvec+k0,dysthe_spec_mean,'k-.','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$k$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\left<\left|\hat{\eta}(k,t_{f})\right|^{2}\right>$','Interpreter','LaTeX','FontSize',30)
    legend({'$Initial$','$NLS$','$Dysthe$','$Full$'},'Interpreter','LaTeX','FontSize',30)
    %legend({'$Initial$','$NLS$','$Dysthe$'},'Interpreter','LaTeX','FontSize',30)
end

