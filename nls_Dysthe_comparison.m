function [nls_avg,std_nls,dys_avg,std_dys,int_avg,std_int] = nls_Dysthe_comparison(Llx,K,ep,tf,dt,om,sig,width,k0,Nens)
    
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
    
    %disp('Spatial mesh size is')
    %disp(Llx/K)
    
    nmax = round(tf/dt); % Step size for time integrator and number of time steps
    
    %disp('Carrier wave number is')
    %disp(k0)
    
    [Om,cg,ad,anl] = param_maker(k0,om,sig);    
    
    KT = 2*K;
    
    Kc = floor(KT/3);
    Kuc = KT - Kc + 1;
    Kc = Kc + 1;    
    
    Kvec = [ 0:K -K+1:-1 ]';
    Kmesh = pi/Llx*Kvec;
    Kmeshf = pi/Llxf*Kvec;
        
    %Xmesh = -Llx:Llx/K:Llx-Llx/K;
    
    uints = zeros(KT,Nens);
    uintsac = zeros(KT,Nens);
    nls_sol = zeros(KT,Nens);
    dysthe_sol = zeros(KT,Nens);
    %act_sol = zeros(KT,Nens);
    
    rvec = exp(-(Kmesh).^2/(2*width^2));
    rvecac = exp(-(Kmeshf-k0).^2/(2*width^2*ep^2));
    kvec = pi/Llx*(-K+1:K)';
    deadinds = abs(Kmesh) > 2*k0;
    %deadindsac = logical(1 - (Kmeshf>(k0-2*k0*ep)).*(Kmeshf<(k0+2*k0*ep)));
    for jj=1:Nens
        arand = exp(2*pi*1i*rand(KT,1));
        uints(:,jj) = KT*sqrt(sqrt(pi)/(Llx*width))*rvec.*arand;   
        uintsac(:,jj) = KT*sqrt(sqrt(pi)/(Llxf*width*ep))*rvecac.*arand;   
        %uints(Kc:Kuc,jj) = 0;
        uints(deadinds,jj) = 0;
        %uintsac(deadindsac,jj) = 0;
    end   
    
    %aval = sqrt(sqrt(pi)/(Llx*width)*sum(rvec.^2));
    %disp('Root mean square NLS amplitude')
    %disp(aval)
    %Since these work over such different time scales.  
    parfor jj=1:Nens
        nls_sol(:,jj) = nls_solver(K,Llx,nmax,ad,anl,dt,uints(:,jj));
        dysthe_sol(:,jj) = vor_Dysthe_solver(K,Llx,nmax,ad,anl,cg,k0,Om,om,sig,ep,dt,uints(:,jj));        
    end    
    %parfor jj=1:Nens
    %   act_sol(:,jj) = afm_dno_solver(K,k0,ep,Llxf,sig,om,Om,tf/ep^2,dt,uintsac(:,jj)) 
    %end
    ekT = (ep/KT)^2;
    nls_spec_plot = ekT*fftshift(abs(mean(nls_sol.*conj(nls_sol),2)));
    dysthe_spec_plot = ekT*fftshift(abs(mean(dysthe_sol.*conj(dysthe_sol),2)));
    mean_int_plot = ekT*fftshift(abs(mean(uints.*conj(uints),2)));
    %mean_act_plot = ekT*fftshift(abs(mean(act_sol.*conj(act_sol),2)));
    
    dk = ep*pi/Llx;
    rkvec = ep*kvec+k0;
    nls_spec_mean = nls_spec_plot/(dk/2*(nls_spec_plot(1)+nls_spec_plot(end)+2*sum(nls_spec_plot(2:end-1))));
    dysthe_spec_mean = dysthe_spec_plot/(dk/2*(dysthe_spec_plot(1)+dysthe_spec_plot(end)+2*sum(dysthe_spec_plot(2:end-1))));
    mean_int = mean_int_plot/(dk/2*(mean_int_plot(1)+mean_int_plot(end)+2*sum(mean_int_plot(2:end-1))));
    %mean_act = mean_act_plot/(dk/2*(mean_act_plot(1)+mean_act_plot(end)+2*sum(mean_act_plot(2:end-1))));
    
    knls = rkvec.*nls_spec_mean;
    kdys = rkvec.*dysthe_spec_mean;
    kint = rkvec.*mean_int;    
    %kact = rkvec.*mean_act;
    
    nls_avg = dk/2*(knls(1)+knls(end)+2*sum(knls(2:end-1)));
    dys_avg = dk/2*(kdys(1)+kdys(end)+2*sum(kdys(2:end-1)));
    int_avg = dk/2*(kint(1)+kint(end)+2*sum(kint(2:end-1)));
    %act_avg = dk/2*(kact(1)+kact(end)+2*sum(kact(2:end-1)));
    
    stnls = (rkvec-nls_avg).^2.*nls_spec_mean;
    stdys = (rkvec-dys_avg).^2.*dysthe_spec_mean;
    stint = (rkvec-int_avg).^2.*mean_int;
    %stact = (rkvec-act_avg).^2.*mean_act;
    std_nls = dk/2*(stnls(1)+stnls(end)+2*sum(stnls(2:end-1)));
    std_dys = dk/2*(stdys(1)+stdys(end)+2*sum(stdys(2:end-1)));
    std_int = dk/2*(stint(1)+stint(end)+2*sum(stint(2:end-1)));
    %std_act = dk/2*(stact(1)+stact(end)+2*sum(stact(2:end-1)));
    disp([nls_avg; dys_avg; int_avg])
    disp([std_nls; std_dys; std_int])
    
    %plot(rkvec,mean_int_plot,'k-',rkvec,nls_spec_plot,'k--',rkvec,dysthe_spec_plot,'k-.','LineWidth',2)
    %plot(rkvec,mean_int_plot,'k-',rkvec,nls_spec_plot,'k--',rkvec,dysthe_spec_plot,'k-.',rkvec,mean_act_plot,'k:','LineWidth',2)
    %h = set(gca,'FontSize',30);
    %set(h,'Interpreter','LaTeX')
    %xlabel('$k$','Interpreter','LaTeX','FontSize',30)
    %ylabel('$\left<\left|\hat{\eta}(k,t_{f})\right|^{2}\right>$','Interpreter','LaTeX','FontSize',30)
    %legend({'$Initial$','$NLS$','$Dysthe$'},'Interpreter','LaTeX','FontSize',30)
    %legend({'$Initial$','$NLS$','$Dysthe$','$Full$'},'Interpreter','LaTeX','FontSize',30)    
end

