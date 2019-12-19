function [nls_avg,nls_std,int_avg,int_std] = nls_stats_only(Llx,K,tf,dt,om,sig,width,k0,Nens)
    
    % [-Llx,Llx] is size of simulation
    % K is number of modes in simulation
    % ep is magnitude of perturbation in Dysthe equation
    % tf is length of time to let simulation run
    % dt is time step in RK4 solver for time simulation
    % k0 is wave number of carrier wave 
    % om is magnitude of vorticity 
    % sig is surface tension
    % Nens is number of ensemble members
    
    nmax = round(tf/dt); % Step size for time integrator and number of time steps
    [Om,cg,ad,anl] = param_maker(k0,om,sig);    
    
    KT = 2*K;
    dk = pi/Llx;
    
    Kc = floor(KT/3);
    Kuc = KT - Kc + 1;
    Kc = Kc + 1;    
    
    Kvec = [ 0:K -K+1:-1 ]';
    Kmesh = dk*Kvec;
    
    rvec = exp(-(Kmesh).^2/(2*width^2));
    kvec = dk*fftshift(Kvec);
    arand = exp(2*pi*1i*rand(KT,Nens));
    uints = KT*sqrt(sqrt(pi)/(Llx*width))*repmat(rvec,1,Nens).*arand;   
    uints(Kc:Kuc,:) = 0;
    
    %aval = sqrt(sqrt(pi)/(Llx*width)*sum(rvec.^2));
    %disp('Root mean square NLS amplitude')
    %disp(aval)
    tmat_nls = zeros(KT,Nens,nmax+1);
    
    parfor jj=1:Nens        
        tmat_nls(:,jj,:) = nls_solver(K,Llx,nmax,ad,anl,dt,uints(:,jj));        
    end
    
    [nls_avg_plot,nls_std_plot] = spec_comp_time(tmat_nls,kvec,dk,nmax);
    [int_avg,int_std] = spec_comp(uints,kvec,dk);
    
    nls_avg = nls_avg_plot(end);
    nls_std = nls_std_plot(end);
    
    disp([nls_avg; int_avg])
    disp([nls_std; int_std])
        
    pdist_nls = pdist_maker(squeeze(tmat_nls(:,:,end)),dk);
    pdist_int = pdist_maker(uints,dk);
    
    figure(2)
    plot(kvec,pdist_int,'k-',kvec,pdist_nls,'k--','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$k$','Interpreter','LaTeX','FontSize',30)
    ylabel('$S(k,\tau_{f},\sigma)$','Interpreter','LaTeX','FontSize',30)
    legend({'$Initial$','$NLS$'},'Interpreter','LaTeX','FontSize',30)
    
    figure(3)
    plot(dt*(0:nmax),nls_avg_plot,'k-','LineWidth',2)
    set(h,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    ylabel('Average','Interpreter','LaTeX','FontSize',30)
    
        
    figure(4)
    plot(dt*(0:nmax),nls_std_plot,'k-','LineWidth',2)
    set(h,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    ylabel('Standard Deviation','Interpreter','LaTeX','FontSize',30)
    
        
end

