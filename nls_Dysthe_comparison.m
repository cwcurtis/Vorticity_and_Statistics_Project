function [nls_avg,nls_std,nls_kts,nls_bfi,dysthe_avg,dysthe_std,dysthe_kts,dysthe_bfi,int_avg,int_std,int_kts] = nls_Dysthe_comparison(Llx,K,ep,tf,dt,om,sig,width,k0,Nens)
    
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
    
    samp_times = linspace(0,tf,20);
    samp_inds = floor(samp_times/dt);
    
    tmat_nls = zeros(KT,Nens,length(samp_inds));
    tmat_dysthe = zeros(KT,Nens,length(samp_inds));
    
    parfor jj=1:Nens        
        tmat_nls(:,jj,:) = nls_solver(K,Llx,nmax,ad,anl,dt,samp_inds,uints(:,jj));
        tmat_dysthe(:,jj,:) = vor_Dysthe_solver(K,Llx,nmax,ad,anl,cg,k0,Om,om,sig,ep,dt,samp_inds,uints(:,jj));                        
    end
    
    [nls_avg_plot,nls_std_plot,nls_skw_plot,nls_kts_plot,nls_rms_plot] = spec_comp_time(tmat_nls,kvec,dk);
    [dysthe_avg_plot,dysthe_std_plot,dysthe_skw_plot,dysthe_kts_plot,dysthe_rms_plot] = spec_comp_time(tmat_dysthe,kvec,dk);
    [int_avg,int_std,int_kts] = spec_comp(uints,kvec,dk);
    
    nls_avg = nls_avg_plot(end);
    nls_std = nls_std_plot(end);
    nls_kts = nls_kts_plot(end);
    nls_rms = nls_rms_plot(end);
    
    dysthe_avg = dysthe_avg_plot(end);
    dysthe_std = dysthe_std_plot(end);
    dysthe_kts = dysthe_kts_plot(end);
    dysthe_rms = dysthe_rms_plot(end);
    
    nls_bfi_plot = abs(sqrt(anl/(2*ad)))*nls_rms_plot./nls_std_plot;
    dysthe_bfi_plot = abs(sqrt(anl/(2*ad)))*dysthe_rms_plot./dysthe_std_plot;
    
    nls_bfi = nls_bfi_plot(end);
    dysthe_bfi = dysthe_bfi_plot(end);
    
    %disp([nls_avg; dysthe_avg; int_avg])
    %disp([nls_std; dysthe_std; int_std])
        
    pdist_nls = pdist_maker(squeeze(tmat_nls(:,:,end)),dk);
    pdist_dysthe = pdist_maker(squeeze(tmat_dysthe(:,:,end)),dk);
    pdist_int = pdist_maker(uints,dk);
    
    figure(1)
    plot(kvec,pdist_int,'k:',kvec,pdist_nls,'k--',kvec,pdist_dysthe,'k-','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlim([-5,5])
    xlabel('$k$','Interpreter','LaTeX','FontSize',30)
    ylabel('$S(k,\tau_{f},\sigma)$','Interpreter','LaTeX','FontSize',30)
    %legend({'Initial','NLS','Dysthe'},'Interpreter','LaTeX','FontSize',30)
    
    figure(2)
    plot(samp_times,nls_avg_plot,'k--',samp_times,dysthe_avg_plot,'k-','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$\tau$','Interpreter','LaTeX','FontSize',30)
    ylabel('$m(\tau)$','Interpreter','LaTeX','FontSize',30)
    %legend({'NLS','Dysthe'},'Interpreter','LaTeX','FontSize',30)
        
    figure(3)
    plot(samp_times,nls_std_plot,'k--',samp_times,dysthe_std_plot,'k-','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$\tau$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\tilde{\sigma}(\tau)$','Interpreter','LaTeX','FontSize',30)
    %legend({'NLS','Dysthe'},'Interpreter','LaTeX','FontSize',30)
    
    figure(4)
    plot(samp_times,nls_skw_plot,'k--',samp_times,dysthe_skw_plot,'k-','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$\tau$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\tilde{s}(\tau)$','Interpreter','LaTeX','FontSize',30)
        
    figure(5)
    plot(samp_times,nls_kts_plot,'k--',samp_times,dysthe_kts_plot,'k-','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$\tau$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\tilde{k}(\tau)$','Interpreter','LaTeX','FontSize',30)
    %legend({'NLS','Dysthe'},'Interpreter','LaTeX','FontSize',30)
    
    figure(6)
    plot(samp_times,nls_bfi_plot,'k--',samp_times,dysthe_bfi_plot,'k-','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$\tau$','Interpreter','LaTeX','FontSize',30)
    ylabel('$BFI(\tau)$','Interpreter','LaTeX','FontSize',30)
    %legend({'NLS','Dysthe'},'Interpreter','LaTeX','FontSize',30)
    
end

