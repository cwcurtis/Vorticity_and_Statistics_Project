function [avg,std,skw,kts,rms] = spec_comp_time(in_sol,kvec,dk)
    
    [KT,~,nmax] = size(in_sol);
    avg = zeros(nmax,1);
    std = zeros(nmax,1);
    skw = zeros(nmax,1);
    kts = zeros(nmax,1);
    rms = zeros(nmax,1);
    for ll = 1:nmax
        csol = squeeze(in_sol(:,:,ll));
        spec_plot = fftshift(abs(mean(csol.*conj(csol),2)));
        rms(ll) = 1/KT*sqrt(sum(spec_plot));
        spec_prob = spec_plot/trap_comp(spec_plot,dk);    
        kavg = kvec.*spec_prob;    
        avg(ll) = trap_comp(kavg,dk);    
        spec_var = (kvec-avg(ll)).^2.*spec_prob;    
        std(ll) = sqrt(trap_comp(spec_var,dk));    
        spec_skw = (kvec-avg(ll)).^3.*spec_prob;
        skw(ll) = trap_comp(spec_skw,dk)/std(ll)^3;
        spec_kts = (kvec-avg(ll)).^4.*spec_prob;    
        kts(ll) = 1/3*trap_comp(spec_kts,dk)/std(ll)^4-1;
    end
    
    