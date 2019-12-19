function [avg,std,kts] = spec_comp(in_sol,kvec,dk)

    spec_plot = fftshift(abs(mean(in_sol.*conj(in_sol),2)));
    spec_prob = spec_plot/trap_comp(spec_plot,dk);    
    kavg = kvec.*spec_prob;    
    avg = trap_comp(kavg,dk);    
    spec_var = (kvec-avg).^2.*spec_prob;    
    std = sqrt(trap_comp(spec_var,dk));
    spec_kts = (kvec-avg).^4.*spec_prob;    
    kts = trap_comp(spec_kts,dk)/std^4;
    
    