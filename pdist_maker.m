function pdist = pdist_maker(in_sol,dk)

    spec_plot = fftshift(abs(mean(in_sol.*conj(in_sol),2)));
    pdist = spec_plot/trap_comp(spec_plot,dk);    
    