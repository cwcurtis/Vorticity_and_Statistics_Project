function ifval = trap_comp(fvals,dx)

    %ifval = dx/2*(fvals(1) + fvals(end) + 2*sum(fvals(2:end-1)));
    
    ifval = dx/3*(fvals(1) + fvals(end-1) + 2*sum(fvals(3:2:end-3)) + 4*sum(fvals(2:2:end-2)));