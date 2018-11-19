function ut = vor_Dysthe_solver_slider(K,Llx,nmax,ad,anl,cg,k0,Om,om,sig,ep,dt,uint)

    % Llx is domain size, i.e. solve on domain [-Llx,Llx]
    % sig is -1 for bright/focusing, +1 for dark/defocusing
    % tf is total simulation time
    % mu is magnitude of damping on envelope of dark initial condition
    Kvec = [ 0:K -K+1:-1 ]';
    Kmesh = pi/Llx*Kvec;    
    Dx = 1i*Kmesh;
    Dx2 = Dx.^2;
    s = sign(k0);
    L0i = (1 - 2*1i*ep*s*cg/(om-2*s*Om)*Dx);
    Lop = (L0i*1i*ad+ ep*sig/(om-2*s*Om)*Dx).*Dx2;
    KT = 2*K;
    Kc = floor(KT/3);
    Kuc = KT - Kc + 1;
    Kc = Kc + 1;
        
    cfac = 475/1440;
    c0 = 1427/1440*1/cfac;
    c1 = -133/240*1/cfac;
    c2 = 241/720*1/cfac;
    c3 = -173/1440*1/cfac;
    c4 = 9/480*1/cfac;

    ind_low = abs(Lop)<=.43/dt;
    ind_mid = abs(Lop)<=1.36/dt & abs(Lop)>.43/dt;
    ind_high = abs(Lop)>1.36/dt;

    Limid = (ones(KT,1)-cfac*dt*Lop).^(-1); % for AM6 method
    Lihigh = (ones(KT,1)-3*dt*Lop/4).^(-1); % for AM2* method

    ul = ind_low.*uint;
    um = ind_mid.*uint;
    uh = ind_high.*uint;

    ulnn3 = ul;
    ulnn2 = ul;
    ulnn1 = ul;

    umnn4 = um;
    umnn3 = um;
    umnn2 = um;
    umnn1 = um;

    uhnn1 = uh; 
    ut = uint;    
    nl = nonlinearity(ut,k0,cg,om,Om,ep,sig,ad,anl,Dx);

    nlnn3 = nl;
    nlnn2 = nl;
    nlnn1 = nl;

    for jj=1:nmax   
        % Compute the nonlinearity
        nl = nonlinearity(ut,k0,cg,om,Om,ep,sig,ad,anl,Dx);
        ab4 = dt*(55/24*nl - 59/24*nlnn1 + 37/24*nlnn2 - 3/8*nlnn3);
    
        % Compute the lows
        ab4l = dt*(55/24*ul - 59/24*ulnn1 + 37/24*ulnn2 - 3/8*ulnn3);
        ulnp1 = ul + Lop.*ab4l + ind_low.*ab4;
    
        ulnn3 = ulnn2;
        ulnn2 = ulnn1;
        ulnn1 = ul;
        ul = ulnp1;
         
        % Compute the mids
        csum = c0*um + c1*umnn1 + c2*umnn2 + c3*umnn3 + c4*umnn4;
        umnp1 = Limid.*(um + csum + ind_mid.*ab4) - csum;
                 
        umnn4 = umnn3;
        umnn3 = umnn2;
        umnn2 = umnn1;
        umnn1 = um;
        um = umnp1;
        
        % Compute the highs
        uhnp1 = Lihigh.*(uh + uhnn1/3 + ind_high.*ab4) - uhnn1/3;
    
        uhnn1 = uh;
        uh = uhnp1;
    
        nlnn3 = nlnn2;
        nlnn2 = nlnn1;
        nlnn1 = nl;
    
        % Put the solution back together
        ut = ul + um + uh;
        ut(Kc:Kuc) = 0;
    end
    
end