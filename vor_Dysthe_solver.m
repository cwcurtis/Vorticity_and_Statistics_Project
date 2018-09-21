function w = vor_Dysthe_solver(K,Llx,nmax,ad,anl,cg,k0,Om,om,sig,ep,dt,uint)

    % solves u_t = dk2Om*u_xx + a3*|u|^2*u 

    % Llx is domain size, i.e. solve on domain [-Llx,Llx]
    % sig is -1 for bright/focusing, +1 for dark/defocusing
    % tf is total simulation time
    % mu is magnitude of damping on envelope of dark initial condition
    Kvec = [ 0:K -K+1:-1 ]';
    Kmesh = pi/Llx*Kvec;    
    Dx = 1i*Kmesh;
    Dx2 = Dx.^2;
    s = sign(k0);
    %L0i = (1 + 2*1i*ep*s*cg/(om-2*s*Om)*Dx).^(-1);
    L0i = (1 - 2*1i*ep*s*cg/(om-2*s*Om)*Dx);
    %Lap = dt*L0i.*(1i*ad + ep*sig/(om-2*s*Om)*Dx).*Dx2;
    Lap = dt*(L0i*1i*ad+ ep*sig/(om-2*s*Om)*Dx).*Dx2;
    KT = 2*K;
    Kc = floor(KT/3);
    Kuc = KT - Kc + 1;
    Kc = Kc + 1;
    
    % Bright/Focusing initial conditions       
    w = uint; 
        
    % Begin setup of RK4 method
    E = exp(Lap/2);
    E2 = E.^2;
                
    % Solve NLS equation in time    
    for nn=1:nmax            
        a = dt*nonlinearity(w,k0,cg,om,Om,ep,sig,ad,anl,Dx);
        af = E.*(w+a/2);
        af(Kc:Kuc) = 0;
        
        b = dt*nonlinearity(af,k0,cg,om,Om,ep,sig,ad,anl,Dx);
        bf = E.*w + b/2;
        bf(Kc:Kuc) = 0;
        
        c = dt*nonlinearity(bf,k0,cg,om,Om,ep,sig,ad,anl,Dx);
        cf = E2.*w + E.*c;
        cf(Kc:Kuc) = 0;
        
        d = dt*nonlinearity(cf,k0,cg,om,Om,ep,sig,ad,anl,Dx);
        d(Kc:Kuc) = 0;
        
        w = E2.*w + (E2.*a + 2*E.*(b+c) + d)/6;    
        w(Kc:Kuc) = 0;
    end
    
end