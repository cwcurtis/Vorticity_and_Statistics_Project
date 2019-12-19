function sol = nls_solver(K,Llx,nmax,ad,anl,dt,samp_inds,uint)

    % solves u_t = dk2Om*u_xx + a3*|u|^2*u 

    % Llx is domain size, i.e. solve on domain [-Llx,Llx]
    % sig is -1 for bright/focusing, +1 for dark/defocusing
    % tf is total simulation time
    % mu is magnitude of damping on envelope of dark initial condition
   
        Kvec = [ 0:K -K+1:-1 ]';
        Kmesh = pi/Llx*Kvec;        
        Dx = 1i*Kmesh; % Fourier transform of first derivative
        Dx2 = Dx.^2; % Fourier transform of second derivative                       
        KT = 2*K;
        Kc = floor(KT/3);
        Kuc = KT - Kc + 1;
        Kc = Kc + 1;
        
        w = uint;             
        
        sol = zeros(KT,length(samp_inds));
        sol(:,1) = w;
    % Begin setup of RK4 method

        Lap = 1i*dt*ad*Dx2;
        E = exp(Lap/2);
        E2 = E.^2;
                
    % Solve NLS equation in time
        cnt = 2;
        for nn=1:nmax            
            wp = ifft(w);
            a = 1i*dt*anl*fft(wp.^2.*conj(wp));
            a(Kc:Kuc) = 0;
            
            ap = ifft(E.*(w+a/2));
            b = 1i*dt*anl*fft(ap.^2.*conj(ap));
            b(Kc:Kuc) = 0;
            
            bp = ifft(E.*w + b/2);
            c = 1i*dt*anl*fft(bp.^2.*conj(bp));
            c(Kc:Kuc) = 0;
            
            cp = ifft(E2.*w + E.*c);
            d = 1i*dt*anl*fft(cp.^2.*conj(cp));
            d(Kc:Kuc) = 0;
            
            w = E2.*w + (E2.*a + 2*E.*(b+c) + d)/6;                                        
            w(Kc:Kuc) = 0;
            if nn == samp_inds(cnt)
                sol(:,cnt) = w;
                cnt = cnt + 1;
            end
        end
        
end