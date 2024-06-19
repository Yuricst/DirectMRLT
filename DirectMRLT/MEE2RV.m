function RV = MEE2RV(mu,MEE)
    % unpack state
    p = MEE(:,1);
    f = MEE(:,2);
    g = MEE(:,3);
    h = MEE(:,4);
    k = MEE(:,5);
    L = MEE(:,6);

    % define variables as in eqn (6.46-6.41)
    q = 1+f.*cos(L)+g.*sin(L);
    r = p./q;
    alpha_sq = h.^2-k.^2;
    chi = sqrt(h.^2+k.^2);
    s_sq = 1+chi.^2;
    
    % Cartesian state (eqn 6.42-6,43)
    r_x = (r./s_sq).*(cos(L)+alpha_sq.*cos(L)+2*h.*k.*sin(L));
    r_y = (r./s_sq).*(sin(L)-alpha_sq.*sin(L)+2*h.*k.*cos(L));
    r_z = (2*r./s_sq).*(h.*sin(L)-k.*cos(L));
    v_x = -(1./s_sq).*sqrt(mu./p).*(sin(L)+alpha_sq.*sin(L)-2*h.*k.*cos(L)+g-2*f.*h.*k+alpha_sq.*g);
    v_y = -(1./s_sq).*sqrt(mu./p).*(-cos(L)+alpha_sq.*cos(L)+2*h.*k.*sin(L)-f+2*g.*h.*k+alpha_sq.*f);
    v_z = (2./s_sq).*sqrt(mu./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);
    RV = [r_x r_y r_z v_x v_y v_z];
end