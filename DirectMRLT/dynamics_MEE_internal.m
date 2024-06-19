%% MEE dynamics - internal version
% 
% Yuri Shimane, yuri.shimane@gatech.edu
%   Created    : 2024/06/19
%   Last edits : 2024/06/19
%
% Syntax:  
%          [dx] = Dynamics(x,u,params,t,vdat)	(Dynamics Only)
%          [dx,g_eq] = Dynamics(x,u,params,t,vdat)   (Dynamics and Eqaulity Path Constraints)
%          [dx,g_neq] = Dynamics(x,u,params,t,vdat)   (Dynamics and Inqaulity Path Constraints)
%          [dx,g_eq,g_neq] = Dynamics(x,u,params,t,vdat)   (Dynamics, Equality and Ineqaulity Path Constraints)
% 
function [dx,g_eq,g_neq] = dynamics_MEE_internal(x,u,params,t,vdat)

    % extract parameters
    GM = vdat.GM;
    c1 = vdat.c1;
    c2 = vdat.c2;

    % unpack state
    p = x(:,1);
    f = x(:,2);
    g = x(:,3);
    h = x(:,4);
    k = x(:,5);
    L = x(:,6);
    m = x(:,7);
    u_r = u(:,1);
    u_t = u(:,2);
    u_n = u(:,3);
    tau = u(:,4);

    % acceleration magnitudes
    accel_thrust_mag = c1 * tau ./ m;

    % aux. variables
    cosL = cos(L);
    sinL = sin(L);
    w = 1 + f.*cosL + g.*sinL;
    hsinL_kcosL = h.*sinL - k.*cosL;
    sqrt_p_GM = sqrt(p/GM);
    s2 = 1 + h.^2 + k.^2;

    % collect acceleration histories
    accel_r = accel_thrust_mag .* u_r;
    accel_t = accel_thrust_mag .* u_t;
    accel_n = accel_thrust_mag .* u_n;

    % terms inside B-matrix
    B12 = sqrt_p_GM .* 2.*p./w;

    B21 = sqrt_p_GM .* sinL;
    B22 = sqrt_p_GM .* ((1+w).*cosL + f)./w;
    B23 = -sqrt_p_GM .* g./w .* hsinL_kcosL;

    B31 = -sqrt_p_GM .* cosL;
    B32 = sqrt_p_GM .* ((1+w).*sinL + g)./w;
    B33 = sqrt_p_GM .* f./w .* hsinL_kcosL;

    B43 = sqrt_p_GM ./ w .* s2/2 .* cosL;

    B53 = sqrt_p_GM ./ w .* s2/2 .* sinL;

    B63 = sqrt_p_GM ./ w .* hsinL_kcosL;

    % state derivatives
    dp = B12 .* accel_t;
    df = B21 .* accel_r + B22 .* accel_t + B23 .* accel_n;
    dg = B31 .* accel_r + B32 .* accel_t + B33 .* accel_n;
    dh = B43 .* accel_n;
    dk = B53 .* accel_n;
    dL = B63 .* accel_n + sqrt(GM./p.^3)./w.^2;
    dm = -c2 .* tau;

    % Return variables
    dx = [dp df dg dh dk dL dm];
    g_eq = [u_r.^2 + u_t.^2 + u_n.^2 - 1];
    g_neq = []; %[u_r.^2 + u_t.^2 + u_n.^2];
end