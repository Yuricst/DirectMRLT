%% Construct problem data struct containing constants
%
%
% 
function data = get_problem_data(GM,LU,MU,thrust,mdot)
    data.MU = MU;
    data.LU = LU;
    data.VU = sqrt(GM/LU);
    data.TU = LU/data.VU;
    data.GM = 1.0;
    data.J2 = 1082.639e-06;
    
    data.c1 = thrust * (1/MU)*(data.TU^2/(1e3*LU));
    data.c2 = abs(mdot) *(data.TU/data.MU);
end