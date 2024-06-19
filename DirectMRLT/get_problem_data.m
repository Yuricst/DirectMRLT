%% Construct problem data struct containing constants
%
% Yuri Shimane, yuri.shimane@gatech.edu
%   Created    : 2024/06/19
%   Last edits : 2024/06/19
%
% This function computes the canonical scales for the problem according to
% the given LU and GM, by defining `VU = sqrt(GM/LU)` and `TU = LU/VU`.
% The thrust and mass-flow rate are then stored in canonical scales as
% `c1` and `c2` respectively. 
% 
% Inputs:
%   GM      : Gravitational constant of the central body [km^3/s^2]
%   LU      : Distance unit [km]
%   MU      : Mass unit [kg]
%   thrust  : Thrust [N]
%   mdot    : Mass flow rate [kg/s]
%   J2      : (optional) J2 constant
%
% Outputs:
%   data    : Struct containing problem constants
% 
function data = get_problem_data(GM,LU,MU,thrust,mdot,options)
    arguments
        GM (1,1) double {mustBeNumeric}
        LU (1,1) double {mustBeNumeric}
        MU (1,1) double {mustBeNumeric}
        thrust (1,1) double {mustBeNumeric}
        mdot (1,1) double {mustBeNumeric}
        % optional arguments
        options.J2 double = 1082.639e-06
    end

    data.MU = MU;
    data.LU = LU;
    data.VU = sqrt(GM/LU);
    data.TU = LU/data.VU;
    data.GM = 1.0;
    data.J2 = options.J2;
    
    data.c1 = thrust * (1/MU)*(data.TU^2/(1e3*LU));
    data.c2 = abs(mdot) *(data.TU/data.MU);
end