%% Orbit transfer problem in MEE
%
% Yuri Shimane, yuri.shimane@gatech.edu
%   Created    : 2024/06/19
%   Last edits : 2024/06/19
%

function [problem, guess] = MEEOrbitTransferProblem(...
    data,MEE_0,MEE_F,m0,t0,tf_bounds,mf_bounds,ICLOCSsettings,options)

    arguments
        data (1,1) struct
        MEE_0 (1,6) double {mustBeNumeric}
        MEE_F (1,6) double {mustBeNumeric}
        m0 (1,1) double {mustBeNumeric}
        t0 (1,1) double {mustBeNumeric}
        tf_bounds (1,2) double {mustBeNumeric}
        mf_bounds (1,2) double {mustBeNumeric}
        ICLOCSsettings (1,1) function_handle
        % optional arguments
        options.objective (1,1) string {mustBeMember(options.objective,["mf","tof"])} = "mf"
        options.max_rev (1,1) double = 100
        options.p_bounds (1,2) double {mustBeNumeric} = ...
            [min(MEE_0(1),MEE_F(1))*0.5 max(MEE_0(1),MEE_F(1))*1.5]
    end

    % store parameters for dynamics
    problem.data = data;
    problem.data.objective = options.objective;

    % Plant model name, used for Adigator
    InternalDynamics = @dynamics_MEE_internal;
    SimDynamics = @dynamics_MEE_internal;

    % Analytic derivative files (optional)
    problem.analyticDeriv.gradCost=[];
    problem.analyticDeriv.hessianLagrangian=[];
    problem.analyticDeriv.jacConst=[];

    % Settings file
    problem.settings = ICLOCSsettings;

    % Initial Time. t0<tf
    problem.time.t0_min = t0;
    problem.time.t0_max = t0;
    guess.t0 = t0;

    % Final time. Let tf_min=tf_max if tf is fixed.
    problem.time.tf_min = tf_bounds(1);     
    problem.time.tf_max = tf_bounds(2); 
    guess.tf = (tf_bounds(1) + tf_bounds(2))/2;

    % Parameters bounds. pl =< p <= pu
    problem.parameters.pl = [];
    problem.parameters.pu = [];
    guess.parameters = [];

    % Initial conditions for system.
    problem.states.x0 = horzcat(MEE_0, m0);

    % Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
    problem.states.x0l = horzcat(MEE_0, m0);
    problem.states.x0u = horzcat(MEE_0, m0);

    % State bounds. xl =< x <= xu
    problem.states.xl = [options.p_bounds(1) -1 -1 -1 -1 MEE_0(6) mf_bounds(1)];
    problem.states.xu = [options.p_bounds(2) 1 1 1 1 2*pi*options.max_rev m0];

    % State error bounds
    problem.states.xErrorTol_local    = [1 1 1 1 1 1 1];
    problem.states.xErrorTol_integral = [1 1 1 1 1 1 1];

    % State constraint error bounds
    problem.states.xConstraintTol = [1 1 1 1 1 1 1];

    % Terminal state bounds. xfl=< xf <=xfu
    problem.states.xfl = [MEE_F(1:5) MEE_0(6) mf_bounds(1)];
    problem.states.xfu = [MEE_F(1:5) 2*pi*options.max_rev mf_bounds(2)];

    % Number of control actions N 
    % Set problem.inputs.N=0 if N is equal to the number of integration steps.  
    % Note that the number of integration steps defined in settings.m has to be divisible 
    % by the  number of control actions N whenever it is not zero.
    problem.inputs.N = 0;       
        
    % Input bounds
    problem.inputs.ul = [-1 -1 -1 0];
    problem.inputs.uu = [ 1  1  1 1];

    problem.inputs.u0l = [-1 -1 -1 0];
    problem.inputs.u0u = [ 1  1  1 1];

    % Input constraint error bounds
    problem.inputs.uConstraintTol = [1 1 1 1];

    % Choose the set-points if required
    problem.setpoints.states = [];
    problem.setpoints.inputs = [];

    % Bounds for path constraint function gl =< g(x,u,p,t) =< gu
    problem.constraints.ng_eq = 1;
    problem.constraints.gTol_eq = [1];

    problem.constraints.gl=[];
    problem.constraints.gu=[];
    problem.constraints.gTol_neq=[];

    problem.constraints.bl=[];
    problem.constraints.bu=[];
    problem.constraints.bTol=[];

    % Obtain guess of states and input sequences with ode solve
    [guess.time,guess.states] = ode45(@(t,x) dynamics_MEE_initialguess(t,x,SimDynamics,problem.data), linspace(guess.t0,guess.tf,1000), problem.states.x0);
    guess.inputs(:,1) = 0 * ones(size(guess.time));  % -sind(0)*ones(size(guess.time));
    guess.inputs(:,2) = 1 * ones(size(guess.time));  %  cosd(0)*cosd(30)*ones(size(guess.time));
    guess.inputs(:,3) = 0 * ones(size(guess.time));  % -cosd(0)*sind(30)*ones(size(guess.time));
    guess.inputs(:,4) = 1 * ones(size(guess.time));  %  cosd(0)*cosd(30)*ones(size(guess.time));
    
    % Get function handles and return to Main.m
    problem.state_names = ['p','f','g','h','k','L','m'];
    problem.input_names = ['u_r','u_t','u_h','tau'];
    problem.data.InternalDynamics   = InternalDynamics;
    problem.data.functionfg         = @fg;
    problem.data.plantmodel         = func2str(InternalDynamics);
    problem.functions               = {@L,@E,@f,@g,@avrc,@b};
    problem.sim.functions           = SimDynamics;
    problem.sim.inputX              = [];
    problem.sim.inputU              = 1:length(problem.inputs.ul);
    problem.functions_unscaled      = {@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc,@b_unscaled};
    problem.data.functions_unscaled = problem.functions_unscaled;
    problem.data.ng_eq              = problem.constraints.ng_eq;
    problem.constraintErrorTol      = [problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];
end


function dx = dynamics_MEE_initialguess(t,x,dyn_func,data)

% initial guess control
u=[0 1 0 1];
params = [];

% Evaluate ODE right-hand side
dx = dyn_func(x',u,params,t,data)';
end


function stageCost = L_unscaled(x,xr,u,ur,p,t,vdat)
    % L_unscaled - Returns the stage cost.
    % The function must be vectorized and
    % xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
    % variable)
    % 
    % Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
    %
    % Inputs:
    %    x  - state vector
    %    xr - state reference
    %    u  - input
    %    ur - input reference
    %    p  - parameter
    %    t  - time
    %    data- structured variable containing the values of additional data used inside
    %          the function
    %
    % Output:
    %    stageCost - Scalar or vectorized stage cost
    %
    %  Remark: If the stagecost does not depend on variables it is necessary to multiply
    %          the assigned value by t in order to have right vector dimesion when called for the optimization. 
    %          Example: stageCost = 0*t;
    stageCost = 0*t;
end

    
function boundaryCost = E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 
    % E_unscaled - Returns the boundary value cost
    %
    % Syntax:  boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 
    %
    % Inputs:
    %    x0  - state at t=0
    %    xf  - state at t=tf
    %    u0  - input at t=0
    %    uf  - input at t=tf
    %    p   - parameter
    %    tf  - final time
    %    data- structured variable containing the values of additional data used inside
    %          the function
    %
    % Output:
    %    boundaryCost - Scalar boundary cost
    %
    if strcmp(data.objective,"mf")
        boundaryCost=-xf(7);    % -(final mass)
    elseif strcmp(data.objective,"tof")
        boundaryCost=tf;        % time of flight
    end
end
    

function bc = b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
    % b_unscaled - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
    %
    % Syntax:  bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
    %
    % Inputs:
    %    x0  - state at t=0
    %    xf  - state at t=tf
    %    u0  - input at t=0
    %    uf  - input at t=tf
    %    p   - parameter
    %    tf  - final time
    %    data- structured variable containing the values of additional data used inside
    %          the function
    %
    %          
    % Output:
    %    bc - column vector containing the evaluation of the boundary function 
    %
    %------------- BEGIN CODE --------------
    varargin = varargin{1};
    bc = [];
    %bc=[xf(2)^2+xf(3)^2;xf(4)^2+xf(5)^2;xf(2)*xf(4)+xf(3)*xf(5);xf(3)*xf(4)-xf(5)*xf(2)];
    %------------- END OF CODE --------------
    % When adpative time interval add constraint on time
    %------------- BEGIN CODE --------------
    if length(varargin) == 2
        option = varargin{1};
        t_segment = varargin{2};
        if ((strcmp(options.discretization,'hpLGR')) || (strcmp(options.discretization,'globalLGR')))  && options.adaptseg==1 
            if size(t_segment,1)>size(t_segment,2)
                bc=[bc;diff(t_segment)];
            else
                bc=[bc,diff(t_segment)];
            end
        end
    end
end
    