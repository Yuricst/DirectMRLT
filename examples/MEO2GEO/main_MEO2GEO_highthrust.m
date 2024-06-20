%% Main script for orbit transfer problem in MEE
% Yuri Shimane, yuri.shimane@gatech.edu
%   Created    : 2024/06/19
%   Last edits : 2024/06/20

% house keeping
clear; close all; clc;
fontsize = 14;

fPath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(fullfile(fPath,'../../DirectMRLT/'));

% define data
GM_SUN = 398600.44;
LU = 42164.0;
MU = 1500;
thrust = 5.0;                  % in N
mdot = thrust/(9.81 * 1500);    % in kg/s
data = get_problem_data(GM_SUN,LU,MU,thrust,mdot);

% initial and final conditions
KEP_0 = [0.5 0.0 deg2rad(23) deg2rad(100) deg2rad(270) deg2rad(23)];
KEP_F = [1.0 0.0 deg2rad(0)  deg2rad(102) deg2rad(34)  deg2rad(5)];
MEE_0 = KEP2MEE(KEP_0);
MEE_F = KEP2MEE(KEP_F);
m0 = 1.0;                           % initial mass, in MU
t0 = 0.0;                           % initial time, in TU
tf_bounds = [5*pi 25*pi];          % bounds on TOF, in TU
mf_bounds = [0.7 m0];               % bounds on final mass, in MU

% initial and final orbit
RV_initial = MEE2RVorbit(data.GM,MEE_0);
RV_final   = MEE2RVorbit(data.GM,MEE_F);

% create problem
objective = "tof";     % "mf" for mass-optimal or "tof" for time-optimal
[problem,guess] = MEEOrbitTransferProblem(...
    data,MEE_0,MEE_F,m0,t0,tf_bounds,mf_bounds,@ICLOCSsettings_GTO2GEO_highthrust, ...
    "objective", objective,"max_rev",25);
RV_guess = MEE2RV(data.GM, guess.states(:,1:6));

% plot initial guess
figure('Position',[600,10,600,500]);
plot3(RV_initial(:,1),RV_initial(:,2),RV_initial(:,3),'-g','LineWidth',1.2);
hold on;
plot3(RV_final(:,1),RV_final(:,2),RV_final(:,3),'-r','LineWidth',1.2);
plot3(RV_guess(:,1),RV_guess(:,2),RV_guess(:,3),'-k','LineWidth',1.2);
xlabel("x, LU");
ylabel("y, LU");
zlabel("z, LU");
grid on; box on; axis equal;
set(gca,'fontsize',fontsize);


% solve problem
options= problem.settings(150);                  % h method
%options= problem.settings(100,4);                  % hp method
[solution,MRHistory] = solveMyProblem(problem, guess, options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.1 );

%% Plots
% plot initial guess
figure('Position',[100,10,1400,700]);
tiledlayout(2,4);
for i = 1:7
    nexttile([1,1]);
    plot(tv, xv(:,i),'-k','LineWidth',1.2);
    xlabel("Time, TU");
    ylabel(strcat(problem.state_names(i)));
    grid on; box on;
    set(gca,'fontsize',fontsize);
end
nexttile;
hold on;
for iu = 1:4
    % plot(tv, uv(:,iu),'LineWidth',1.2);
    plot(solution.T, solution.U(:,iu),'LineWidth',1.2);
end
legend("u_R","u_T","u_N","||u||");
xlabel("Time, TU");
ylabel('u');
grid on; box on;
set(gca,'fontsize',fontsize);
saveas(gcf,fullfile(fPath,strcat("GTO2GEO_statehistory_",objective,".png")));

% solved transfer, initial and final orbit
RV = MEE2RV(data.GM, xv(:,1:6));

figure('Position',[600,10,600,500]);
plot3(RV_initial(:,1),RV_initial(:,2),RV_initial(:,3),'-g','LineWidth',1.2);
hold on;
plot3(RV_final(:,1),RV_final(:,2),RV_final(:,3),'-r','LineWidth',1.2);
plot3(RV(:,1),RV(:,2),RV(:,3),'-k','LineWidth',1.2);
xlabel("x, LU");
ylabel("y, LU");
zlabel("z, LU");
grid on; box on; axis equal;
set(gca,'fontsize',fontsize);
saveas(gcf,fullfile(fPath,strcat("GTO2GEO_trajectory_",objective,".png")));
