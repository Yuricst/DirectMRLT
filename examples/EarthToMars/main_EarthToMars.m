%% Main script for orbit transfer problem in MEE
% Yuri Shimane, yuri.shimane@gatech.edu
%   Created    : 2024/06/19
%   Last edits : 2024/06/20

% house keeping
clear; close all; clc;

fPath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(fullfile(fPath,'../../DirectMRLT/'));

% define data
GM_SUN = 1.3271244004193938E+11;
LU = 149.6e6;
MU = 1500;
thrust = 0.8;
mdot = thrust/(9.81 * 3500);
data = get_problem_data(GM_SUN,LU,MU,thrust,mdot);

% initial and final conditions, in LU & non-dim & rad
MEE_0 = KEP2MEE([1.00000011 0.01671022 deg2rad(0.00005) ...
                 deg2rad(-11.26064) deg2rad(102.94719) deg2rad(100.46435)]);
MEE_F = KEP2MEE([1.52366231 0.09341233 deg2rad(1.85061) ...
                 deg2rad(49.57854) deg2rad(336.04084) deg2rad(355.45332)]);
m0 = 1.0;                           % initial mass, in MU
t0 = 0.0;                           % initial time, in TU
tf_bounds = [0.5*pi 3*pi];          % bounds on TOF, in TU
mf_bounds = [0.3 m0];               % bounds on final mass, in MU

% create problem
objective = "mf";     % "mf" for mass-optimal or "tof" for time-optimal
[problem,guess] = MEEOrbitTransferProblem(...
    data,MEE_0,MEE_F,m0,t0,tf_bounds,mf_bounds,@ICLOCSsettings_Earth2Mars, ...
    "objective", objective);

% solve problem
% options= problem.settings(150);                  % h method
options= problem.settings(100,4);                  % hp method
[solution,MRHistory] = solveMyProblem(problem, guess, options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.1 );

%% Plots
% plot initial guess
fontsize = 14;
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
saveas(gcf,fullfile(fPath,strcat("EarthToMars_statehistory_",objective,".png")));

% solved transfer, initial and final orbit
RV = MEE2RV(data.GM, xv(:,1:6));
RV_initial = MEE2RVorbit(data.GM,MEE_0);
RV_final   = MEE2RVorbit(data.GM,MEE_F);

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
saveas(gcf,fullfile(fPath,strcat("EarthToMars_trajectory_",objective,".png")));
