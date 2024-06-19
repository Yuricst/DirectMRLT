% Main script for orbit transfer problem in MEE

% house keeping
clear; close all; clc;
addpath('../DirectMRLT/')

% dynamics
GM = 1.0;
MEE_0 = [1 0 0 0 0 0];
MEE_F = [1.52 0 0 0.04 0.02 0];
m0 = 1.0;
t0 = 0.0;
tf_bounds = [0.1*pi 5*pi];
mf_bounds = [0.3 m0];

data = get_problem_data();

% create problem
[problem,guess] = MEEOrbitTransferProblem(...
    data,MEE_0,MEE_F,m0,t0,tf_bounds,mf_bounds,@settings_MEE);
RV = MEE2RV(GM, guess.states(:,1:6));

% [dx,g_eq] = dynamics_MEE_internal(guess.states,guess.inputs,[],...
%     guess.time,data);

% solve problem
options= problem.settings(150);                  % h method
% options= problem.settings(100,4);                  % hp method
[solution,MRHistory] = solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.1 );


%% Plots
% plot initial guess
fontsize = 14;
figure('Position',[100,10,1000,600]);
tiledlayout(2,4);
for i = 1:7
    nexttile([1,1]);
    plot(guess.time, guess.states(:,i),'-k','LineWidth',1.2);
    xlabel("Time, TU");
    ylabel(strcat(problem.state_names(i)));
    grid on; box on;
    set(gca,'fontsize',fontsize);
end
nexttile;
plot(guess.time, guess.inputs(:,4),'-k','LineWidth',1.2);
xlabel("Time, TU");
ylabel('u');
grid on; box on;
set(gca,'fontsize',fontsize);

% initial and final orbit
RV_initial = MEE2RVorbit(GM,MEE_0);
RV_final   = MEE2RVorbit(GM,MEE_F);

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