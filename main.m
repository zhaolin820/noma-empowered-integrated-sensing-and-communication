clc
clear all
close all


addpath('./function/');
% Parameters
cvx_solver mosek; % MOSEK solver is used to accelerate the optimization
rho_c = 10; % regularization parameter of communication rate 
rho_r = 1; % regularization parameter of effective sensing power 
epsilon = 0.4; % reduction factor of the penalty parameter in each iteration

para = para_init(); % please set the system parameters in the function para_init()
theta_degree = -90:90;

% Generate user location
[user_loc, user_angle, d_BU] = generate_user_location(para);
plot_location(para,user_loc);

% Path loss
path_loss = para.pathloss_direct(d_BU)';
path_loss = sqrt(10.^((-para.noise - path_loss)/10));

% Generate channels
[h] = generate_channel(para, path_loss);

% successive convex approximization algorithm
[W,R_curr,P,rank_diff_all, R_all, P_all] = SCA_algorithm(para, h, rho_c, rho_r, epsilon);

% calculate beampattern
[pattern] = beampattern(sum(W,3),para.N,theta_degree);


% plot convergence
figure;
subplot(2,1,1); hold on; box on;
obj_all = rho_c*R_all + rho_r*P_all;
n = 0:length(R_all)-1;
yyaxis left; plot(n, obj_all, '-o', 'LineWidth',1.5); ylabel("Objective Value");
yyaxis right; plot(n, rank_diff_all,'-s', 'LineWidth',1.5); ylabel("Penalty Term");
xlabel("Number of Outer Iterations");
title("Convergence of the proposed algorithm");
legend('Objective Value', 'Penalty Term');

% plot beampattern
subplot(2,1,2);
plot(theta_degree, pattern, '-b', LineWidth=1.5);
xlim([-90,90]);
xlabel("Azimuth Angle (degree)"); ylabel("Beampattern");
title(['Achieved beampattern when R = ' num2str(sum(min(R_curr,[],2))) ' bps/Hz']);
    
