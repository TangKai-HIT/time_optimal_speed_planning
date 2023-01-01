% A simple cubic spline example
addpath(genpath(pwd));
%% Init 
t_wp = [0, 5, 7, 8, 10, 15, 18];
q_wp = [3; -2; -5; 0; 6; 12; 8];
init_vel = 2;
end_vel = -3;
numSample = 100;

%% Compute trajectory
[q, dq, ddq, dddq, t_sample, splineCoeff]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, numSample);

%% Plot
figure('Position',[500,100, 1000,800]);
subplot(2,2,1)
plot(t_sample, q, '-r', 'LineWidth', 0.8); hold on;
scatter(t_wp, q_wp, 'black', 'filled' ,'o', SizeData=20);
title('pos');

subplot(2,2,2)
plot(t_sample, dq, '-r', 'LineWidth', 0.8); 
title('vel');

subplot(2,2,3)
plot(t_sample, ddq, '-r', 'LineWidth', 0.8); 
title('accel');

subplot(2,2,4)
plot(t_sample, dddq, '-r', 'LineWidth', 0.8); 
title('jerk');
