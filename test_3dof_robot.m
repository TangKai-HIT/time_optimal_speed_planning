%Testing time optimal speed planning algorithm using a 3-DOF robot
clear; clc; close all;
addpath(genpath(pwd));

%% Init  robot & parameters
[robot,importInfo] = importrobot('robot_3R.slx','DataFormat','row');
T_tip = trvec2tform([0,0.325,0]);
% show(robot); showdetails(importInfo);
addCollision(robot.Bodies{1}, "cylinder", [0.02, 0.16]);
addCollision(robot.Bodies{2}, "cylinder", [0.02, 0.28], eul2tform([-pi/2,0,0],'XYZ')*trvec2tform([0,0,0.14]));
addCollision(robot.Bodies{3}, "cylinder", [0.02, 0.325], eul2tform([-pi/2,0,0],'XYZ')*trvec2tform([0,0,0.325/2]));
% show(robot,'Collisions','on','Visuals','off');
%Waypoints & Boundaris
s_wp = [0,0.25,0.5,0.75,1];

q_wp = [   0,           0,             0;
                            1.288, -0.2864, -0.2982;
                            2.59,    -0.03045,  -0.5995;
                            4.374,  -0.04647,  -0.582;
                            5.334,  -0.1657,    -0.4504];
Phi = [2, 2, 2]; %velocity constraint
Alpha = [1.5, 1.5, 1.5]; %acceleration constraint
Mu = [9, 9, 9]; %torque constraint

%% Interpolate using cubic spline
init_vel = [0,0,0];
end_vel = [0,0,0];
numSample = 100;
h = 1 / (numSample - 1);
[q_s, dq_s, ddq_s, dddq_s, s_sample, splineCoeff]=cubicSplineTraj_wp(q_wp, s_wp, init_vel, end_vel, numSample);

%% Speed planning
constraints = getSpeedPlanConstraints(Phi, Alpha, Mu, 0, 0);
M_fn = @(q) massMatrix(robot, q);
C_fn = @(q, dq) velocityProduct(robot, q, dq);
G_fn = @(q) gravityTorque(robot, q);

%Use nonlinear convex optimization: 
% options = optimoptions(@fmincon, "MaxFunctionEvaluations", 5000, "PlotFcn", {'optimplotfvalconstr'});
% options = optimoptions(@fmincon, "MaxFunctionEvaluations", 5000);
% [dq_t, ddq_t, Tau, t_sample]=timeOptimSpeedPlan(q_s, dq_s, ddq_s, h, M_fn, C_fn, G_fn, constraints, options);

%Use LP
% [dq_t, ddq_t, Tau, t_sample]=timeOptimSpeedPlan_LP(q_s, dq_s, ddq_s, h, M_fn, C_fn, G_fn, constraints);

%Use fast LP
[dq_t, ddq_t, Tau, t_sample]=timeOptimSpeedPlan_fastLP(q_s, dq_s, ddq_s, h, M_fn, C_fn, G_fn, constraints);

fprintf("Optimal Time: %0.3f\n", t_sample(end));
% plot results
plotSpeedPlanResult(dq_t, ddq_t, Tau, t_sample, constraints);

%% Get robot motion joint trajectory PD control simulation results
%define joint Space Motion Model
jointMotionModel = jointSpaceMotionModel('RigidBodyTree',robot,'MotionType','PDControl');

%numerical simulation of the dynamic mdel using ode15s
timeInterval = [0, t_sample(end)];
jointStateTraj = [q_s'; dq_t'];
[t_sim, jointState] = ode15s(@(t, state) jointTrajTrackModel(jointMotionModel, state, jointStateTraj, t, t_sample), ...
                                       timeInterval, jointStateTraj(:,1));

%% Plot motion
figure();
show(robot, q_s(1, :), 'Collisions','on', 'PreservePlot', false);
hold on;
axis([-0.5 0.5 -0.5 0.5 -0.1 0.5]);

for i=1:numSample
    % Current time 
    tNow= t_sample(i);
    % Interpolate simulated joint positions to get configuration at current time
    configNow = interp1(t_sim, jointState(:,1: 3), tNow);
    poseNow = getTransform(robot, configNow, "Body3");
    show(robot,configNow, 'Collisions','on', 'PreservePlot',false);
    jointSpaceMarker = plot3(poseNow(1,4),poseNow(2,4),poseNow(3,4),'r.','MarkerSize',20);
    drawnow;
end

% Add a legend and title
legend(jointSpaceMarker, 'endeffector trajectory');
title('Manipulator Trajectories')