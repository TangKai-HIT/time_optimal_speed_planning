%Testing time optimal speed planning algorithm using a 3-DOF robot
clear; clc; close all;
addpath(genpath(pwd));

%% Init  robot
[robot,importInfo] = importrobot('robot_3R.slx','DataFormat','row');
T_tip = trvec2tform([0,0.325,0]);
% show(robot); showdetails(importInfo);
addCollision(robot.Bodies{1}, "cylinder", [0.02, 0.16]);
addCollision(robot.Bodies{2}, "cylinder", [0.02, 0.28], eul2tform([-pi/2,0,0],'XYZ')*trvec2tform([0,0,0.14]));
addCollision(robot.Bodies{3}, "cylinder", [0.02, 0.325], eul2tform([-pi/2,0,0],'XYZ')*trvec2tform([0,0,0.325/2]));
show(robot,'Collisions','on','Visuals','off');
%Waypoints & Boundaris
waypoints.s = [0,0.25,0.5,0.75,1];

waypoints.q = [   0,           0,             0;
                            1.288, -0.2864, -0.2982;
                            2.59,    -0.03045,  -0.5995;
                            4.374,  -0.04647,  -0.582;
                            5.334,  -0.1657,    -0.4504];
Phi = [2, 2, 2];
Alpha = [1.5, 1.5, 1.5];
Mu = [9, 9, 9];

%% 