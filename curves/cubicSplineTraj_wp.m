function [q, dq, ddq, dddq, t_sample, splineCoeff]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, numSample)
% CUBICSPLINETRAJ_WP plan cubic spline trajectory through way points with specified initial and final velocities, 
% intermediate way point velocities are computed by assuming continuous acceleration.
%   Inputs:
%       q_wp: waypoints, n X dim
%       t_wp: time at waypoints 1 X n
%       init_vel: 1 X dim
%       end_vel: 1 X dim
%       numSample: number of samples
%   Outputs:
%       t_sample: 1 X numSample
%       splineCoeff: 1 X 4 cell for spline coefficients

N = length(t_wp);
dim = size(q_wp, 2);
t_sample = linspace(t_wp(1), t_wp(end), numSample);

q = zeros(numSample, dim);
dq = zeros(numSample, dim);
ddq = zeros(numSample, dim);
dddq = zeros(numSample, dim);

[splineCoeff, ~]=getCubicSpline_wp(q_wp, t_wp, init_vel, end_vel);

for k=1:N-1
    id = find(t_sample>=t_wp(k) & t_sample<=t_wp(k+1));
    
    time_k = t_sample(id)';
    a0_k = splineCoeff{1}(k, :);
    a1_k = splineCoeff{2}(k, :);
    a2_k = splineCoeff{3}(k, :);
    a3_k = splineCoeff{4}(k, :);

    q(id, :) =  a0_k + (time_k - t_wp(k)) .* a1_k + (time_k - t_wp(k)).^2 .* a2_k + (time_k - t_wp(k)).^3 .* a3_k;
    dq(id, :) =  a1_k + 2*(time_k - t_wp(k)) .* a2_k + 3*(time_k - t_wp(k)).^2 .* a3_k;
    ddq(id, :) =  2*a2_k + 6*(time_k - t_wp(k)) .* a3_k;
    dddq(id, :) =  6 * a3_k;
end