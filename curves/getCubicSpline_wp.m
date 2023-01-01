function splineParams=getCubicSpline_wp(q_wp, t_wp, init_vel, end_vel)
% GETCUBICSPLINE_WP get cubic spline parameters for way point interpolation with specified initial and final velocities, 
% intermediate way point velocities are computed by assuming continuous acceleration.
%   q_wp: waypoints, n X dim
%   t_wp: time at waypoints 1 X n
%   init_vel: 1 X dim
%   end_vel: 1 X dim
%   splineParams: parameter struct

dim = size(q_wp, 2);
N = length(t_wp);

T = diff(t_wp);
%% Solve intermediate velocities
tridiag_v = zeros(N-2, 3); %store tridiagonal elements of the sparse tridiagonal matrix(N-2 X N-2)
c = zeros(N-2, dim); 

c(1, :) = -T(2)*init_vel;
c(end, :) = -T(end-1)*end_vel;

tridiag_v(1, :) = [0, 2*(T(1)+T(2)), T(1)];
tridiag_v(end, :) = [T(end), 2*(T(end-1)+T(end)), 0];

for n=1:N-2
    if n>1 && n<N-2
        tridiag_v(n, 1) = T(n+1);
        tridiag_v(n, 2) = 2*(T(n)+T(n+1));
        tridiag_v(n, 3) = T(n);
    end

    c(n, :) = c(n, :) + 3/(T(n)*T(n+1)) * (T(n)^2 * (q_wp(n+2,:) - q_wp(n+1,:)) + T(n+1)^2 * (q_wp(n+1,:) - q_wp(n,:)));
end


