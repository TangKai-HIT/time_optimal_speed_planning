function [dq_t, ddq_t, Tau, t_sample, exitFlag]=timeOptimSpeedPlan_toppra(q_s, dq_s, ddq_s, h, M_fn, C_fn, G_fn, constraints)
% TIMEOPTIMSPEEDPLAN_TOPPRA time optimal speed planning with kinematic and dynamic constraints 
% by reachability analysis. 
% *note: reparameterized constraints be like: d_i*a + c_i*b + g_i (or a_i*u + b_i*x + c_i, a=u=dds_ddt; b=x=ds_dt^2) 
%        b_{i+1} = b_i + 2*h*a_i (x_{i+1} = x_i + 2*h*u_i)

%   Inputs:
%       q_s, dq_s, ddq_s: N X dim, sampled trajectory on parameterized curves (e.g. cubic spline)
%       h: sample parameter gap--s_f/(N-1)
%       M_fn, C_fn, G_fn: function handles for evaluating mass term M,
%                                       velocity product C(q, dq)*dq and gravity term G respectively (with row vector input)
%       constraints: constraints struct
%   Outputs:
%       dq_t: sampled velocity w.r.t time (N X dim)
%       ddq_t: sampled  acceleration w.r.t time (N-1 X dim)
%       Tau: sampled generalized force term (N-1 X dim)
%       t_sample: corresponding time samples
%       exitFlag: 1--success, 0--failed

%% Init variables
[N, dim] = size(q_s);

Tau = zeros(N-1, dim);
t_sample = zeros(1, N);

%% Get linear inequalities
%init
A = zeros(2*2*dim, 2, N-1); %2*dynamic*acceleration
b = zeros(2*2*dim, 1, N-1);

Ub = ones(2, N-1) * inf; %[x; u]
Lb = - ones(2, N-1) * inf;
Lb(1, :) = zeros(1, N-1); % x=ds_dt^2>=0

% dynamic (-torque <= d*u + c*x + g <= torque)
D_d = zeros(N-1, dim); %matrix of d_{i,j}
D_c =zeros(N-1, dim); %matrix of c_{i,j}
D_g = zeros(N-1, dim); %matrix of g_{i,j}

% kinematic
D_torque = constraints.Mu'; %dim X 1, torque constraints 
D_accel = constraints.Alpha'; %dim X 1, acceleration constraints

% generate inequality matrices & upper bound
 for i=1:N-1
    %dynamic inequality constraints (-torque <= d*u + c*x + g <= torque)
    D_d(i,:) = M_fn(q_s(i, :)) * dq_s(i, :)';
    D_c(i,:) = M_fn(q_s(i, :)) * ddq_s(i, :)' + C_fn(q_s(i, :), dq_s(i, :))';
    D_g(i,:) = G_fn(q_s(i, :));
    
    A(1:dim, :, i) = [D_c(i,:)', D_d(i,:)'];
    b(1:dim, :, i) = D_torque - D_g(i,:)';

    A(dim+1 : 2*dim, :, i) = [-D_c(i,:)', -D_d(i,:)'];
    b(dim+1 : 2*dim, :, i) = D_torque + D_g(i,:)';

    %acceleration inequality constraints (-alpha <= dq_ds*u + ddq_dds*x <= alpha)
    A(2*dim+1 : 3*dim, :, i) = [ddq_s(i, :)', dq_s(i, :)']; %[x, u]
    b(2*dim+1 : 3*dim, :, i) = D_accel;
    
    A(3*dim+1 : 4*dim, :, i) = [-ddq_s(i, :)', -dq_s(i, :)']; %[x, u]
    b(3*dim+1 : 4*dim, :, i) = D_accel;

    %set ub
    Ub(1, i) = min((constraints.Phi.^2) ./ (dq_s(i, :).^2)); %speed constraints
 end

%% Solve toppra
exitFlag = 1;
disp("Doing Speed Planning (solving toppra)......");

tic;

%backward
[controllable_sets, backward_exitFlag] = toppra_backward_pass(A, b, Lb, Ub, constraints.b_end, h, N);
if ~backward_exitFlag
    exitFlag = 0;
    return;
end

%forward
if constraints.b_init<controllable_sets(1,1) || constraints.b_init>controllable_sets(1,2)
    disp("Initial condition not in controllable set! Forward passed failed!");
    exitFlag = 0;
    return;
end

[b_optim, a_optim, forward_exitFlag] = toppra_forward_pass(A, b, Lb(2, :), Ub(2, :), controllable_sets, constraints.b_init, h, N);
if ~forward_exitFlag
    exitFlag = 0;
    return;
end

solveTime = toc;

fprintf("Finished! Run Time:%.4f\n", solveTime);

%% Compute outputs
v_optim = sqrt(b_optim);
dq_t = v_optim .* dq_s;
ddq_t = a_optim .* dq_s(1:end-1, :) + b_optim(1:end-1) .* ddq_s(1:end-1, :);

for i=1:N-1
    t_sample(i+1) = t_sample(i) + 1/(v_optim(i) + v_optim(i+1));
    Tau(i, :) = D_d(i, :) * a_optim(i) + D_c(i, :) * b_optim(i) + D_g(i, :);
end
t_sample = t_sample .* (2*h);
