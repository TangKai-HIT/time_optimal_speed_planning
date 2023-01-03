function [dq_t, ddq_t, Tau, t_min]=timeOptimSpeedPlan(q_s, dq_s, ddq_s, h, M_fn, C_fn, G_fn, constraints)
% TIMEOPTIMSPEEDPLAN time optimal speed planning with kinematic and dynamic constraints 
% by solving a discretized convex optimization problem
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
%       t_min: optimal time

%% Init variables
[N, dim] = size(q_s);

Tau = zeros(N-1, dim);
t_sample = zeros(1, N);

b_optim = zeros(N, 1);
b_optim(1) = constraints.b_init;
b_optim(end) = constraints.b_end;

%% Set linear inequality constraints
% dynamic
D_d = zeros(dim*(N-1), N-1);
D_c = zeros(dim*(N-1), N-1);
g = zeros(dim*(N-1), 1);

% kinematic
D_dgamma = zeros(dim*(N-1), N-1);
D_ddgamma = zeros(dim*(N-1), N-1);

Mu = kron(ones(N-1, 1), constraints.Mu');
Alpha = kron(ones(N-1, 1), constraints.Alpha');
Lb = zeros(N-2, 1);
Ub = zeros(N-2, 1);

% generate inequality matrices & upper bound
[K, e]=getDiffMatrices(N-2, 1, constraints.b_init, constraints.b_end);

 for i=1:N-1
    D_d((i-1)*dim+1 : i*dim, i) = M_fn(q_s(i, :)) * dq_s(i, :)';
    D_c((i-1)*dim+1 : i*dim, i) = M_fn(q_s(i, :)) * ddq_s(i, :)' + C_fn(q_s(i, :), dq_s(i, :))';
    g((i-1)*dim+1 : i*dim, 1) = G_fn(q_s(i, :));

    D_dgamma((i-1)*dim+1 : i*dim, i) = dq_s(i, :);
    D_ddgamma((i-1)*dim+1 : i*dim, i) = ddq_s(i, :);
    %set ub
    if i>1
        Ub(i-1) = min((constraints.Phi.^2) ./ (dq_s(i, :).^2));
    end
 end

 A1_t = D_d * K ./(2*h) + D_c(:, 2:end);
 C1_t = D_d * e ./(2*h) + g + constraints.b_init *  D_c(:, 1);

 A2_t = D_dgamma * K ./(2*h) + D_ddgamma(:, 2:end);
 C2_t = D_dgamma * e ./(2*h) + constraints.b_init *  D_ddgamma(:, 1);

 A = [A1_t; -A1_t; A2_t; -A2_t];
 b = [Mu - C1_t; 
        Mu + C1_t; 
        Alpha - C2_t;
        Alpha + C2_t];

%% Solve LP
f = -ones(N-2,1);
options = optimoptions(@linprog);

disp("Doing Speed Planning (solving LP using dual-simplex)......");
tic;
b_optim(2 : N-1) = linprog(f, A, b, [], [], Lb, Ub, options);
solveTime = toc;
fprintf("Finished! Run Time:%.4f\n", solveTime);

%% Compute outputs
a_optim = diff(b_optim) ./ (2*h);
v_optim = sqrt(b_optim);
dq_t = v_optim .* dq_s;
ddq_t = a_optim .* dq_s(1:end-1, :) + b_optim(1:end-1) .* ddq_s(1:end-1, :);

t_min = 0;
for i=1:N-1
    if (v_optim(i) + v_optim(i+1)) > 0
        t_min = t_min + 1/(v_optim(i) + v_optim(i+1));
    end
    Tau(i, :) = D_d((i-1)*dim+1 : i*dim, i) * a_optim(i) + D_c((i-1)*dim+1 : i*dim, i) * b_optim(i) + g((i-1)*dim+1 : i*dim, 1);
end
t_min = t_min * (2*h);
