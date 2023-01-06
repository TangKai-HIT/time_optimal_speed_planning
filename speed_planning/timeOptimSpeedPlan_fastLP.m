function [dq_t, ddq_t, Tau, t_sample]=timeOptimSpeedPlan_fastLP(q_s, dq_s, ddq_s, h, M_fn, C_fn, G_fn, constraints)
% TIMEOPTIMSPEEDPLAN_FASTLP time optimal speed planning with kinematic and dynamic constraints 
% by solving fast LP. The problem should satisfy several assumptions to form the structure: 
%   (1) max torque>= g term;

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

%% Init variables
[N, dim] = size(q_s);

Tau = zeros(N-1, dim);
t_sample = zeros(1, N);

%% Get fast LP linear inequalities
% dynamic
D_d = zeros(N-1, dim); %matrix of d_{i,j}
D_c =zeros(N-1, dim); %matrix of c_{i,j}
D_g = zeros(N-1, dim); %matrix of g_{i,j}

% kinematic
D_dgamma = zeros(N-1, dim); %matrix of dgamma_{i,j}
D_ddgamma = zeros(N-1, dim); %matrix of ddgamma_{i,j}

D_mu = kron(ones(N-1, 1), constraints.Mu);
D_alpha = kron(ones(N-1, 1), constraints.Alpha);
Ub = zeros(N, 1);

% generate inequality matrices & upper bound
 for i=1:N
    if i<N
        D_d(i, :) = M_fn(q_s(i, :)) * dq_s(i, :)';
        D_c(i, :) = M_fn(q_s(i, :)) * ddq_s(i, :)' + C_fn(q_s(i, :), dq_s(i, :))';
        D_g(i, :) = G_fn(q_s(i, :));
    
        D_dgamma(i, :) = dq_s(i, :);
        D_ddgamma(i, :) = ddq_s(i, :);
    end
    %set ub
    Ub(i) = min((constraints.Phi.^2) ./ (dq_s(i, :).^2));
 end

 %get slopes and intersections for fast LP
 [slopes, intersections]=get_fastLPIneq(D_d, D_c, D_g, D_mu, D_dgamma, D_ddgamma, D_alpha, h, Ub);

%% Solve fast LP
disp("Doing Speed Planning (solving fast LP)......");
tic;
b_optim = fastLP_solve(Ub, slopes, intersections, constraints.b_init, constraints.b_end);
solveTime = toc;
fprintf("Finished! Run Time:%.4f\n", solveTime);

%% Compute outputs
a_optim = diff(b_optim) ./ (2*h);
v_optim = sqrt(b_optim);
dq_t = v_optim .* dq_s;
ddq_t = a_optim .* dq_s(1:end-1, :) + b_optim(1:end-1) .* ddq_s(1:end-1, :);

for i=1:N-1
    t_sample(i+1) = t_sample(i) + 1/(v_optim(i) + v_optim(i+1));
    Tau(i, :) = D_d(i, :) * a_optim(i) + D_c(i, :) * b_optim(i) + D_g(i, :);
end
t_sample = t_sample .* (2*h);
