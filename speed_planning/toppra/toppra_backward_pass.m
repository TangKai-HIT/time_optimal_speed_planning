function [K_set, exitFlag] = toppra_backward_pass(A, b, lb, ub, x_end, N)
%BACKWARD_PASS backward pass phase of TOPPRA
%   (x_i = dot_s_i^2, u_i = ddot_s_i)
%   Inputs: 
%       A: k X 2 X N-1, inequality LHS(left hand side) set of [x_i, u_i]
%       b: k X 1 X N-1, inequality RHS set of [x_i, u_i]
%       lb: 2 X N-1, lower bound of [x_i, u_i]
%       ub: 2 X N-1, upper bound of [x_i, u_i]
%       x_end: end condition of x, x_N
%   Outputs:
%       K_set: N X 2, controllable set K (allowed range of x_i = dot_s_i^2)
%       exitFlag: 1--success, 0--no solution

k = size(A, 1);
K_set = zeros(N, 2);
K_set(end, :) = ones(1, 2) * x_end;

A_temp = zeros(k+2, 2);
b_temp = zeros(k+2, 1);

for i=1:N-1
    if i < N-1
        A_temp(1:k, :) = A(:, :, i);
        linprog
    else

    end
end

end

