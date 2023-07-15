function [x_optim, u_optim, exitFlag] = toppra_forward_pass(A, b, lb, ub, K_set, x_init, h, N)
%FORWARD_PASS forward pass phase of TOPPRA
%   (x_i = dot_s_i^2, u_i = ddot_s_i)
%   Inputs: 
%       A: k X 2 X N-1, inequality LHS(left hand side) set of [x_i, u_i]
%       b: k X 1 X N-1, inequality RHS set of [x_i, u_i]
%       lb: 1 X N-1 or N-1 X 1, lower bound of u_i
%       ub: 1 X N-1 or N-1 X 1, upper bound of u_i
%       K_set: N X 2, controllable set K (allowed range of x_i = dot_s_i^2)
%       h: scalar, step -- s_{i+1} - s_i
%       x_init: init condition of x, x_0
%   Outputs:
%       x_optim, u_optim
%       exitFlag: 1--success, 0--no solution

k = size(A, 1);
% x_end = K_set(end, 1);

x_optim = zeros(N, 1);
x_optim(1) = x_init;

u_optim = zeros(N-1, 1);

%greedily select the controls u_i
exitFlag = 1;

for i=1:N-1
    lb_i = lb(i);
    ub_i = ub(i);
    x_i = x_optim(i);

    %Ax* <= b
     for j = 1:k
         %mind singularity point may exist: A(j,1)*x + A(j,2)*u <= b, A(j,2)==0
         if abs(A(j, 2, i)) > 1e-8
             if A(j, 2, i) > 0 
                ub_i = min(ub_i, (b(j, 1, i) - A(j, 1, i)*x_i) / A(j, 2, i));
             else
                lb_i = max(lb_i, (b(j, 1, i) - A(j, 1, i)*x_i) / A(j, 2, i));
             end
         end
     end

    %x*_i + 2h*u_i <=x*_{i+1}_max
    ub_i = min(ub_i, (K_set(i+1, 2) - x_i) / (2*h));

    %x*_i + 2h*u_i >= x*_{i+1}_min
    lb_i = max(lb_i, (K_set(i+1, 1) - x_i) / (2*h));
    
    %update u*_i
    if lb_i < ub_i || abs(ub_i - lb_i)<1e-8
        u_optim(i) = ub_i;
    else
        fprintf("Forward pass failed! Infeasible range in step %d\n", i);
        exitFlag = 0;
        return;
    end

    %update x*_{i+1}
    x_optim(i+1) = x_optim(i) + 2*h*u_optim(i);
end

end

