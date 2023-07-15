function [K_set, exitFlag] = toppra_backward_pass(A, b, lb, ub, x_end, h, N)
%BACKWARD_PASS backward pass phase of TOPPRA
%   (x_i = dot_s_i^2, u_i = ddot_s_i)
%   Inputs: 
%       A: k X 2 X N-1, inequality LHS(left hand side) set of [x_i, u_i]
%       b: k X 1 X N-1, inequality RHS set of [x_i, u_i]
%       lb: 2 X N-1, lower bound of [x_i, u_i]
%       ub: 2 X N-1, upper bound of [x_i, u_i]
%       h: scalar, step -- s_{i+1} - s_i
%       x_end: end condition of x, x_N
%   Outputs:
%       K_set: N X 2, controllable set K (allowed range of x_i = dot_s_i^2)
%       exitFlag: 1--success, 0--no solution

k = size(A, 1);
K_set = -ones(N, 2); %init all as -1
K_set(end, :) = ones(1, 2) * x_end;

A_temp = zeros(k+2, 2);
b_temp = zeros(k+2, 1);

exitFlag = 1;

options = optimoptions('linprog','Display','off');

for i = N-1 : -1 : 1
    if i < N-1
        A_temp(1:k, :) = A(:, :, i);
        b_temp(1:k, :) = b(:, :, i);

        %x_i+1_min <= x_i + 2*h*u_i
        A_temp(k+1, :) = [-1, -2*h];
        b_temp(k+1) = - K_set(i+1, 1);

        %x_i+1_max >= x_i + 2*h*u_i
        A_temp(k+2, :) = [1, 2*h];
        b_temp(k+2) = K_set(i+1, 2);

        [~, negative_max_x, max_exitFlag, ~] = linprog([-1; 0],A_temp,b_temp,[],[],lb(:, i),ub(:, i),options); %min -x

        [~, min_x, min_exitFlag, ~] = linprog([1; 0],A_temp,b_temp,[],[],lb(:, i),ub(:, i),options); %min x
    else
        Aeq = [1, 2*h];
        beq = x_end;

        [~, negative_max_x, max_exitFlag, ~] = linprog([-1; 0],A(:, :, i),b(:, :, i),Aeq,beq,lb(:, i),ub(:, i),options); %min -x

        [~, min_x, min_exitFlag, ~] = linprog([1; 0],A(:, :, i),b(:, :, i),Aeq,beq,lb(:, i),ub(:, i),options); %min x
    end
    
    if (max_exitFlag == -2 || max_exitFlag == 3 || max_exitFlag == -9)
        fprintf("Backward pass failed! Infeasible max x at step %d, LP exitFlag=%d\n", i, max_exitFlag);
        exitFlag = 0;
        return;
    end

    if (min_exitFlag == -2 || min_exitFlag == 3 || min_exitFlag == -9)
        fprintf("Backward pass failed! Infeasible min x at step %d, LP exitFlag=%d\n", i, min_exitFlag);
        exitFlag = 0;
        return;
    end
    
    if min_x > (-negative_max_x)
        fprintf("Backward pass failed! Empty controllable set at step %d\n", i);
        exitFlag = 0;
        return;
    end

    K_set(i, 1) = min_x;
    K_set(i, 2) = - negative_max_x;
end

end

