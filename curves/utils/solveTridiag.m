function x = solveTridiag(A, d)
% SOLVETRIDIAG solve linear equation with tridiagonal matrix A, with A being invertible
%   A: n X 3 matrix with tridiagonal elements stored in each row -- [a_n, b_n, c_n]
%   d: n X dim

%% Perform Thomas algorithm
N = size(A, 1);
dim = size(d, 2);
x = zeros(N, dim);

if abs(A(1, 2)) <1e-8 % if b_1==0
    x(2, :) = d(1, :)/A(1, 3);

    A_rec = A(2:end, :);
    b2 = A_rec(1, 2); A_rec(1, 2) = A_rec(1, 1);
    A_rec(1, 1) = 0;
    a3 = A_rec(2, 1); A_rec(2, 1) = 0;

    d_rec = d(2:end, :);
    d_rec(1,:) = d_rec(1,:) - b2 * x(2, :);
    d_rec(2,:) = d_rec(2,:) - a3 * x(2, :);

    x_rec = solveTridiag(A_rec, d_rec);
    x([1, 3:N], :) =  x_rec;
else
    %Forward Elimination
    for k=2:N
        m = A(k, 1)/A(k-1, 2);
        A(k, 2) = A(k, 2) - m * A(k-1, 3);
        d(k, :) = d(k, :) - m * d(k-1, :);
    end
    
    % Backward substitution
    x(end, :) = d(end, :)/A(end, 2);
    for k = N-1: -1: 1
        x(k, :) = (d(k, :) - A(k, 3)*x(k+1, :))/A(k, 2);
    end
end
