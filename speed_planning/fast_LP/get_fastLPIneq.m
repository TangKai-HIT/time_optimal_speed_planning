function [slopes, intersections]=get_fastLPIneq(D_d, D_c, D_g, D_mu, D_dgamma, D_ddgamma, D_alpha, h, ub)
%GET_FASTLPINEQ return slops and intersections of forward &backward
%linear inequalities of the time optimal speed planning problem with special structure.
%   Output:
%       slopes:  1X2 cell -- slops{1}->b_k(x_i) (decreasing order in dim);  
%                                        slops{2}->f_k(x_{i+1}) (decreasing order in dim, then 1/f_k is in increasing order)
%       intersections: 1X2 cell -- intersections{1}->b_k(x_i) (correspond to elements in slops);  
%                                                     intersections{2}->f_k(x_{i+1}) (correspond to elements in slops)
%   *Note that -1 elements in slops & intersections are preserved (useless) positions

[n, dim] = size(D_d); %n=N-1

slopes = {-ones(n, 2*dim+1), -ones(n, 2*dim)};
intersections = {-ones(n, 2*dim+1), -ones(n, 2*dim)};

for i=1:n
    F_ki = -ones(1, 2*dim); %slope of f^i_k
    F_bi = -ones(1, 2*dim); %intersection of f^i_k
    B_ki = -ones(1, 2*dim+1);  %slope of b^i_k
    B_bi = -ones(1, 2*dim+1); %intersection of b^i_k

    for j=1:dim
        %d_{i,j} c_{i,j}
        if D_d(i, j)>0 && D_c(i,j)>=0
            B_ki(1, j) = D_d(i, j)/(D_d(i, j)+2*h*D_c(i, j));
            B_bi(1, j) = 2*h*(D_mu(i, j) - D_g(i, j))/(D_d(i, j)+2*h*D_c(i, j));
            F_ki(1, j) = (D_d(i, j)+2*h*D_c(i, j))/D_d(i, j);
            F_bi(1, j) = 2*h*(D_mu(i, j) + D_g(i, j))/D_d(i, j);

        elseif D_d(i, j)<0 && D_c(i,j)<=0
            B_ki(1, j) = D_d(i, j)/(D_d(i, j) + 2*h*D_c(i, j));
            B_bi(1, j) = - 2*h*(D_mu(i, j) + D_g(i, j)) / (D_d(i, j)+2*h*D_c(i, j));
            F_ki(1, j) = (D_d(i, j)+2*h*D_c(i, j)) / D_d(i, j);
            F_bi(1, j) = - 2*h*(D_mu(i, j) - D_g(i, j)) / D_d(i, j);
        
        elseif D_d(i, j)>0 && D_c(i,j)<=0
            B_ki(1, j) = (D_d(i, j) - 2*h*D_c(i, j)) / D_d(i, j);
            B_bi(1, j) = 2*h*(D_mu(i, j) - D_g(i, j)) / D_d(i, j);
            F_ki(1, j) = D_d(i, j) / (D_d(i, j) - 2*h*D_c(i, j));
            F_bi(1, j) = 2*h*(D_mu(i, j) + D_g(i, j)) / (D_d(i, j) - 2*h*D_c(i, j));
        
        elseif D_d(i, j)<0 && D_c(i,j)>=0
            B_ki(1, j) = (D_d(i, j) - 2*h*D_c(i, j)) / D_d(i, j);
            B_bi(1, j) = - 2*h*(D_mu(i, j) + D_g(i, j)) / D_d(i, j);
            F_ki(1, j) = D_d(i, j) / (D_d(i, j) - 2*h*D_c(i, j));
            F_bi(1, j) = - 2*h*(D_mu(i, j) - D_g(i, j)) / (D_d(i, j) - 2*h*D_c(i, j));

        elseif D_d(i, j)==0 && D_c(i,j)~=0
            B_ki(1, j) = 0;
            B_bi(1, j) = (D_mu(i,j) - D_g(i, j))/abs(D_c(i, j));
        end

        %gamma_{i,j} dgamma_{i,j}
        if D_dgamma(i, j) *  D_ddgamma(i, j) >0 || (D_dgamma(i, j)~=0 && D_ddgamma(i, j)==0)
            B_ki(1, j+dim) = D_dgamma(i, j)/(D_dgamma(i, j) + 2*h*D_ddgamma(i, j));
            B_bi(1, j+dim) = abs(2*h*D_alpha(i, j)/(D_dgamma(i, j) + 2*h*D_ddgamma(i, j)));
            F_ki(1, j+dim) = (D_dgamma(i, j) + 2*h*D_ddgamma(i, j))/D_dgamma(i, j);
            F_bi(1, j+dim) = abs(2*h*D_alpha(i, j)/D_dgamma(i, j));

        elseif D_dgamma(i, j) *  D_ddgamma(i, j) <0
            B_ki(1, j+dim) = (D_dgamma(i, j) - 2*h*D_ddgamma(i, j)) / D_dgamma(i, j);
            B_bi(1, j+dim) = abs(2*h*D_alpha(i, j)/D_dgamma(i, j));
            F_ki(1, j+dim) =  D_dgamma(i, j) / (D_dgamma(i, j) - 2*h*D_ddgamma(i, j));
            F_bi(1, j+dim) = abs(2*h*D_alpha(i, j)/(D_dgamma(i, j) - 2*h*D_ddgamma(i, j)));

        elseif D_dgamma(i, j)==0 && D_ddgamma(i, j)~=0
            B_ki(1, j+dim) = 0;
            B_bi(1, j+dim) = abs(D_alpha(i, j)/D_ddgamma(i, j));
        end
    end

    %u_{i+1}
    B_ki(1, end) = 0;
    B_bi(1, end) = ub(i+1);

    %b_k (decreasing)
    [slopes{1}(i, :), index] = sort(B_ki, 'descend');
    intersections{1}(i, :) = B_bi(index);
    %f_k (decreasing)
    [slopes{2}(i, :), index] = sort(F_ki, 'descend');
    intersections{2}(i, :) = F_bi(index);
end
