function [A, b, ub]=get_LPIneq(D_d, D_c, D_g, D_mu, D_dgamma, D_ddgamma, D_alpha, h, ub)
%GET_LPINEQ return linear inequalities constraints of the time optimal speed planning problem with special structure (LP).

dim = size(D_d, 2);
N = length(ub);

A = -ones(2*(N-1)*(2*dim), N);
b = -ones(2*(N-1)*(2*dim), 1);

index = 1;
for i=1:N-1
    for j=1:dim
        F_k_ij = -ones(1,2); %slope of f^{i,j}_k
        F_b_ij = -ones(1,2); %intersection term of f^{i,j}_k
        B_k_ij = -ones(1,2);  %slope of b^{i,j}_k
        B_b_ij = -ones(1,2); %intersection term of b^{i,j}_k

        %d_{i,j} c_{i,j}
        if D_d(i, j)>0 && D_c(i,j)>=0
            B_k_ij(1) = D_d(i, j)/(D_d(i, j)+2*h*D_c(i, j));
            B_b_ij(1) = 2*h*(D_mu(i, j) - D_g(i, j))/(D_d(i, j)+2*h*D_c(i, j));
            F_k_ij(1)= (D_d(i, j)+2*h*D_c(i, j))/D_d(i, j);
            F_b_ij(1) = 2*h*(D_mu(i, j) + D_g(i, j))/D_d(i, j);

        elseif D_d(i, j)<0 && D_c(i,j)<=0
            B_k_ij(1) = D_d(i, j)/(D_d(i, j) + 2*h*D_c(i, j));
            B_b_ij(1) = - 2*h*(D_mu(i, j) + D_g(i, j)) / (D_d(i, j)+2*h*D_c(i, j));
            F_k_ij(1) = (D_d(i, j)+2*h*D_c(i, j)) / D_d(i, j);
            F_b_ij(1) = - 2*h*(D_mu(i, j) - D_g(i, j)) / D_d(i, j);
        
        elseif D_d(i, j)>0 && D_c(i,j)<=0
            B_k_ij(1) = (D_d(i, j) - 2*h*D_c(i, j)) / D_d(i, j);
            B_b_ij(1) = 2*h*(D_mu(i, j) - D_g(i, j)) / D_d(i, j);
            F_k_ij(1) = D_d(i, j) / (D_d(i, j) - 2*h*D_c(i, j));
            F_b_ij(1) = 2*h*(D_mu(i, j) + D_g(i, j)) / (D_d(i, j) - 2*h*D_c(i, j));
        
        elseif D_d(i, j)<0 && D_c(i,j)>=0
            B_k_ij(1) = (D_d(i, j) - 2*h*D_c(i, j)) / D_d(i, j);
            B_b_ij(1) = - 2*h*(D_mu(i, j) + D_g(i, j)) / D_d(i, j);
            F_k_ij(1) = D_d(i, j) / (D_d(i, j) - 2*h*D_c(i, j));
            F_b_ij(1) = - 2*h*(D_mu(i, j) - D_g(i, j)) / (D_d(i, j) - 2*h*D_c(i, j));

        elseif D_d(i, j)==0 && D_c(i,j)~=0
            ub(i+1) = min(ub(i+1), (D_mu(i,j) - D_g(i, j))/abs(D_c(i, j)));
        end

        %gamma_{i,j} dgamma_{i,j}
        if D_dgamma(i, j) *  D_ddgamma(i, j) >0 || (D_dgamma(i, j)~=0 && D_ddgamma(i, j)==0)
            B_k_ij(2) = D_dgamma(i, j)/(D_dgamma(i, j) + 2*h*D_ddgamma(i, j));
            B_b_ij(2) = abs(2*h*D_alpha(i, j)/(D_dgamma(i, j) + 2*h*D_ddgamma(i, j)));
            F_k_ij(2) = (D_dgamma(i, j) + 2*h*D_ddgamma(i, j))/D_dgamma(i, j);
            F_b_ij(2) = abs(2*h*D_alpha(i, j)/D_dgamma(i, j));

        elseif D_dgamma(i, j) *  D_ddgamma(i, j) <0
            B_k_ij(2) = (D_dgamma(i, j) - 2*h*D_ddgamma(i, j)) / D_dgamma(i, j);
            B_b_ij(2) = abs(2*h*D_alpha(i, j)/D_dgamma(i, j));
            F_k_ij(2) =  D_dgamma(i, j) / (D_dgamma(i, j) - 2*h*D_ddgamma(i, j));
            F_b_ij(2) = abs(2*h*D_alpha(i, j)/(D_dgamma(i, j) - 2*h*D_ddgamma(i, j)));

        elseif D_dgamma(i, j)==0 && D_ddgamma(i, j)~=0
            ub(i+1) = min(ub(i+1), abs(D_alpha(i, j)/D_ddgamma(i, j)));
        end
        
        %remove -1s
        B_k_ij = B_k_ij(B_k_ij>=0);
        B_b_ij = B_b_ij(B_b_ij>=0);
        F_k_ij = F_k_ij(F_k_ij>=0);
        F_b_ij = F_b_ij(F_b_ij>=0);

        % Add (b_{i+1}<= B_k_ij*b_i+B_b_ij) to A, b
        if ~isempty(B_k_ij)
            for k=1:length(B_k_ij)
                A_newRow = zeros(1, N);
                A_newRow(i : i+1) = [-B_k_ij(k), 1];
                A(index, :) = A_newRow;
    
                b(index) = B_b_ij(k);
                index = index +1; %next row
            end
        end

        % Add (b_{i}<= F_k_ij*b_{i+1}+F_b_ij) to A, b
        if ~isempty(F_k_ij)
            for k=1:length(F_k_ij)
                A_newRow = zeros(1, N);
                A_newRow(i : i+1) = [1, -F_k_ij(k)];
                A(index, :) = A_newRow;
    
                b(index) = F_b_ij(k);
                index = index +1; %next row
            end
        end

    end
end

A = A(1:index-1, :);
b = b(1:index-1);
