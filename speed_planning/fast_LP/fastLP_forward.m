function u_bar = fastLP_forward(ub, slopes, intersections, u_init, u_end)
%FASTLP_FORWARD perform forward phase of fast LP

N = length(ub);
u_bar = ub;
u_bar(1) = u_init;
u_bar(end) = u_end;

%define evaluation functions for b and f_inv & fixed point solver
b_eval = @(t, k, b) k*t+b;
f_inv_eval = @(t, k, b) t/k - b/k;
f_b_fixPt = @(k_f, k_b, b_f, b_b) (k_f*b_b+b_f)/(1-k_f*k_b); %f comp b fixed point solver

for i=1:N-1
    %load slopes & intersections of linear inequalities, remove -1s
    B_ki = slopes{1}(i, :);   B_ki=B_ki(B_ki>=0);
    F_ki = slopes{2}(i, :);   F_ki=F_ki(F_ki>=0);
    B_bi = intersections{1}(i, :);    B_bi=B_bi(B_bi>=0);
    F_bi = intersections{2}(i, :);    F_bi=F_bi(F_bi>=0);
    
    eps = length(B_ki);
    phi = length(F_ki);

    %do update loops
    x = u_bar(i);
    y = -inf;
    z = inf;

    while y<z && abs(y-z)>1e-8
        %update k, y
        b_i_now = b_eval(x, B_ki(eps), B_bi(eps));
        b_i_next = b_eval(x, B_ki(eps-1), B_bi(eps-1));
        while eps>1 && b_i_next<b_i_now
            eps = eps-1;
            b_i_now = b_i_next; 
            if eps >1
                b_i_next = b_eval(x, B_ki(eps-1), B_bi(eps-1)); 
            end
        end
        k = eps;
        y = b_i_now;
        
        %update j, z
        if phi > 0
            f_inv_i_now = f_inv_eval(x, F_ki(phi), F_bi(phi));
            f_inv_i_next = f_inv_eval(x, F_ki(phi-1), F_bi(phi-1));
            while phi>1 && f_inv_i_next>f_inv_i_now
                phi = phi-1;
                f_inv_i_now = f_inv_i_next; 
                if phi >1
                    f_inv_i_next = f_inv_eval(x, F_ki(phi-1), F_bi(phi-1));
                end
            end
            j = phi;
            z = f_inv_i_now;

            %get fixed point
            x = f_b_fixPt(F_ki(j), B_ki(k), F_bi(j), B_bi(k));
        else
            z = -inf;
        end
    end
    
    if i>1
        u_bar(i) = x;
    end

    if i<N-1
        u_bar(i+1) = y;
    end
end