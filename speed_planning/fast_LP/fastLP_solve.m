function u_bar = fastLP_solve(ub, slopes, intersections, u_init, u_end)
%FASTLP_SOLVE   solve special structured problem of fast LP   

N = length(ub);
f_eval = @(t, k, b) t*k + b;

%% Perform forward phase
u_bar = fastLP_forward(ub, slopes, intersections, u_init, u_end);

%% Perform backward phase
for i=N-1 : 1
    %load slopes & intersections of linear inequalities, remove -1s
    F_ki = slopes{2}(i, :);   F_ki=F_ki(F_ki>=0);
    F_bi = intersections{2}(i, :);    F_bi=F_bi(F_bi>=0);
    r_i = length(F_ki);
    %backward restriction
    B_i = zeros(1, r_i);
    for j=1:r_i
        B_i(j) = f_eval(u_bar(i+1), F_ki(j), F_bi(j));
    end
    u_bar(i) = min(u_bar(i), min(B_i), ub(i));
end