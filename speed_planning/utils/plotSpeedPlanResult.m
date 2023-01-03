function plotSpeedPlanResult(dq_t, ddq_t, Tau, t_sample, constraints)
% PLOTSPEEDPLANRESULT plot speed planning results

dim = size(dq_t, 2);
cmap = cool(dim); %define color map

figure('Position',[500,50, 1200,900]);

subplot(3,1,1)
hold on;
for i=1:dim
    p(i) = plot(t_sample, dq_t(:, i), '-', 'LineWidth', 1.5, 'Color', cmap(i, :)); 
    tag{i} = sprintf('$\\dot{q}_%d$', i);
    plot(t_sample, constraints.Phi(i) * ones(size(t_sample)), '--', 'LineWidth', 0.8, 'Color', cmap(i, :));
    plot(t_sample, - constraints.Phi(i) * ones(size(t_sample)), '--', 'LineWidth', 0.8, 'Color', cmap(i, :));
end
title('Velocity');
legend(p, tag, "Interpreter","latex","FontSize", 12);

subplot(3,1,2)
hold on;
for i=1:dim
    p(i) = plot(t_sample(1:end-1), ddq_t(:, i), '-', 'LineWidth', 1.5, 'Color', cmap(i, :)); 
    tag{i} = sprintf('$\\ddot{q}_%d$', i);
    plot(t_sample(1:end-1), constraints.Alpha(i) * ones(size(t_sample(1:end-1))), '--', 'LineWidth', 0.8, 'Color', cmap(i, :));
    plot(t_sample(1:end-1), - constraints.Alpha(i) * ones(size(t_sample(1:end-1))), '--', 'LineWidth', 0.8, 'Color', cmap(i, :));
end
title('Acceleration');
legend(p, tag, "Interpreter","latex","FontSize", 12);

subplot(3,1,3)
 hold on;
for i=1:dim
    p(i) = plot(t_sample(1:end-1), Tau(:, i), '-', 'LineWidth', 1.5, 'Color', cmap(i, :)); 
    tag{i} = sprintf('$\\tau_%d$', i);
    plot(t_sample(1:end-1), constraints.Mu(i) * ones(size(t_sample(1:end-1))), '--', 'LineWidth', 0.8, 'Color', cmap(i, :));
    plot(t_sample(1:end-1), - constraints.Mu(i) * ones(size(t_sample(1:end-1))), '--', 'LineWidth', 0.8, 'Color', cmap(i, :));
end
title('Force');
legend(p, tag, "Interpreter","latex","FontSize", 12);
