function adaptive_exper

    Fehlberg = struct();
    Fehlberg.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
    Fehlberg.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
    25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    Fehlberg.A = [0,0,0,0,0,0;...
    1/4, 0,0,0,0,0;...
    3/32, 9/32, 0,0,0,0;...
    1932/2197, -7200/2197, 7296/2197, 0,0,0;...
    439/216, -8, 3680/513, -845/4104, 0,0;...
    -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];

    tspan = [0 30];
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1;
    V0 = [x0;y0;dxdt0;dydt0];
    V_true = compute_planetary_motion(tspan(end), V0, orbit_params);
    rate_func = @(t,V) gravity_rate_func(t,V,orbit_params);

    p = 5;
    h_ref = 0.05;
    h_list = logspace(-3, -1, 15);    
    error_list = logspace(-12, -3, 15); 

    nH = numel(h_list);
    nE = numel(error_list);

    global_err_fixed = zeros(1, nH);
    h_avg_list_fixed = zeros(1, nH);
    num_evals_list_fixed = zeros(1, nH);
    Fehlberg_fixed = Fehlberg;
    Fehlberg_fixed.B = Fehlberg.B(1,:)';

    for i = 1:nH
        [~, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration( ...
            rate_func, tspan, V0, h_list(i), Fehlberg_fixed);
        h_avg_list_fixed(i) = h_avg;
        num_evals_list_fixed(i) = num_evals;
        global_err_fixed(i) = norm(X_list(end,:)' - V_true);
    end

    global_err_adapt = zeros(1, nE);
    h_avg_list_adapt = zeros(1, nE);
    num_evals_list_adapt = zeros(1, nE);

    fail_rate_list = zeros(1, nE);

    for i = 1:nE
        [~, X_list, h_avg, num_evals, fail_rate] = explicit_RK_variable_step_integration(rate_func, tspan, V0, h_ref, Fehlberg, p, error_list(i));
        h_avg_list_adapt(i) = h_avg;
        num_evals_list_adapt(i) = num_evals;
        fail_rate_list(i) = fail_rate;
        global_err_adapt(i) = norm(X_list(end,:)' - V_true);
    end 

    idx_fixed = global_err_fixed > 1e-13;
    idx_adapt = global_err_adapt > 1e-13;
    mdl_fixed = polyfit(log(h_avg_list_fixed(idx_fixed)), log(global_err_fixed(idx_fixed)), 1);
    mdl_adapt = polyfit(log(h_avg_list_adapt(idx_adapt)), log(global_err_adapt(idx_adapt)), 1);
    p_fixed = mdl_fixed(1);
    p_adapt = mdl_adapt(1);

    % Plot global vs avg step size

    h_fit = logspace(log10(min([h_avg_list_fixed, h_avg_list_adapt])), ...
                     log10(max([h_avg_list_fixed, h_avg_list_adapt])), 50);
    figure();
    loglog(h_avg_list_fixed, global_err_fixed, 'b.', 'MarkerSize', 20);
    hold on;
    loglog(h_avg_list_adapt, global_err_adapt, 'm.', 'MarkerSize', 20);

    title('Global Truncation Error vs Step Size', 'Interpreter', 'Latex', 'FontSize', 22);
    xlabel('Average Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 18);
    ylabel('Global Truncation Error (-)', 'Interpreter', 'Latex', 'FontSize', 18);
    legend('Fixed step', 'Adaptive step', 'Interpreter', 'Latex', 'Location', 'southeast');
    set(gca, 'FontSize', 14)
    hold off;

    % Plot global vs num function

    figure();
    loglog(num_evals_list_fixed, global_err_fixed, 'b.', 'MarkerSize', 20);
    hold on;
    loglog(num_evals_list_adapt, global_err_adapt, 'm.', 'MarkerSize', 20);

    title('Global Truncation Error vs Function Evaluations', 'Interpreter', 'Latex', 'FontSize', 22);
    xlabel('Number of Function Evaluations (-)', 'Interpreter', 'Latex', 'FontSize', 18);
    ylabel('Global Truncation Error (-)', 'Interpreter', 'Latex', 'FontSize', 18);
    legend('Fixed step', 'Adaptive step', 'Interpreter', 'Latex', 'Location', 'northeast');
    set(gca, 'FontSize', 14)
    hold off

    % Plot fail rate vs avg step size 
    figure();
    semilogx(h_avg_list_adapt, fail_rate_list, 'c.', 'MarkerSize', 20)

    title('Step Failure Rate vs Average Step Size', 'Interpreter', 'Latex', 'FontSize', 22)
    xlabel('Average Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Failure Rate (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    set(gca, 'FontSize', 14)

    % Position, Velocity
    error = 1e-4;
    [t_list, X_list, h_avg, num_evals, fail_rate] = explicit_RK_variable_step_integration( ...
        rate_func, tspan, V0, h_ref, Fehlberg, p, error);
    t_range = linspace(tspan(1), tspan(2), 1000);
    V = compute_planetary_motion(t_range, V0, orbit_params);

    % Plot Position
    % x
    figure();
    subplot(2,1,1); 
    plot(t_range, V(:,1), '--', 'LineWidth', 2);
    hold on
    plot(t_list, X_list(:,1), 'ro-', 'markerfacecolor','k','markeredgecolor','k','markersize', 2);
    
    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$x(t)$', 'Interpreter', 'Latex'); 
    title('Position vs Time', 'Interpreter', 'Latex')
    legend('Exact Solution', 'Fehlburg Adaptive')
    hold off;

    % y
    subplot(2,1,2); 
    plot(t_range, V(:,2), '--')
    hold on;
    plot(t_list, X_list(:,2), 'ro-', 'markerfacecolor','k','markeredgecolor','k','markersize',3);
    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$y(t)$', 'Interpreter', 'Latex')
    legend('Exact Solution', 'Fehlburg Adaptive')
    hold off;
    
    % Plot Velocity
    figure();   
    subplot(2,1,1); hold on
    plot(t_range, V(:,3), '--')
    plot(t_list, X_list(:,3), 'ro-', 'markerfacecolor','k','markeredgecolor','k','markersize',2);
    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$v_x(t)$', 'Interpreter', 'Latex'); 
    title('Velocity vs Time', 'Interpreter', 'Latex');
    legend('Exact Solution', 'Fehlburg Adaptive')
    hold off;

    subplot(2,1,2); 
    plot(t_range, V(:,4), '--');
    hold on;
    plot(t_list, X_list(:,4), 'ro-', 'markerfacecolor','k','markeredgecolor','k','markersize',2);
    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$v_y(t)$', 'Interpreter', 'Latex')
    legend('Exact Solution', 'Fehlburg Adaptive')
    hold off;
    
    % step size vs dist
    h_list = diff(t_list);
    r = sqrt(X_list(1:end-1,1).^2 + X_list(1:end-1,2).^2); 
    
    figure();
    loglog(r, h_list, 'k.', 'MarkerSize', 12)

    xlabel('Distance from Sun: r (-)', 'Interpreter', 'Latex', 'FontSize', 18);
    ylabel('Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 18);
    title('Step Size vs Distance from Sun', 'Interpreter', 'Latex', 'FontSize', 22);

end