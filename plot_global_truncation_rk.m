function plot_global_truncation_rk
    % Plotting parameters
    tspan = [0 10];
    h_list = logspace(-4, -1, 20);

 
    %system params
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];

    V_true = compute_planetary_motion(tspan(end), V0, orbit_params);
    rate_func = @(t,V) gravity_rate_func(t,V,orbit_params);

    % Forward Euler Method
    BT_eul.A = 0;
    BT_eul.B = 1;
    BT_eul.C = 0;


    % Explicit midpoint method
    BT_midpt.A = [0, 0; 1/2, 0];
    BT_midpt.B = [0; 1];
    BT_midpt.C = [0; 1/2];

    % Kutta's third-order method
    BT_third.A = [0, 0, 0; 1/2, 0, 0; -1, 2, 0];
    BT_third.B = [1/6; 2/3; 1/6];
    BT_third.C = [0; 1/2; 1];

    BT_list = {BT_eul, BT_midpt, BT_third};
    p_expected = [1; 2; 3];
    names = {'Forward Euler'; 'Explicit Midpoint'; "Kutta's Third-Order"};

    %initilaize variables
    nM = numel(BT_list);
    nH = numel(h_list);
    global_err = zeros(nM, nH);
    h_avg_list = zeros(nM, nH);
    num_evals_list = zeros(nM, nH);
    p_global = zeros(nM, 1);
    k_global = zeros(nM, 1);

    for m = 1:nM
        for i = 1:nH
            [~, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration( ...
                rate_func, tspan, V0, h_list(i), BT_list{m});
            h_avg_list(m,i) = h_avg;
            num_evals_list(m,i) = num_evals;
            global_err(m,i) = norm(X_list(end,:)' - V_true);
        end
    end
    
    fit_range = [1e-3, 1e-1;    % Forward Euler
                 1e-3, 1e-1;    % Explicit midpoint
                 1e-3, 1e-1];   % Kutta's third-order
    for m = 1:nM
        idx = h_avg_list(m,:) >= fit_range(m,1) & h_avg_list(m,:) <= fit_range(m,2) ...
              & global_err(m,:) > 1e-13;
        mdl = polyfit(log(h_avg_list(m,idx)), log(global_err(m,idx)), 1);
        p_global(m) = mdl(1);
        k_global(m) = exp(mdl(2));
    end

    % Table of p values
    p_table = table(names, p_global, p_expected, ...
        'VariableNames', {'Method', 'p_measured', 'p_expected'})
    
       % Global error vs step size, with fit lines
    h_fit = logspace(log10(min(h_avg_list(:))), log10(max(h_avg_list(:))), 50);
 
    figure();
    loglog(h_avg_list(1,:), global_err(1,:), 'b.', 'MarkerSize', 20)
    hold on;
    loglog(h_avg_list(2,:), global_err(2,:), 'c.', 'MarkerSize', 20)
    loglog(h_avg_list(3,:), global_err(3,:), 'm.', 'MarkerSize', 20)
 
    % Fit lines (k*h^p), one per method
    loglog(h_fit, k_global(1)*h_fit.^p_global(1), 'b-', 'LineWidth', 1.5, 'HandleVisibility', 'off')
    loglog(h_fit, k_global(2)*h_fit.^p_global(2), 'c-', 'LineWidth', 1.5, 'HandleVisibility', 'off')
    loglog(h_fit, k_global(3)*h_fit.^p_global(3), 'm-', 'LineWidth', 1.5, 'HandleVisibility', 'off')
 
    title('Global Truncation Error vs Step Size', 'Interpreter', 'Latex', 'FontSize', 22)
    xlabel('Average Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Global Truncation Error (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    legend('Forward Euler', 'Explicit Midpoint', "Kutta's Third-Order", ...
        'Interpreter', 'Latex', 'Location', 'southeast')
    set(gca, 'FontSize', 14)
    hold off
    
    
   