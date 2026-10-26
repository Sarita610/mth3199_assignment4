function plotting

    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1;
    V0 = [x0;y0;dxdt0;dydt0];
    
    tspan = [0 30];
    h = 0.25;

    t_range = linspace(tspan(1), tspan(2), 100);
    V = compute_planetary_motion(t_range, V0, orbit_params);
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

    [t_eul, X_eul] = explicit_RK_fixed_step_integration(rate_func, tspan, V0, h, BT_eul);
    [t_midpt, X_midpt] = explicit_RK_fixed_step_integration(rate_func, tspan, V0, h, BT_midpt);
    [t_third, X_third] = explicit_RK_fixed_step_integration(rate_func, tspan, V0, h, BT_third);
    

    %Plot - both x and y 
    figure();
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on

    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V(:,1), V(:,2), 'k', 'LineWidth', 2);

    plot(X_eul(:,1), X_eul(:,2), '--', 'LineWidth', 1.5);
    plot(X_midpt(:,1), X_midpt(:,2), ':', 'LineWidth', 1.5);
    plot(X_third(:,1), X_third(:,2), '-.', 'LineWidth', 1.5);

    xlabel('$x$ position', 'Interpreter', 'Latex');
    ylabel('$y$ position', 'Interpreter', 'Latex');
    title('Planetary Motion Plot', 'Interpreter', 'Latex', 'FontSize', 22)

    legend('Sun', 'Exact solution', 'Forward Euler', 'Explicit midpoint', "Kutta's third-order", 'Location', 'best', 'Interpreter', 'Latex');
    hold off;

    % Plot subfigures
    figure();
    
    subplot(2,1,1);
    plot(t_range, V(:,1), 'k-', 'LineWidth', 2);
    hold on;
    plot(t_eul, X_eul(:,1), '--', 'LineWidth', 1.5);
    plot(t_midpt, X_midpt(:,1), ':', 'LineWidth', 1.5);
    plot(t_third, X_third(:,1), '-.', 'LineWidth', 1.5);

    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$x(t)$', 'Interpreter', 'Latex');
    title('Planetary Motion x-component');
    legend('Exact solution', 'Forward Euler', 'Explicit midpoint', "Kutta's third-order", 'Location', 'best', 'Interpreter', 'Latex');

    subplot(2,1,2);
    plot(t_range, V(:,2), 'k-', 'LineWidth', 2);
    hold on;
    plot(t_eul, X_eul(:,2), '--', 'LineWidth', 1.5);
    plot(t_midpt, X_midpt(:,2), ':', 'LineWidth', 1.5);
    plot(t_third, X_third(:,2), '-.', 'LineWidth', 1.5);

    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$y(t)$', 'Interpreter', 'Latex');
    title('Planetary Motion y-component');
    legend('Exact solution', 'Forward Euler', 'Explicit midpoint', "Kutta's third-order", 'Location', 'best', 'Interpreter', 'Latex');



end