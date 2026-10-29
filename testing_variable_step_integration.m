function testing_variable_step_integration

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
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    
    t_range = linspace(tspan(1),tspan(2),500);
    V = compute_planetary_motion(t_range,V0,orbit_params);
    rate_func = @(t,V) gravity_rate_func(t,V,orbit_params);

    p = 5;                 
    error_desired = 1e-5;
    h_ref = 0.05;

    [t_Fehl, X_Fehl, ~, ~, ~] = explicit_RK_variable_step_integration(rate_func, tspan, V0, h_ref, Fehlberg, p, error_desired);  

    % Plot - both x and y 
    figure(); 
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on

    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    hold on;
    plot(V(:,1),V(:,2),'k','LineWidth',2);

    plot(X_Fehl(:,1),X_Fehl(:,2),'--','LineWidth', 1.5);

    xlabel('$x$ position', 'Interpreter', 'Latex');
    ylabel('$y$ position', 'Interpreter', 'Latex');
    title('Planetary Motion Plot', 'Interpreter', 'Latex', 'FontSize', 22)

    legend('Sun', 'Exact solution', 'Fehlburg', 'Location', 'best', 'Interpreter', 'Latex');
    hold off;

    % Plot subfigures 
    figure();

    subplot(2,1,1); 
    plot(t_range, V(:,1), 'k-', 'LineWidth', 2);
    hold on;
    plot(t_Fehl, X_Fehl(:,1),'--', 'LineWidth', 1.5);

    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$x(t)$', 'Interpreter', 'Latex');
    title('Planetary Motion x-component');
    legend('Exact solution', 'Fehlburg', 'Location', 'best', 'Interpreter', 'Latex');

    subplot(2,1,2); 
    plot(t_range,V(:,2),'k','LineWidth',2);
    hold on; 
    plot(t_Fehl,X_Fehl(:,2),'--','markerfacecolor','k','markeredgecolor','k','markersize',3);
    xlabel('Time', 'Interpreter', 'Latex');
    ylabel('$y(t)$', 'Interpreter', 'Latex');
    title('Planetary Motion y-component');
    legend('Exact solution', 'Fehlburg', 'Location', 'best', 'Interpreter', 'Latex');
end