function eval_embedded_RK

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

    t_ex = 0.492;
    X_ex = solution02(t_ex);
    %tspan = [0 30];

    h_list = logspace(-5, 2, 100);

    error_list_fehlberg1 = zeros(size(h_list));
    error_list_fehlberg2 = zeros(size(h_list));
    difference_list_fehlberg = zeros(size(h_list));
    exact_difference_list_fehlberg = zeros(size(h_list));

for i = 1:length(h_list)

        [XB1, XB2, num_evals] = RK_step_embedded(@rate_func02, t_ex, X_ex, h_list(i), Fehlberg);
        X_fehlberg = solution02(t_ex+h_list(i));
        X_fehlberg = X_fehlberg(:);
        local_error_1 = norm(XB1 - X_fehlberg);
        error_list_fehlberg1(i) = local_error_1;
        local_error_2 = norm(XB2 - X_fehlberg);
        error_list_fehlberg2(i) = local_error_2;
        difference = norm(XB2 - XB1);
        difference_list_fehlberg(i) = difference;
        exact_difference_list_fehlberg(i) = norm(X_fehlberg - X_ex);
end
    idx = h_list>= 1e-5 & h_list<=3e-1;
    floor_val = 1e-13;
    hwin = h_list >= 2e-2 & h_list <= 3e-1;

    idx1 = hwin & error_list_fehlberg1 > floor_val;
    idx2 = hwin & error_list_fehlberg2 > floor_val;
    idxd = hwin & difference_list_fehlberg > floor_val;

    mdl_local_fehlberg1 = polyfit(log(h_list(idx1)), log(error_list_fehlberg1(idx1)), 1);
    mdl_local_fehlberg2 = polyfit(log(h_list(idx2)), log(error_list_fehlberg2(idx2)), 1);
    mdl_difference      = polyfit(log(h_list(idxd)), log(difference_list_fehlberg(idxd)), 1);
    mdl_true_difference = polyfit(log(h_list(idx)), log(exact_difference_list_fehlberg(idx)), 1);
    p_local_fehlberg1 = mdl_local_fehlberg1(1)
    p_local_fehlberg2 = mdl_local_fehlberg2(1)
    p_difference      = mdl_difference(1)
    p_true_difference = mdl_true_difference(1);

% Fehlberg local truncation plots
    figure();
    loglog(h_list, error_list_fehlberg1, h_list, error_list_fehlberg2, ...
           h_list, difference_list_fehlberg, h_list, exact_difference_list_fehlberg, 'LineWidth', 1.5);
    hold on
    loglog(h_list(idx), exp(polyval(mdl_local_fehlberg1, log(h_list(idx)))), 'k--');
    loglog(h_list(idx), exp(polyval(mdl_local_fehlberg2, log(h_list(idx)))), 'k--');
    xlabel('step size $h$','Interpreter','Latex');
    ylabel('error','Interpreter','Latex');
    title('Fehlberg local truncation error');

% Fehlberg local truncation --
    figure();
    loglog(h_list, error_list_fehlberg1, '.', 'MarkerSize', 16);
    hold on
    loglog(h_list, error_list_fehlberg2, '.', 'MarkerSize', 14);
    loglog(h_list, difference_list_fehlberg, '.', 'MarkerSize', 10);
    loglog(h_list, exact_difference_list_fehlberg, '.', 'MarkerSize', 10);
    % loglog(h_list(idx), exp(polyval(mdl_local_fehlberg1, log(h_list(idx)))), 'k--');
    % loglog(h_list(idx), exp(polyval(mdl_local_fehlberg2, log(h_list(idx)))), 'k--');
    xlabel('Step Size $h$ (-)','Interpreter','Latex', 'FontSize',16);
    ylabel('Local Truncation Error (-)','Interpreter','Latex', 'FontSize',16);
    title('Fehlberg Local Truncation Error', 'FontSize', 20);
    legend('Fehlberg: 5th-Order Error', 'Fehlberg: 4th-Order Error', '$X_{\mathrm{diff}}$ : $|XB1-XB2|$', '$|f(t_{\mathrm{ref}}+h)-f(t_{\mathrm{ref}})|$', 'Interpreter', 'Latex', 'Location', 'northwest', 'FontSize', 12 )

end 