function string_simulation(string_params)
    my_rate_func = @(t, V) string_rate_func01(t, V, string_params);

    U0 = zeros([string_params.n, 1]);
    dUdt0 = zeros([string_params.n, 1]);
    V0 = [U0; dUdt0];  % 8×1 vector
    
    n = ceil(string_params.omega_Uf);
    while mod(n, 4) ~= 1
        n = n + 1;  % Increment n to find the next valid value
    end
    n*pi / (2*string_params.omega_Uf)
    tspan = [0, n*pi / (2*string_params.omega_Uf)];

    % run simulation
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0);

    % extract displacements
    Ulist = Vlist(:, 1:string_params.n);

    % plot
    % figure();
    % plot(tlist, Ulist, 'LineWidth', 1.5);
    % xlabel('Time (s)');
    % ylabel('Vertical Displacement U_i');
    % title('Masses on String');
    % legend('U_1','U_2','U_3','U_4');


    %anime
    x = linspace(0, string_params.L, string_params.n+2);
    
    % initial Y-values
    U_full = [0, Ulist(1,:), string_params.Uf_func(0)];
    
    figure
    hold on;

    modeShape = plot(x, [0; string_params.U_r; 0], 'b:', 'LineWidth', 1.5);
    
    hLine = plot(x, U_full, 'k:', 'LineWidth', 1.5);
    
    hPtsRed = plot(x(2:end-1), Ulist(1,:), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    
    hPtsBlue = plot([x(1), x(end)], [U_full(1), U_full(end)], 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    
    xlim([-0.25 string_params.L+0.25])
    ylim([-0.3 0.3])
    %ylim([-max((abs(string_params.U_r))), max(abs(string_params.U_r))])
    
    %for k = 1:length(tlist)
    %for k = length(tlist)-1:length(tlist)
    k = 1;
    while k <= length(tlist)
        if abs(Ulist(k,end)) >= abs(string_params.U_r(end)) && sign(Ulist(k,end)) == sign(string_params.U_r(end))
            break
        end
        modeShape.YData = [0; string_params.U_r; 0];

        % update full string with boundary points
        U_full = [0, Ulist(k,:), string_params.Uf_func(tlist(k))];
        
        hLine.YData = U_full;
        
        % update interior points
        hPtsRed.YData  = Ulist(k,:);
        
        % update boundary points
        hPtsBlue.YData = [U_full(1), U_full(end)];
        title(sprintf('Frequency: %.3f', string_params.omega_Uf))
        subtitle(sprintf('t = %.3f s', tlist(k)));
        %pause(0.005);
        drawnow;
        k = k + 1;  % Increment k to move to the next time step
    end

end
