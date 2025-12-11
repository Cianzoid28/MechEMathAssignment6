function nat_freq_string_simulation(string_params)
    close all;
    % num_masses = 4;
    % total_mass = 1;
    % tension_force = 15;   % N
    % string_length = 10;  % m
    % damping_coeff = string_length / tension_force;
    % 
    % dx = string_length / (num_masses + 1);
    % 
    % amplitude_Uf = 0.01;
    % omega_Uf = nat_freq;
    % 
    % % boundary forcing at right end
    % Uf_func    = @(t) amplitude_Uf * cos(omega_Uf * t);
    % dUfdt_func = @(t) -omega_Uf * amplitude_Uf * sin(omega_Uf * t);
    % 
    % % pack params
    % string_params = struct();
    % string_params.n = num_masses;
    % string_params.M = total_mass;
    % string_params.Uf_func = Uf_func;
    % string_params.dUfdt_func = dUfdt_func;
    % string_params.Tf = tension_force;
    % string_params.L = string_length;
    % string_params.c = damping_coeff;
    % string_params.dx = dx;

    my_rate_func = @(t, V) string_rate_func01(t, V, string_params);

    % initial conditions
    %U0     = [-0.01; 0.01; -0.01; 0.01];
    %dUdt0  = [-0.01; 0.01; -0.01; 0.01];
    U0     = zeros([string_params.n, 1]);
    dUdt0     = zeros([string_params.n, 1]);
    V0 = [U0; dUdt0];  % 8×1 vector

    tspan = [0 128];

    % run simulation
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0);

    % extract displacements
    Ulist = Vlist(:, 1:string_params.n);

    % plot
    figure;
    plot(tlist, Ulist, 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Vertical Displacement U_i');
    title('Masses on String');
    legend('U_1','U_2','U_3','U_4');


    %anime
    x = linspace(0, string_params.L, string_params.n+2);
    
    % initial Y-values
    U_full = [0, Ulist(1,:), string_params.Uf_func(0)];
    
    % figure
    % hold on;
    % 
    % hLine = plot(x, U_full, 'k:', 'LineWidth', 1.5);
    % 
    % hPtsRed = plot(x(2:end-1), Ulist(1,:), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    % 
    % hPtsBlue = plot([x(1), x(end)], [U_full(1), U_full(end)], 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    % 
    % xlim([-0.25 string_params.L+0.25])
    % ylim([-0.15 0.15])
    % 
    % for k = 1:length(tlist)
    %     % update full string with boundary points
    %     U_full = [0, Ulist(k,:), string_params.Uf_func(tlist(k))];
    % 
    %     hLine.YData = U_full;
    % 
    %     % update interior points
    %     hPtsRed.YData  = Ulist(k,:);
    % 
    %     % update boundary points
    %     hPtsBlue.YData = [U_full(1), U_full(end)];
    % 
    %     title(sprintf('t = %.3f s', tlist(k)));
    %     pause(0.01);
    %     drawnow;
    % end

    figure
    x
    string_params.Ur(:, 1)
    plot(x, string_params.Ur(:, 1))

end
