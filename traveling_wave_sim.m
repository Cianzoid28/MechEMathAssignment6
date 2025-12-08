function traveling_wave_sim();
    
    %Parameters
    num_masses = 300;
    string_length = 12; %M
    total_mass = 1; %kg
    tension_force = 100; %N
    damping_coeff = 0;
    dx = string_length/(num_masses+1);

    %Pulse Parameters
    pulse_width = 0.05; %duration of pulse
    pulse_height = 0.01; %peak displacement

    %Pulse type
    use_triangle = true; %Change to false for B-Spline

    rho = total_mass/string_length; %linear density kg/m
    c = sqrt(tension_force / rho); %wave speed m/s
    time_one = string_length / c; %time it takes for one wave
    T_final = 4 * time_one; %4 round trips

    %Time step
    dt = 0.2 * dx/c;
    tspan = 0:dt:T_final;

    if use_triangle == true
        Uf_func = @(t_in) triangle_pulse(t_in,pulse_width,pulse_height);
        dUfdt_func = @(t_in) triangle_pulse_derivative(t_in,pulse_width,pulse_height);
        pulse_name = 'Triangle';
    else
        Uf_func = @(t_in) b_spline_pulse(t_in, pulse_width, pulse_height);
        dUfdt_func = @(t_in) b_spline_pulse_derivative(t_in,pulse_width,pulse_height);
        pulse_name = 'B-spline';
    end

    string_params = struct();
    string_params.n            = num_masses;
    string_params.M            = total_mass;
    string_params.Tf           = tension_force;
    string_params.L            = string_length;
    string_params.c            = damping_coeff;
    string_params.dx           = dx;
    string_params.Uf_func      = Uf_func;
    string_params.dUfdt_func   = dUfdt_func;
 


