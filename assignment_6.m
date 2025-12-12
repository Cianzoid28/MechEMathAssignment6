function assignment_6()
    close all;
    %experiment_1();
    %experiment_2();
    %experiment_3();
    experiment_4();
end

function experiment_1()

    n=5;
    mode_index = 4;

    string_params.n = n; %number of masses
    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = .0001; %damping coefficient
    string_params.dx = string_params.L/(string_params.n+1); %horizontal spacing between masses

    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
    %Use MATLAB to solve the generalized eigenvalue problem
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);

    mode_shape_LA = Ur_mat(:,mode_index);

    omega_n = sqrt(-lambda_mat(mode_index,mode_index));

    omega = omega_n; %rads/sec
    A = 3;

    string_params.Uf_func = @(t_in) A*sin(omega*t_in); %function describing motion of end point
    string_params.dUfdt_func = @(t_in) -omega*A*cos(omega*t_in); %time derivative of Uf
    string_params.Uf_amplitude = A;

    U0 = zeros(string_params.n,1);
    dUdt0 = zeros(string_params.n,1);

    V0 = [U0;dUdt0];

    tlist = linspace(0,20*(2*pi)/omega, 10000+1);

    my_rate_func = @(t_in,V_in) string_rate_func02(t_in,V_in,string_params);
    
    [tlist, Vresult] = ode45(my_rate_func,tlist,V0);

    animate_string_with_mode_shape(tlist,Vresult,string_params,mode_shape_LA);

end

function experiment_2()

    n=200;

    string_params.n = n; %number of masses
    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = .0001; %damping coefficient
    string_params.dx = string_params.L/(n+1); %horizontal spacing between masses


    string_params.Uf_func = @(t_in) 0; %function describing motion of end point
    string_params.dUfdt_func = @(t_in) 0; %time derivative of Uf
    string_params.Uf_amplitude = 0;

    rho = string_params.M/string_params.L;

    c = sqrt(string_params.Tf/rho);

    x_list = linspace(0,string_params.L,n+2)';
    x_list = x_list(2:end-1);

    w_pulse = string_params.L/5;
    h_pulse = 7;

    x_shift = 0;
    U0 = b_spline_pulse(x_list-x_shift,w_pulse,h_pulse);
    dUdt0 = -c*b_spline_pulse_derivative(x_list-x_shift,w_pulse,h_pulse);

    V0 = [U0;dUdt0];

    tlist = linspace(0,6*string_params.L/c,10000+1);

    my_rate_func = @(t_in,V_in) string_rate_func02(t_in,V_in,string_params);
    
    [tlist, Vresult] = ode45(my_rate_func,tlist,V0);

    animate_string_travelling_wave(tlist,Vresult,string_params,x_shift+w_pulse/2, c);

end

function experiment_3()
    

      n=200;
    mode_index = 10;

    string_params.n = n; %number of masses
    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = .0001; %damping coefficient
    string_params.dx = string_params.L/(string_params.n+1); %horizontal spacing between masses

    rho = string_params.M/string_params.L;
    c = sqrt(string_params.Tf/rho);
   
    omega_n_spatial = pi*mode_index/string_params.L;
    omega_n = c*omega_n_spatial;

    x_list = linspace(0,string_params.L,n+2)';
    x_list = x_list(2:end-1);

    mode_shape_WE = sin(omega_n_spatial*x_list);

    omega = omega_n; %rads/sec
    A = 3;

    string_params.Uf_func = @(t_in) A*sin(omega*t_in); %function describing motion of end point
    string_params.dUfdt_func = @(t_in) -omega*A*cos(omega*t_in); %time derivative of Uf
    string_params.Uf_amplitude = A;

    U0 = zeros(string_params.n,1);
    dUdt0 = zeros(string_params.n,1);

    V0 = [U0;dUdt0];

    tlist = linspace(0,20*(2*pi)/omega, 5000+1);

    my_rate_func = @(t_in,V_in) string_rate_func02(t_in,V_in,string_params);
    
    [tlist, Vresult] = ode45(my_rate_func,tlist,V0);

    animate_string_with_mode_shape(tlist,Vresult,string_params,mode_shape_WE);

end

function experiment_4()
    
    mode_index = 3;

    string_params.M = 10; %total mass attached to the string
    string_params.Tf = 2; %tension in string
    string_params.L = 7; %length of string
    string_params.c = .0001; %damping coefficient

    rho = string_params.M/string_params.L;
    c = sqrt(string_params.Tf/rho);

    omega_n_spatial = pi*mode_index/string_params.L;
    omega_n = c*omega_n_spatial;

    x_list_continuous = linspace(0,string_params.L,1000);
    mode_shape_WE = sin(omega_n_spatial*x_list_continuous);


    n_list = [5,10,30,50];

    for k = 1:length(n_list)
        n = n_list(k);

        string_params.n = n; %number of masses
        string_params.dx = string_params.L/(string_params.n+1); %horizontal spacing between masses
    
        [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
        
        %Use MATLAB to solve the generalized eigenvalue problem
        [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
    
        mode_shape_LA = [0;Ur_mat(:,n+1-mode_index);0];
        mode_shape_LA = mode_shape_LA/max(abs(mode_shape_LA));
    
        omega_n = sqrt(-lambda_mat(n+1-mode_index,n+1-mode_index));
    
        x_list = linspace(0,string_params.L,n+2)';
    
        subplot(length(n_list),1,k);
        hold on
        plot(x_list_continuous, mode_shape_WE,'b');
        plot(x_list,mode_shape_LA,'o-', 'Color','k', 'markerfacecolor','r','markersize',4);
    end

end

%b-spline pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = b_spline_pulse(t,w,h)
    t = 4*t/w;
    b3 = (0<=t).*(t<1).*(t.^3)/4;
    t = t-1;
    b2 = (0<=t).*(t<1).*(-3*t.^3+3*t.^2+3*t+1)/4;
    t = t-1;
    b1 = (0<=t).*(t<1).*(3*t.^3-6*t.^2+4)/4;
    t = t-1;
    b0 = (0<=t).*(t<1).*(-t.^3+3*t.^2-3*t+1)/4;
    res = h*(b0+b1+b2+b3);
end
%b-spline pulse function (derivative)
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: derivative of pulse evaluated at t
function res = b_spline_pulse_derivative(t,w,h)
    t = 4*t/w;
    b3 = (0<=t).*(t<1).*(3*t.^2)/4;
    t = t-1;
    b2 = (0<=t).*(t<1).*(-9*t.^2+6*t+3)/4;
    t = t-1;
    b1 = (0<=t).*(t<1).*(9*t.^2-12*t)/4;
    t = t-1;
    b0 = (0<=t).*(t<1).*(-3*t.^2+6*t-3)/4;
    res = (4*h/w)*(b0+b1+b2+b3);
    end

%triangle pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=w)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = triangle_pulse(t,w,h)
t = t*(2/w);
res = 1-min(1*abs(t-1),1);
res = h*res;
end

%triangle pulse function (derivative)
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=w)
%h: height of pulse
%OUTPUTS:
%res: derivative of pulse evaluated at t
function res = triangle_pulse_derivative(t,w,h)
    t = t*(2/w);
    res = -sign(t-1).*(abs(t-1)<1);
    res = (2*h/w)*res;
end


function animate_string_travelling_wave(tlist,Vresult,string_params,x_centroid0, c)
    n = string_params.n;
    L = string_params.L;

    x_list = linspace(0,L,n+2);

    maxU = max(max(Vresult(:,1:n)));
    minU = min(max(Vresult(:,1:n)));

    h = max([abs(maxU),abs(minU), string_params.Uf_amplitude]);
    axis([0,L,-1.1*h,1.1*h]);
    hold on

    string_plot = plot(0,0,'o-', 'Color','k', 'markerfacecolor','r','markersize',4);

    centroid_plot = plot(0,0,'k',LineStyle=':');

    xlabel('x');
    ylabel('U(t,x)');
    title('Plot of vibrating string');

    for k = 1:length(tlist)

        t = tlist(k);
        U = Vresult(k,1:n);
        Uf = string_params.Uf_func(t);

        U_padded = [0,U,Uf];

        %short blurb showing how to find x-coord of tracking line
        %x = x-coord of tracking line, t = current time
        %c = wave speed, w = pulse width (in time), L = string length
        x_centroid = c*t+x_centroid0;
        x_centroid = mod(x_centroid,2*L);
        if x_centroid > L
            x_centroid = 2*L - x_centroid;
        end

        set(centroid_plot,'xdata',[x_centroid,x_centroid],'ydata',[-1.1*h,1.1*h])
        

        set(string_plot,'xdata', x_list, 'ydata', U_padded);
        drawnow;
    end



end

function animate_string_with_mode_shape(tlist,Vresult,string_params,mode_shape)
    n = string_params.n;
    L = string_params.L;

    x_list = linspace(0,L,n+2);

    maxU = max(max(Vresult(:,1:n)));
    minU = min(max(Vresult(:,1:n)));

    h = max([abs(maxU),abs(minU), string_params.Uf_amplitude]);
    axis([0,L,-1.1*h,1.1*h]);
    hold on

    string_plot = plot(0,0,'o-', 'Color','k','linewidth', 2, 'markerfacecolor','r','markersize',4);

    mode_shape_plot = plot(0,0,'o-', 'Color','b','linewidth', 2, 'markerfacecolor','k','markersize',4);

    mode_shape_padded = [0;mode_shape;0];

    scale_factor = max(abs(maxU),abs(minU))/max(abs(mode_shape_padded));
    mode_shape_padded = mode_shape_padded*scale_factor;
    set(mode_shape_plot,'xdata',x_list,'ydata',mode_shape_padded);

    xlabel('x');
    ylabel('U(t,x)');
    title('Plot of vibrating string');

    for k = 1:length(tlist)
        t = tlist(k);
        U = Vresult(k,1:n);
        Uf = string_params.Uf_func(t);

        U_padded = [0,U,Uf];
        

        set(string_plot,'xdata', x_list, 'ydata', U_padded);
        drawnow;
    end



end

%build the mass and stiffness matrices that describe the 2nd order system.
%INPUTS
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
%OUTPUTS
%M_mat: the n x n mass (inertia) matrix
%K_mat: the n x n stiffness matrix
function [M_mat,K_mat] = construct_2nd_order_matrices(string_params)
    
    n = string_params.n; %number of masses
    M = string_params.M; %total mass attached to the string    
    Tf = string_params.Tf; %tension in string
    dx = string_params.dx; %horizontal spacing between masses
    n = string_params.n;
    I_n = eye(n);
    M_left = circshift(I_n,[0,-1]);
    M_right = circshift(I_n,[0,1]);
    my_laplacian = M_left-2*I_n+M_right;

    my_laplacian(1,n) = 0;
    my_laplacian(n,1) = 0;

    K_mat = (Tf/dx) * my_laplacian;
    M_mat = (M/n) *I_n;

end

%INPUTS
%t: current time
%V: system state. V = [U;dUdt] where
% U and dUdt are n x 1 column vectors
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
function dVdt = string_rate_func02(t,V,string_params)
    n = string_params.n; %number of masses
    M = string_params.M; %total mass attached to the string
    Uf_func = string_params.Uf_func; %function describing motion of end point
    dUfdt_func = string_params.dUfdt_func; %time derivative of Uf

    Tf = string_params.Tf; %tension in string
    L = string_params.L; %length of string
    c = string_params.c; %damping coefficient
    dx = string_params.dx; %horizontal spacing between masses
    
    %unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);

   U_left = [0;U(1:end-1)];
   U_right = [U(2:end);Uf];

   dUdt_left = [0;dUdt(1:end-1)];
   dUdt_right = [dUdt(2:end);dUfdt];

   term1 = (Tf/dx)*(U_left-2*U+U_right);
   term2 = (c/dx)*(dUdt_left-2*dUdt+dUdt_right);

    
    %compute acceleration
    d2Udt2 = (term1+term2)/(M/n);
    %assemble state derivative
    dVdt = [dUdt;d2Udt2];
end 