%INPUTS
%t: current time
%V: system state. V = [U;dUdt] where
%   U and dUdt are n x 1 column vectors
%string_params: a struct containing the system parameters describing the string
%   string_params.n: number of masses
%   string_params.M: total mass attached to the string
%   string_params.Uf_func: function describing motion of end point
%   string_params.dUfdt_func: time derivative of Uf
%   string_params.Tf: %tension in string
%   string_params.L: %length of string
%   string_params.c: %damping coefficient
%   string_params.dx: %horizontal spacing between masses
function dVdt = string_rate_func01(t,V,string_params)
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

    %compute acceleration
    %% for loop way
    % d2Udt2 = zeros(1, n);
    % for i = 1:n
    %     if i == 1
    %         U_im1 = [0; 0];
    %     else
    %         U_im1 = [U(i-1); dUdt(i-1)];
    %     end
    %     if i == n
    %         U_ip1 = [Uf; dUfdt];
    %     else
    %         U_ip1 = [U(i+1); dUdt(i+1)];
    %     end
    %     d2Udt2(i) = (n / M) * ((Tf / dx) * (U_im1(1) - 2*U(i) + U_ip1(1)) + (c / dx) * (U_im1(2) - 2*dUdt(i) + U_ip1(2)));
    % end

    %% laplacian way
    %construct the nxn discrete laplacian matrix
    laplacian = circshift(eye(n),-1,2) + -2*eye(n) + circshift(eye(n),1,2);
    laplacian(1,end) = laplacian(1,end)-1; %delete unwanted 1 in top right corner
    laplacian(end,1) = laplacian(end,1)-1; %delete unwanted 1 in bottom right corner

    B1 = zeros(n, 1);
    B1(end) = Uf;
    B2 = zeros(n, 1);
    B2(end) = dUfdt;

    d2Udt2 = ((n / M) * ((Tf / dx) * (sum(laplacian .* U, 2) + B1) + (c / dx) * (sum(laplacian .* dUdt, 2) + B2)))';

    %assemble state derivative
    dVdt = [dUdt;d2Udt2];
end