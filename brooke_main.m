close all; clc; clear;
string_params = struct();

string_params.n = 4;    % number of masses
string_params.M = 1;    % total mass (kg)
string_params.Tf = 10000;  % tension force (N)
string_params.L = 10;   % string length (m)
string_params.c = string_params.L / string_params.Tf;   % damping coefficient
string_params.dx = string_params.L / (string_params.n+1);   % distance between masses

amplitude_Uf = 0.01;

% get natural frequencies
[M_mat, K_mat] = construct_2nd_order_matrices(string_params);
[U_r, omega_r] = find_mode_shape_and_resonant_frequencies(M_mat, K_mat);

for i = 1:length(U_r)
    omega_Uf = omega_r(i);
    string_params.omega_Uf = omega_Uf;
    string_params.U_r = U_r(:, i) * 0.1;
    string_params.Uf_func = @(t) amplitude_Uf * cos(omega_Uf * t);
    string_params.dUfdt_func = @(t) -omega_Uf * amplitude_Uf * sin(omega_Uf * t);
    string_simulation(string_params)
end
