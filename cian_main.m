close all;clc; clear;
num_masses = 29;
total_mass = 4;
tension_force = 10;   % N
string_length = 10;  % m
amplitude_Uf = 0.01;

damping_coeff = 0.15; %string_length / tension_force
dx = string_length / (num_masses + 1);

string_params = struct();
string_params.n = num_masses;
string_params.M = total_mass;
string_params.Tf = tension_force;
string_params.L = string_length;
string_params.c = damping_coeff;
string_params.dx = dx;

%Get nat freqs
[M_mat, K_mat] = construct_2nd_order_matrices(string_params);
[Ur_mat,lambda_mat] = eig(K_mat,M_mat);

nat_freqs = sqrt(diag(lambda_mat));
string_params.Ur = Ur_mat;
string_params.nat_freqs = nat_freqs;

omega_Uf = nat_freqs(1);

% boundary forcing at right end
Uf_func    = @(t) amplitude_Uf * cos(omega_Uf * t);
dUfdt_func = @(t) -omega_Uf * amplitude_Uf * sin(omega_Uf * t);
string_params.Uf_func = Uf_func;
string_params.dUfdt_func = dUfdt_func;

%Sim
nat_freq_string_simulation(string_params)
