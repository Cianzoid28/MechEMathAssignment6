string_params = struct();

t = 1;
V = [
    -.02, 0
    .02, 0
    -.02, 0
    ];
string_params.n = 3;
string_params.M = 2;
string_params.Uf_func = @(t)0.01*sin(t);
string_params.dUfdt_func = @(t)0.01*cos(t);
string_params.Tf = 9;
string_params.L = 12;
string_params.c = string_params.L / string_params.Tf;
string_params.dx = string_params.L / (string_params.n+1);

string_rate_func01(t, V, string_params)

[M_mat, K_mat] = construct_2nd_order_matrices(string_params);
[U_r, lambda] = find_mode_shape_and_resonant_frequencies(M_mat, K_mat)
