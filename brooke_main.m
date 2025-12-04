string_params = struct();

t = 4;
V = [
    .03, -.02
    .01, 0
    -.05, .02
    ];
string_params.n = 3;
string_params.M = 2;
string_params.Uf_func = @(t)3*t;
string_params.dUfdt_func = @(t)3;
string_params.Tf = 9;
string_params.L = 12;
string_params.c = string_params.L / string_params.Tf;
string_params.dx = string_params.L / (string_params.n+1);

string_rate_func01(t, V, string_params)
