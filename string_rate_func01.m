function dVdt = string_rate_func01(t, V, string_params)

    n  = string_params.n;
    M  = string_params.M;
    Uf = string_params.Uf_func(t);
    dUfdt = string_params.dUfdt_func(t);
    Tf = string_params.Tf;
    c  = string_params.c;
    dx = string_params.dx;

    % unpack state vector
    U     = V(1:n);
    dUdt  = V(n+1:2*n);

    % Laplacian
    Q = my_laplacian(n);

    % boundary forcing at right end
    B1 = zeros(n,1); 
    B1(end) = Uf;
    B2 = zeros(n,1);
    B2(end) = dUfdt;

    a1 = (Tf/dx) * (Q*U + B1);
    a2 = (c/dx)  * (Q*dUdt + B2);

    %deriv
    d2Udt2 = (n/M) * (a1 + a2);
    dVdt = [dUdt; d2Udt2];
end
