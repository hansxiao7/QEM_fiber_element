function Q = gravity_force(rho_g, N, L, A)

[xn,wn] = lobatto_points(N);

Q = zeros(2*N+2, 1);

for i = 1:N
    x = xn(i);
    w = wn(i);
    S = shape_fun(N, x, 0);
    Q = Q + w*L/2* S' * [0; -rho_g *A];
end

end