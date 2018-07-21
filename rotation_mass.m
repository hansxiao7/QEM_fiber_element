function M = rotation_mass(rho, I, L, N)

[xn, wn] = lobatto_points(N);

M = zeros(2*N+2,2*N+2);

for i = 1:N
    w = wn(i);
    x = xn(i);
    Sx = shape_fun(N, x, 1) * 2/L;
    Sx(1,:) = zeros(1,2*N + 2);
    M = M + w* (Sx') *Sx;
end

M = M * rho * I * L /2;


end