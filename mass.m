function M = mass(rho, A, L, N)

[xn, wn] = lobatto_points(N);

M = zeros(2*N+2,2*N+2);

for i = 1:N
    w = wn(i);
    x = xn(i);
    S = shape_fun(N, x, 0);
    M = M + w* (S') *S;
end

M = M * rho * A * L /2;


end
