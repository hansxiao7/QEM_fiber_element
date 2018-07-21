function phi = curvature(e, x, L, N)

Sxx = shape_fun(N, x, 2)*4/L^2;

phi = Sxx(2,:) * e;