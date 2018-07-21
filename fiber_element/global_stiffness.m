function k = global_stiffness(e, x, L, N, section)

k_local = local_stiffness(section);

Sx = shape_fun(N, x, 1) * 2/ L;
Sxx = shape_fun(N, x, 2) * 4/ L^2;

S1x = Sx(1, :);
S2x = Sx(2, :);
S2xx = Sxx(2, :);

A = [S1x'+ S1x'*S1x*e + S2x'*S2x*e, S2xx'];

k = real(L / 2 * A*k_local*A');

end