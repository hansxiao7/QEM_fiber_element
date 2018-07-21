function k = k_axial_elastic(x, L, N, E, A)

Sx = shape_fun(N, x, 1);

Sx1 = Sx(1,:);


k = E*A*L/2 * (Sx1'*Sx1 * 4/L^2);

end