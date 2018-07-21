function k = k_axial(e, x, L, N, E, A)

% -1 <= x <= 1

Sx = shape_fun(N,x, 1);

Sx1 = Sx(1,:);
Sx2 = Sx(2,:);

k = E*A*L/2*((4/L^2*(Sx1')*Sx1 + 4/L^2*(Sx2')*Sx2) * axial_strain(e, x, L, N)...
    +(Sx1' * 2/L + 4/L^2*(Sx1')*Sx1 * e + 4/L^2*(Sx2')*Sx2 * e) * (Sx1 + e'* (Sx1' * Sx1 + Sx2' * Sx2)));

end
