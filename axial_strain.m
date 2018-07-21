function epislon = axial_strain(e, x, L, N)
% -1 <= x <= 1

Sx = shape_fun(N, x, 1);

B = Sx * e * 2/L;
pupx = B(1);
pvpx = B(2);

epislon = pupx + 1/2* (pupx^2 + pvpx^2);

end