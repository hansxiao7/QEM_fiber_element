function k = k_f(e, x, L, N, E, I)

Sxx = shape_fun(N, x, 2);
Sxx2 = Sxx(2,:);

f = constant_f(e, x,L,  N);

k = E*I*L/2 * 16/ (L^4) *(Sxx2')*Sxx2 ;

end