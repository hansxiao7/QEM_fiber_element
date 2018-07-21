function SS = shape_fun(N, x, der)

S = zeros(2, 2*N + 2);
Sx = zeros(2, 2*N + 2);
Sxx = zeros(2, 2*N + 2);

C = lagrange_interp(N,x);
Cx = lagrange_interp_x(N,x);
Cxx = lagrange_interp_xx(N,x);

for i = 1:N
    S(1, 2*i) = C(i);
    Sx(1, 2*i) = Cx(i);
    Sxx(1, 2*i) = Cxx(i);
end

[h11, hn1, hi0] = H(N, x);
[h11_x, hn1_x, hi0_x] = H_x(N, x);
[h11_xx, hn1_xx, hi0_xx] = H_xx(N, x);

S(2,1) = h11;
Sx(2,1) = h11_x;
Sxx(2,1) = h11_xx;
S(2, 2*N +2) = hn1;
Sx(2, 2*N +2) = hn1_x;
Sxx(2, 2*N +2) = hn1_xx;

for j = 1:N
    S(2, 2*j+1) = hi0(j);
    Sx(2, 2*j+1) = hi0_x(j);
    Sxx(2, 2*j+1) = hi0_xx(j);
end

switch der
    case 0
        SS = S;
    case 1
        SS = Sx;
    case 2
        SS = Sxx;
end

SS(1,1) = SS(1,2);
SS(1,2) = 0;
SS(2,2) = SS(2,3);
SS(2,3) = SS(2,1);
SS(2,1) = 0;


end
