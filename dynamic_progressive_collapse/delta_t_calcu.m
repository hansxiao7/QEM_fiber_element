clc
clear
%calculate the delta_t
N = 7;
[xn, wn] = lobatto_points(N);

e = zeros(2*N + 2, 1);
A = 26*34;
E = 4110;
L = 206;
I = 34*26^3/12;
rho = 2.17e-7;

K_origin = zeros(2*N+2, 2*N+2);
for i = 1:N
    x = xn(i);
    w = wn(i);
    K_origin = K_origin + w * (k_axial(e, x, L, N, E, A)+k_f(e, x, L, N, E, I));
end

M = mass(rho, A, L, N);
%additional mass from column
M(2*N+1, 2*N+1) = M(2*N+1, 2*N+1) + 0.0128;
M(2*N, 2*N) = M(2*N, 2*N) + 0.0128;
Minv = zeros(2*N+2, 2*N+2);
for i = 1:2*N+2
    if M(i, i)~= 0 
        Minv(i,i) = 1/M(i,i);
    end
end

omega = eig(Minv * K_origin);
omega = sqrt(max(omega));

delta_t = 2 / omega * 0.8;
