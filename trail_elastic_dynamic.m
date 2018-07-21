clc
clear

%choose the number of points along the beam
N = 13;
[xn, wn] = lobatto_points(N);

%initialization
e = zeros(2*N + 2, 1);

A = 26*34;
E = 4110;
L = 206;
I = 34*26^3/12;
F= 3*E*I/L^2/8;
rho = 84e-6;

%external force vector
Qext = zeros(2*N+2, 1);
Qext(2*N+1,1) = -F;

% K_origin = zeros(2*N+2, 2*N+2);
% for i = 1:N
%     x = xn(i);
%     w = wn(i);
%     K_origin = K_origin + w * (k_axial(e, x, L, N, E, A)+k_f(e, x, L, N, E, I));
% end
% 
% M = mass(rho, A, L, N) + rotation_mass(rho, I , L, N);
% Minv = inv(M);
% omega = eig(Minv * K_origin);
% omega = sqrt(max(omega));

% delta_t = 2 / omega * 0.8*100;
delta_t = 0.001;
tend = 2;

% data = central_diff_elastic(delta_t, tend, Qext, rho, A, L, N, E, I);
data = implicit_elastic(delta_t, tend, Qext, rho, A, L, N, E, I);