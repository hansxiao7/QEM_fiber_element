clc
clear

%choose the number of points along the beam
N = 13;
[xn, wn] = lobatto_points(N);

%initialization
e = zeros(2*N + 2, 1);

A = 0.1*0.1;
E = 2.07e11;
L = 2;
I = 0.1*0.1^3/12;
F= 3*E*I/L^2/4;

%external force vector
Qext = zeros(2*N+2, 1);
Qext(2*N+1,1) = -F;

Kf = zeros(2*N+2, 2*N+2);
Kl = zeros(2*N+2, 2*N+2);

Q = zeros(2*N+2, 1);
K = zeros(2*N+2, 2*N+2);


for j = 1:100
    Q = zeros(2*N+2, 1);
    K = zeros(2*N+2, 2*N+2);
    for i = 1:N
        x = xn(i);
        w = wn(i);
        Q = Q + w * (Q_axial(e, x, L, N, E, A)+Q_f(e, x, L, N, E, I));
        K = K + w * (k_axial(e, x, L, N, E, A)+k_f(e, x, L, N, E, I));
    end
    Q(1:3,1) = zeros(3,1);
    k22 = K(4:2*N+2, 4:2*N+2);
%     KK(4:2*N+2, 4:2*N+2, j) = k22;
    delta_e = k22\(Qext(4:2*N+2,1) - Q(4:2*N+2,1));
    e = e+ [0;0;0;delta_e];
end

    