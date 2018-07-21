clc
clear
%multielement

%choose the number of points along the beam
N = 7;
[xn, wn] = lobatto_points(N);

%initialization
e = zeros(16*N-5, 1);

A = 0.1*0.1;
E = 2.07e11;
L = 2;
I = 0.1*0.1^3/12;
F= 3*E*I/(16^2);

%external force vector
Qext = zeros(16*N-5, 1);
Qext(16*N-6,1) = -F;


for j = 1:10
    K_all = zeros(16*N-5, 16*N-5);
    Q_all = zeros(16*N-5, 1);
    for m = 1:8
        e_ele = e((m-1)*2*N + 2 -m : m*2*N+3-m, 1);
        Q = zeros(2*N+2, 1);
        K = zeros(2*N+2, 2*N+2);
        for i = 1:N
            x = xn(i);
            w = wn(i);
            Q = Q + w * (Q_axial(e_ele, x, L, N, E, A)+Q_f(e_ele, x, L, N, E, I));
            K = K + w * (k_axial(e_ele, x, L, N, E, A)+k_f(e_ele, x, L, N, E, I));
        end
        Q_all((m-1)*2*N + 2 -m : m*2*N+3-m, 1) =Q_all((m-1)*2*N + 2 -m : m*2*N+3-m, 1) +  Q;
        K_all((m-1)*2*N + 2 -m : m*2*N+3-m, (m-1)*2*N + 2 -m : m*2*N+3-m) = K_all((m-1)*2*N + 2 -m : m*2*N+3-m, (m-1)*2*N + 2 -m : m*2*N+3-m) + K;
    end
  
    Q_all(1:3,1) = zeros(3,1);
    k22 = K_all(4:16*N-5, 4:16*N-5);
%     KK(4:16*N-5, 4:16*N-5, j) = k22;
    delta_e = k22\(Qext(4:16*N-5,1) - Q_all(4:16*N-5,1));
    e = e+ [0;0;0;delta_e];
end

    