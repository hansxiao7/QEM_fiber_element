clc
clear

N = 5;
[xn, wn] = lobatto_points(N);

%initialization
e = zeros(2*N + 2, 1);

A = 26*34;
L = 206;
E = 4110;
I = 34*26^3/12;
rho_g = 8.4e-5;
Qg = gravity_force(rho_g, N, L, A);

%progressive collapse beam condition
Qg(1:3,1) = zeros(3,1);
Qg(2*N, 1) = 0;
Qg(2*N+2, 1) = 0;

%external force from column
Q_column = zeros(2*N+2, 1);
Q_column(2*N + 1, 1) = -4.95;

%total external force = beam gravity + column gravity
Qext = Qg + Q_column;

sections = sections_init(N);

dd =[];
% for k=1:60
%     F= k;
% 
%     %external force vector
%     Qext = zeros(2*N+2, 1);
%     Qext(2*N+1,1) = -F;
% 

for i = 1:1000
    [Q,sections] = internal_force_inelastic(e, L, N, sections);
    if isreal(Q)==0
        error('aaa');
    end
    K = zeros(2*N+2, 2*N+2);
    for j = 1:N
        x = xn(j);
        w = wn(j);
        K = K + w * global_stiffness(e, x, L, N, sections(j));
    end
    Q(1:3,1) = zeros(3,1);
    Q(2*N, 1) = 0;
    Q(2*N+2, 1) = 0;
    k22 = K(4:2*N-1, 4:2*N-1);
    k23 = K(2*N+1, 4:2*N-1);
    k32 = K(4:2*N-1, 2*N+1);
    k33 = K(2*N+1, 2*N+1);
    KK = [k22, k32;k23, k33];

%     KK(4:2*N+2, 4:2*N+2, j) = k22;
    delta_e = KK\([Qext(4:2*N-1,1); Qext(2*N+1,1)] - [Q(4:2*N-1,1); Q(2*N+1,1)]);
    e = e+ [0;0;0;delta_e(1:max(size(delta_e))-1);0;delta_e(max(size(delta_e)));0];
end
   dd = [dd,e(2*N+1,1)]; 
% end