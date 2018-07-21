clc
clear

N= 5;
[xn, wn] = lobatto_points(N);

rho_g = 84e-6;
rho_g = rho_g + 80.88e-6;
A = 34*26;
L = 206;

I = 34*26^3/12;

Qg = gravity_force(rho_g, N, L, A);

%fixed-fixed beam
Qg(1:3,1) = zeros(3,1);
Qg(2*N:2*N+2, 1) = zeros(3,1);

sections = sections_init(N);
%initialization
e = zeros(2*N + 2, 1);

result = zeros(2*N+2, 10);

for m = 1:10
    Qext = m*0.1*Qg;
        for i = 1:1000
            [Q,sections] = internal_force_inelastic(e, L, N, sections);
 
            K = zeros(2*N+2, 2*N+2);
            for j = 1:N
                x = xn(j);
                w = wn(j);
                K = K + w * global_stiffness(e, x, L, N, sections(j));
            end
            Q(1:3,1) = zeros(3,1);
            Q(2*N:2*N+2,1) = zeros(3,1);
            k22 = K(4:2*N-1, 4:2*N-1);
            delta_e = k22\(Qext(4:2*N-1,1) - Q(4:2*N-1,1));
            e = e+ [0;0;0;delta_e;0;0;0];
            if norm(delta_e)<1e-16
                break;
            end
        end
    result(:, m) = e;
end
