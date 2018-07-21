function data = implicit_elastic(delta_t, tend, Qe, rho, A, L, N, E, I)

[xn, wn] = lobatto_points(N);

nout = floor(tend/delta_t) + 1;

data.t = zeros(1, nout);
data.e = zeros(2*N + 2, nout);
data.ed = zeros(2*N + 2, nout);
data.a = zeros(2*N + 2, nout);

% Initial conditions (t=0)
e = zeros(2*N +2,1);
ed = zeros(2*N + 2, 1);
a = zeros(2*N + 2, 1);

% Store data at t0
data.t(1) = 0;
data.e(:,1) = e;
data.ed(:,1) = ed;
data.a(:,1) = a;

% Mass matrix (constant)
M = mass(rho, A, L, N) + rotation_mass(rho, I , L, N);

M(1:3, :) = zeros(3, 2*N+2);
M(:, 1:3) = zeros(2*N+2, 3);

% Loop over all output times.
hw = waitbar(0,'Initializing waitbar...');

for i = 1:nout -1 
    
    delta_e = zeros(2*N-1, 1);    
    old_e = e;
    
    for m = 1:5000
        Qt = zeros(2*N + 2, 1);
        e = e + [0;0;0;delta_e];
        for k = 1:N
            x = xn(k);
            w = wn(k);
            Qt = Qt + w * (Q_axial(e, x, L, N, E, A)+Q_f(e, x, L, N, E, I));
        end
        Qt(1:3, 1) = zeros(3,1);
        Kt = zeros(2*N + 2, 2*N + 2);
        for j = 1:N
            x = xn(j);
            w = wn(j);
            Kt = Kt + w * (k_axial(e, x, L, N, E, A)+k_f(e, x, L, N, E, I));
        end
        Keff = Kt + 4/delta_t^2 * M;
        k22 = Keff(4:2*N+2, 4:2*N+2);
        k22inv = real(inv(k22));
        Feff = Qe - Qt - M * ( 4/delta_t^2*(e - old_e) - 4/ delta_t*ed - a);
        delta_e = (k22inv)*Feff(4:2*N+2,1);
    end
    
    e = e + [0;0;0;delta_e];
    a = 4/delta_t^2*(e - old_e) - 4/delta_t*ed - a;
    ed = -ed + 2/delta_t*(e - old_e);
    
    data.e(:,i+1) = e;
    data.ed(:, i+1) = ed;
    data.a(:,i+1) = a;
    data.t(i+1) = data.t(i) + delta_t;
    t_str = sprintf('t = %.4f', i*delta_t);
    waitbar(i/(nout-1),hw,t_str);
end
close(hw);
end
    
    


