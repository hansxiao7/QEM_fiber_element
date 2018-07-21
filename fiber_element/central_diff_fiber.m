function data = central_diff_fiber(delta_t, tend, Qe, rho, A, L, N, I, sections)

% [xn, wn] = lobatto_points(N);

nout = floor(tend/delta_t) + 1;

data.t = zeros(1, nout);
data.e = zeros(2*N + 2, nout);
data.ed = zeros(2*N + 2, nout);

% Initial conditions (t=0)
e = zeros(2*N +2,1);
ed = zeros(2*N + 2, 1);

% Store data at t0
data.t(1) = 0;
data.e(:,1) = e;
data.ed(:,1) = ed;

% Mass matrix (constant)
M = mass(rho, A, L, N) + rotation_mass(rho, I , L, N);
Minv = inv(M);
%left end fiexed
Minv(1:3, :) = zeros(3, 2*N+2);
Minv(:, 1:3) = zeros(2*N+2, 3);


old_e = e;

% Loop over all output times.
hw = waitbar(0,'Initializing waitbar...');

for j = 1:nout-1
    [Qint,sections] = internal_force_inelastic(e, L, N, sections);
    Qint(1:3, 1) = zeros(3,1);
    temp = e;
    e = delta_t^2 * Minv * (Qe - Qint) +2*e - old_e;
    e(1:3,1) = zeros(3,1);
    old_e = temp;
    data.e(:,j+1) = e;
    data.t(j+1) = data.t(j) + delta_t;
    t_str = sprintf('t = %.4f', j*delta_t);
    waitbar(j/(nout-1),hw,t_str);
end
close(hw);
end

