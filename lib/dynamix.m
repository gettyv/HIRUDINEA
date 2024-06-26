function out = dynamix(t, state, W, psi, C_S, Y, gamma, J, m, t_mu, mu, ol_f, ol_phase, ol_mag, ol_fn)
%dynamix This ODE solver function calculates the accelerations of the
%tensegrity bars
%   This takes in the state vector at each time step along with some other
%   parameters to determine what the updated acceleration of the bars would
%   be, then returns those accelerations.
%   ol_fn is a function that can be called in the form f(f*t + phi)


% Rotation matrices
% Hardcoded for computation time
rotz_120 = [[-0.5, -0.866, 0],
            [0.8660, -0.5, 0],
            [0, 0, 1]];
rotz_n120 = [[-0.5, 0.866, 0],
            [-0.8660, -0.5, 0],
            [0, 0, 1]];
rotz_90 = [[0, -1, 0],
            [1, 0, 0],
            [0, 0, 1]];



% Extract # of vertebrae from W size
n_v = size(W, 2) / 3;

% Reshape state vector to recover structure
state = reshape(state, 3, 4 * n_v);

% The Q matrix needs to be reconstructed from the state
R = zeros([3, n_v]);
B = zeros([3, n_v * 3]);
R_dot = zeros(size(R));
B_dot = zeros(size(B));
for v = 1:n_v
    R(:,v) = state(:,v);

    R_dot(:,v) = state(:, n_v*2 + v);

    B(:,(3*v)-2) = state(:,n_v + v);
    B(:,(3*v)-1) = rotz_120 * state(:,n_v + v);
    B(:,(3*v)-0) = rotz_n120 * state(:,n_v + v);

    B_dot(:,(3*v)-2) = state(:,n_v*3 + v);
    B_dot(:,(3*v)-1) = rotz_120 * state(:,n_v*3 + v);
    B_dot(:,(3*v)-0) = rotz_n120 * state(:,n_v*3 + v);
end
Q = [R B];
Q_dot = [R_dot B_dot];

ol_map = ones(1, 4 * (n_v-1));
for v = 1:(n_v-1)
    % Calculate open loop control signal and apply to map
    ol_cmd = ol_mag * ol_fn(ol_f * t + (v-1)*ol_phase);
    ol_map(4*(v-1) + 1) = 1+ol_cmd;
    ol_map(4*(v-1) + 2) = 1-ol_cmd;
end

% Apply open loop command to gamma values
gamma_array = diag(gamma)';
gamma_ol = diag(ol_map .* gamma_array);

% Now the forces on each of the points can be calculated
F_Q = (W - (Q*transpose(psi)+Y)*transpose(C_S)*gamma_ol*C_S)*psi;

f_r = F_Q(:, 1:n_v);
f_r_friction = zeros(size(f_r));
f_b = F_Q(:, n_v+1:end);

% Calculate forces and accelerations of r, b
tau = zeros([3, n_v]);
tau_friction = zeros(size(tau));
r_ddot_m = zeros([3, n_v]);
b_ddot_m = zeros([3, n_v]);

for v = 1:n_v
    index = (3*v)-2;
    tau(:,v) = cross(B(:,index), f_b(:,index)) + cross(B(:,index+1), f_b(:,index+1)) + cross(B(:,index+2), f_b(:,index+2));

    % Apply friction to sliding and rotating actions
    %tau_friction(:,v) = frix(tau(:,v), 0);
    %f_r_friction(:,v) = frix(f_r(:,v), 0);

    cross_p = cross(B(:,index), tau(:,v));
    b = B(:,index);
    b_dot = B_dot(:,index);    

    tangent_unit_vector = B(:,index+2) / norm(B(:,index+2));
    perp_unit_vector = rotz_90 * B(:,index+2) / norm(B(:,index+2));

    tangent_mu = 0.1;
    perp_mu = 1.6;

    r_ddot_m(:,v) = (f_r(:,v) - tangent_mu * dot(R_dot(:,v), tangent_unit_vector) * tangent_unit_vector - perp_mu * dot(R_dot(:,v), perp_unit_vector) * perp_unit_vector)/ m;
    b_ddot_m(:,v) = -(cross_p/(J*(norm(b))^2))-b*(norm(b_dot)^2)/(norm(b)^2) - t_mu * b_dot;
end

accels = zeros(3, 2*n_v);
vels = zeros(size(accels));
for v = 1:n_v
    index = (3*v)-2;
    accels(:,v) =       r_ddot_m(:,v);
    accels(:,n_v+v) =   b_ddot_m(:,v);
    vels(:,v) =         R_dot(:,v);
    vels(:,n_v+v) =     B_dot(:,index);
end

out = reshape([vels accels], [], 1);
end

