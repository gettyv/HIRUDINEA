function [N, C, C_B, C_S, M, M_B, M_S] = sys_def(r, theta, n_v, b, phi, alpha)
%sys_def Creates the N, C, and M matrices for given system variables
%   Detailed explanation goes here

% Base structure of vertebrae
V_base = [
    r*cos(pi-theta/2),r*cos(pi-theta/2),r;
    r*sin(pi-theta/2),-r*sin(pi-theta/2),0;
    1,1,1;
];

% Construct N Matrix
N = [V_base zeros(3, 3*n_v-3)];

% Iterate through vertebrae as groups of three columns
for i = 3:3:3*(n_v-1)

    % Find center of previous vertebrae
    prev_ctr = [
        mean(N(1, i-2:i));
        mean(N(2, i-2:i));
        mean(N(3, i-2:i));
    ];

    % Create a new vertebrae by taking the base one and translating it
    starting_translation = lintr((2*r-b) + prev_ctr(1), prev_ctr(2));
    N(:,1+i) = starting_translation * V_base(:,1);
    N(:,2+i) = starting_translation * V_base(:,2);
    N(:,3+i) = starting_translation * V_base(:,3);

    % Rotate the new vertebrae around the center of the previous vertebrae
    phi_rot = rotate_abt_point(prev_ctr(1), prev_ctr(2), phi*(i/3));
    N(:,1+i) = phi_rot * N(:,1+i);
    N(:,2+i) = phi_rot * N(:,2+i);
    N(:,3+i) = phi_rot * N(:,3+i);

    % Find the new center of the vertebrae
    average_after_phi = [
        mean(N(1, i+1:i+3));
        mean(N(2, i+1:i+3));
        mean(N(3, i+1:i+3));
    ];

    % Rotate the vertebrae around its own center
    alpha_rot = rotate_abt_point(average_after_phi(1), average_after_phi(2), alpha*(i/3));
    N(:,1+i) = alpha_rot * N(:,1+i);
    N(:,2+i) = alpha_rot * N(:,2+i);
    N(:,3+i) = alpha_rot * N(:,3+i); 
end

% Construct C Matrix
E = eye(3*n_v);
D_B = zeros(3*n_v,3*n_v);
D_S = zeros(3*n_v,4*(n_v-1));
for i = 1:3:3*n_v
    D_B(:,i)   = E(:,(i+2)) - E(:,(i+1));
    D_B(:,i+1) = E(:,i) - E(:,(i+2));
    D_B(:,i+2) = E(:,(i+1)) - E(:,i);
end

connection_array = 1:3:3*(n_v-1);
D_index_array = 1:4:4*(n_v-1);

for i = 1:numel(connection_array)
    % TO FIX: Maybe keep current i independent from saving variable as they
    % already work well for accessing the E matrix
    j = D_index_array(i);
    k = connection_array(i);
    D_S(:,j)   = E(:,(k+4)) - E(:,(k+1));
    D_S(:,j+1) = E(:,(k+3)) - E(:,k);
    D_S(:,j+2) = E(:,(k+4)) - E(:,(k+2));
    D_S(:,j+3) = E(:,(k+3)) - E(:,(k+2));
end

D = [D_B,D_S];

C_B = D_B';
C_S = D_S';
C = D';

M = N * C';

M_S = N * D_S;
M_B = N * D_B;


end

% Helper Fn
function TR = lintr(tx, ty)
    TR = [1 0 tx; 0 1 ty; 0 0 1];
end

function RTR = rotate_abt_point(tx, ty, theta)
    ROT = rotz(theta);
    RTR = lintr(tx, ty) * ROT * lintr(-tx, -ty);
end

