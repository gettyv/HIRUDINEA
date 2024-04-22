%r b1 b2 b3 method
clc;
clear all;
close all;

% length of the radius of circle 
r = 1;
%syms r

% distance between two spines % distance between two circles
b = 1.5;
%syms b

% angle between two bar (radian)
theta = 2*pi/3;
%syms theta

% angle of rotation
phi = 0*pi/180;
%phi=0;
%syms phi

% rotation matrix
rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];

% define the nodes (currently 6 nodes, 10 members)
n_1 = [r*cos(pi-theta/2);r*sin(pi-theta/2)];
n_2 = [r*cos(pi-theta/2);-r*sin(pi-theta/2)];
n_3 = [r;0];
n_4 = rot*[r*cos(pi-theta/2)+(2*r-b);r*sin(pi-theta/2)];
n_5 = rot*[r*cos(pi-theta/2)+(2*r-b);-r*sin(pi-theta/2)];
n_6 = rot*[3*r-b;0];
figure(1);
plot(n_1(1),n_1(2),'-s',n_2(1),n_2(2),'-o',n_3(1),n_3(2),'-p',n_4(1),n_4(2),'-s',n_5(1),n_5(2),'-o',n_6(1),n_6(2),'-p'); hold on;
plot([n_3(1),n_2(1)],[n_3(2),n_2(2)],'k-', 'LineWidth', 1); hold on;
plot([n_1(1),n_3(1)],[n_1(2),n_3(2)],'k-', 'LineWidth', 1); hold on;
plot([n_2(1),n_1(1)],[n_2(2),n_1(2)],'k-', 'LineWidth', 1); hold on;
plot([n_6(1),n_5(1)],[n_6(2),n_5(2)],'b-', 'LineWidth', 1); hold on;
plot([n_4(1),n_6(1)],[n_4(2),n_6(2)],'b-', 'LineWidth', 1); hold on;
plot([n_5(1),n_4(1)],[n_5(2),n_4(2)],'b-', 'LineWidth', 1); hold on;
axis square
axis([-3 3 -3 3])
% define the members: 6 bars and 4 cables
s_1 = n_5 - n_2;
s_2 = n_4 - n_1;
s_3 = n_5 - n_3;
s_4 = n_4 - n_3;

b_1 = n_3 - n_2;
b_2 = n_1 - n_3;
b_3 = n_2 - n_1;
b_4 = n_6 - n_5;
b_5 = n_4 - n_6;
b_6 = n_5 - n_4;

%consturct D matrix
E = eye(6);
D_B = zeros(6,6);
D_S = zeros(6,4);
D_B(:,1) = E(:,3)- E(:,2);
D_B(:,2) = E(:,1)- E(:,3);
D_B(:,3) = E(:,2)- E(:,1);
D_B(:,4) = E(:,6)- E(:,5);
D_B(:,5) = E(:,4)- E(:,6);
D_B(:,6) = E(:,5)- E(:,4);

D_S(:,1) = E(:,5)- E(:,2);
D_S(:,2) = E(:,4)- E(:,1);
D_S(:,3) = E(:,5)- E(:,3);
D_S(:,4) = E(:,4)- E(:,3);
D = [D_B,D_S];

%construct C matrix
C_B = transpose(D_B);
C_S = transpose(D_S);
C = transpose(D);
% eliminate one row of C, not necessary for null space
%C = C(:,1:5);

%construct N,M,B,S
N = [n_1,n_2,n_3,n_4,n_5,n_6];
B = [b_1,b_2,b_3,b_4,b_5,b_6];
S = [s_1,s_2,s_3,s_4];
M = [B,S];

% equilibrium
syms lambda_1 lambda_2 lambda_3 lambda_4 lambda_5 lambda_6;
LAMBDA = transpose(repmat([lambda_1; lambda_2; lambda_3; lambda_4; lambda_5; lambda_6],[1 2]));


syms gamma_1 gamma_2 gamma_3 gamma_4;
GAMMA = transpose(repmat([gamma_1; gamma_2; gamma_3; gamma_4],[1 2]));

F = [-B.*LAMBDA, S.*GAMMA]*C;
%solve(F(1,1)==0,F(2,1)==0,F(1,2)==0,F(2,2)==0,F(1,3)==0,F(2,3)==0,F(1,4)==0,F(2,4)==0,F(1,5)==0,F(2,5)==0,F(1,6)==0,F(2,6)==0)

% solve F = 0;
vars = [lambda_1,lambda_2,lambda_3,lambda_4,lambda_5,lambda_6,gamma_1,gamma_2,gamma_3,gamma_4];
[A,b0] = equationsToMatrix(F,vars);

rankA = rank(A);
redA = rref(A);

% Reduced row echelon form 
% Au = [A b0];
% R = rref(Au); %linsolve

nullA = null(A)
lambda_1_soln = cast(nullA(1),'double');
lambda_2_soln = cast(nullA(2),'double');
lambda_3_soln = cast(nullA(3),'double');
lambda_4_soln = cast(nullA(4),'double');
lambda_5_soln = cast(nullA(5),'double');
lambda_6_soln = cast(nullA(6),'double');
gamma_1_soln = cast(nullA(7),'double');
gamma_2_soln = cast(nullA(8),'double');
gamma_3_soln = cast(nullA(9),'double');
gamma_4_soln = cast(nullA(10),'double');


%% dynamics
%% initial values for Q: R and B

% angle of rotation
phi = 10*pi/180;
% rotation matrix
rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];
% rotate second spine to perturb it from equilibrium found above
n_4_rot = rot*[r*cos(pi-theta/2)+(2*r-b);r*sin(pi-theta/2)];
n_5_rot = rot*[r*cos(pi-theta/2)+(2*r-b);-r*sin(pi-theta/2)];
n_6_rot = rot*[3*r-b;0];
figure(1)
plot([n_6_rot(1),n_5_rot(1)],[n_6_rot(2),n_5_rot(2)],'b--', 'LineWidth', 1); hold on;
plot([n_4_rot(1),n_6_rot(1)],[n_4_rot(2),n_6_rot(2)],'b--', 'LineWidth', 1); hold on;
plot([n_5_rot(1),n_4_rot(1)],[n_5_rot(2),n_4_rot(2)],'b--', 'LineWidth', 1); hold on;
r1 = (n_1+n_2+n_3)/3; %r is center of mass, r1 is for the first spine
r2 = (n_4_rot+n_5_rot+n_6_rot)/3; %r2 for the second spine

% rotation matrix
% rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];
rotp120 = [cos(2*pi/3) -sin(2*pi/3); sin(2*pi/3) cos(2*pi/3)];
rotn120 = [cos(-2*pi/3) -sin(-2*pi/3); sin(-2*pi/3) cos(-2*pi/3)];
b1 = (n_1 - (n_1+n_2+n_3)/3)/norm(n_1 - (n_1+n_2+n_3)/3);
b2 = rotp120*b1;
b3 = rotn120*b1;
b4 = (n_4_rot - (n_4_rot+n_5_rot+n_6_rot)/3)/norm(n_4_rot - (n_4_rot+n_5_rot+n_6_rot)/3);
b5 = rotp120*b4;
b6 = rotn120*b4;

el = norm(n_1 - (n_1+n_2+n_3)/3); %length from centroid to node(all same)

%Q = [r1, b1 ,b2, b3, r2, b4, b5, b6];
Q = [r1,r2,b1,b2,b3,b4,b5,b6];

%{
%% psi and Y
% psi_t = [1, 1, 1, 0, 0, 0;
%         l, 0, 0, 0, 0, 0;
%         0, l, 0, 0, 0, 0;
%         0, 0, l, 0, 0, 0;
%         0, 0, 0, 1, 1, 1;
%         0, 0, 0, l, 0, 0;
%         0, 0, 0, 0, l, 0;
%         0, 0, 0, 0, 0, l];
psi_t = [1, 1, 1, 0, 0, 0;
        0, 0, 0, 1, 1, 1;    
        l, 0, 0, 0, 0, 0;
        0, l, 0, 0, 0, 0;
        0, 0, l, 0, 0, 0;
        0, 0, 0, l, 0, 0;
        0, 0, 0, 0, l, 0;
        0, 0, 0, 0, 0, l];
psi = transpose(psi_t);        
Y = zeros(2,6);
%}

%% fix the first spine
psi_t = [0, 0, 0, 0, 0, 0;
        0, 0, 0, 1, 1, 1;    
        0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0;
        0, 0, 0, el, 0, 0;
        0, 0, 0, 0, el, 0;
        0, 0, 0, 0, 0, el];
psi = transpose(psi_t);        
Y = [n_1,n_2,n_3,zeros(2,3)];

N
(Q * psi_t + Y)

error = N - (Q * psi_t + Y)

%% W is external force: all zero for now
W = zeros(2,6);

%% Gamma
% lastA =[
%  
%                                                                                                                                                     -((r*cos(theta/2) - b*cos(phi) + 2*r*cos(phi) - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi))*(b^2*sin(phi)^3 + 4*r^2*sin(phi)^3 - 3*r^2*cos(theta/2)^2*sin(phi)^3 + r^2*cos(theta/2)^3*sin(phi)^3 - 2*r^2*sin(theta/2)^2*sin(phi)^3 + b^2*cos(phi)^2*sin(phi) + 4*r^2*cos(phi)^2*sin(phi) + r^2*sin(theta/2)*cos(phi) - 4*b*r*sin(phi)^3 + b^2*cos(theta/2)*sin(phi)^3 + b^2*sin(theta/2)*cos(phi)^3 - 5*r^2*sin(theta/2)*cos(phi)^2 + 6*r^2*sin(theta/2)*cos(phi)^3 - 2*r^2*sin(theta/2)*sin(phi)^2 - 2*r^2*cos(phi)*sin(phi) + 2*b*r*cos(theta/2)^2*sin(phi)^3 + b*r*sin(theta/2)^2*sin(phi)^3 - r^2*cos(theta/2)*cos(phi)*sin(phi) - 4*b*r*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)^2*sin(phi)^3 + b^2*cos(theta/2)*cos(phi)^2*sin(phi) + b^2*sin(theta/2)*cos(phi)*sin(phi)^2 + r^2*cos(theta/2)^2*cos(phi)*sin(phi) + 6*r^2*sin(theta/2)*cos(phi)*sin(phi)^2 + r^2*sin(theta/2)^2*cos(phi)*sin(phi) + 2*b*r*sin(theta/2)*cos(phi)^2 - 2*b*r*cos(theta/2)*sin(phi)^3 - 5*b*r*sin(theta/2)*cos(phi)^3 + b*r*sin(theta/2)*sin(phi)^2 - 3*r^2*cos(theta/2)^2*cos(phi)^2*sin(phi) + r^2*cos(theta/2)^3*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)*cos(phi)^2 - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)^3 + b*r*cos(phi)*sin(phi) - 2*r^2*sin(theta/2)^2*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)*sin(phi)^2 + b*r*cos(theta/2)*cos(phi)*sin(phi) + r^2*cos(theta/2)*sin(theta/2)^2*cos(phi)^2*sin(phi) - 2*b*r*cos(theta/2)*cos(phi)^2*sin(phi) - 5*b*r*sin(theta/2)*cos(phi)*sin(phi)^2 - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2 + 2*b*r*cos(theta/2)^2*cos(phi)^2*sin(phi) + b*r*cos(theta/2)*sin(theta/2)*cos(phi)^3 + b*r*sin(theta/2)^2*cos(phi)^2*sin(phi) + b*r*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2))/(r*(cos(theta/2) + 1)*(b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2 - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi))*(b*sin(phi) - r*sin(theta/2) - 2*r*sin(phi) + r*cos(theta/2)^2*sin(phi) + r*sin(theta/2)^2*sin(phi) + b*cos(theta/2)*sin(phi) - b*sin(theta/2)*cos(phi) - r*cos(theta/2)*sin(phi) + 3*r*sin(theta/2)*cos(phi)))
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           -(b^2*cos(phi)^3 - 2*r^2*cos(phi)^2 + 4*r^2*cos(phi)^3 - r^2*cos(theta/2)^2*sin(phi)^2 + b^2*cos(phi)*sin(phi)^2 - r^2*cos(theta/2)*cos(phi) + 4*r^2*cos(phi)*sin(phi)^2 + b*r*cos(phi)^2 - 4*b*r*cos(phi)^3 + 3*r^2*cos(theta/2)*cos(phi)^2 - 4*r^2*cos(theta/2)*cos(phi)^3 + 2*r^2*cos(theta/2)*sin(phi)^2 - 2*r^2*sin(theta/2)*sin(phi)^3 - r^2*cos(theta/2)^2*cos(phi)^2 + r^2*cos(theta/2)^2*cos(phi)^3 + r^2*sin(theta/2)*cos(phi)*sin(phi) - 4*b*r*cos(phi)*sin(phi)^2 - 4*r^2*cos(theta/2)*cos(phi)*sin(phi)^2 - 2*r^2*sin(theta/2)*cos(phi)^2*sin(phi) - b*r*cos(theta/2)*cos(phi)^2 + 2*b*r*cos(theta/2)*cos(phi)^3 - b*r*cos(theta/2)*sin(phi)^2 + b*r*sin(theta/2)*sin(phi)^3 + r^2*cos(theta/2)^2*cos(phi)*sin(phi)^2 + r^2*cos(theta/2)*sin(theta/2)*sin(phi)^3 + 2*b*r*cos(theta/2)*cos(phi)*sin(phi)^2 + b*r*sin(theta/2)*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)*cos(phi)^2*sin(phi))/(r*(cos(theta/2) + 1)*(b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2 - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi)))
%                                                                                                                                                                                                                                                                                                                                                                                                                                                         (b^2*sin(phi)^3 + 4*r^2*sin(phi)^3 - 3*r^2*cos(theta/2)^2*sin(phi)^3 + r^2*cos(theta/2)^3*sin(phi)^3 - 2*r^2*sin(theta/2)^2*sin(phi)^3 + b^2*cos(phi)^2*sin(phi) + 4*r^2*cos(phi)^2*sin(phi) + r^2*sin(theta/2)*cos(phi) - 4*b*r*sin(phi)^3 + b^2*cos(theta/2)*sin(phi)^3 + b^2*sin(theta/2)*cos(phi)^3 - 5*r^2*sin(theta/2)*cos(phi)^2 + 6*r^2*sin(theta/2)*cos(phi)^3 - 2*r^2*sin(theta/2)*sin(phi)^2 - 2*r^2*cos(phi)*sin(phi) + 2*b*r*cos(theta/2)^2*sin(phi)^3 + b*r*sin(theta/2)^2*sin(phi)^3 - r^2*cos(theta/2)*cos(phi)*sin(phi) - 4*b*r*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)^2*sin(phi)^3 + b^2*cos(theta/2)*cos(phi)^2*sin(phi) + b^2*sin(theta/2)*cos(phi)*sin(phi)^2 + r^2*cos(theta/2)^2*cos(phi)*sin(phi) + 6*r^2*sin(theta/2)*cos(phi)*sin(phi)^2 + r^2*sin(theta/2)^2*cos(phi)*sin(phi) + 2*b*r*sin(theta/2)*cos(phi)^2 - 2*b*r*cos(theta/2)*sin(phi)^3 - 5*b*r*sin(theta/2)*cos(phi)^3 + b*r*sin(theta/2)*sin(phi)^2 - 3*r^2*cos(theta/2)^2*cos(phi)^2*sin(phi) + r^2*cos(theta/2)^3*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)*cos(phi)^2 - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)^3 + b*r*cos(phi)*sin(phi) - 2*r^2*sin(theta/2)^2*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)*sin(phi)^2 + b*r*cos(theta/2)*cos(phi)*sin(phi) + r^2*cos(theta/2)*sin(theta/2)^2*cos(phi)^2*sin(phi) - 2*b*r*cos(theta/2)*cos(phi)^2*sin(phi) - 5*b*r*sin(theta/2)*cos(phi)*sin(phi)^2 - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2 + 2*b*r*cos(theta/2)^2*cos(phi)^2*sin(phi) + b*r*cos(theta/2)*sin(theta/2)*cos(phi)^3 + b*r*sin(theta/2)^2*cos(phi)^2*sin(phi) + b*r*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2)/(2*sin(theta/2)*(r + r*cos(theta/2))*(b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2 - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi)))
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    -(r*sin(theta/2) + b*sin(phi) - 2*r*sin(phi) + r*cos(theta/2)^2*sin(phi) + r*sin(theta/2)^2*sin(phi) + b*cos(theta/2)*sin(phi) + b*sin(theta/2)*cos(phi) - r*cos(theta/2)*sin(phi) - 3*r*sin(theta/2)*cos(phi))/(2*sin(theta/2)*(b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2 - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi)))
%                                                                                                                                                                                                                                                                             -(b^2*sin(phi)^3 + 4*r^2*sin(phi)^3 - 3*r^2*cos(theta/2)^2*sin(phi)^3 + r^2*cos(theta/2)^3*sin(phi)^3 - 2*r^2*sin(theta/2)^2*sin(phi)^3 + b^2*cos(phi)^2*sin(phi) + 4*r^2*cos(phi)^2*sin(phi) + r^2*sin(theta/2)*cos(phi) - 4*b*r*sin(phi)^3 + b^2*cos(theta/2)*sin(phi)^3 + b^2*sin(theta/2)*cos(phi)^3 - 5*r^2*sin(theta/2)*cos(phi)^2 + 6*r^2*sin(theta/2)*cos(phi)^3 - 2*r^2*sin(theta/2)*sin(phi)^2 - 2*r^2*cos(phi)*sin(phi) + 2*b*r*cos(theta/2)^2*sin(phi)^3 + b*r*sin(theta/2)^2*sin(phi)^3 - r^2*cos(theta/2)*cos(phi)*sin(phi) - 4*b*r*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)^2*sin(phi)^3 + b^2*cos(theta/2)*cos(phi)^2*sin(phi) + b^2*sin(theta/2)*cos(phi)*sin(phi)^2 + r^2*cos(theta/2)^2*cos(phi)*sin(phi) + 6*r^2*sin(theta/2)*cos(phi)*sin(phi)^2 + r^2*sin(theta/2)^2*cos(phi)*sin(phi) + 2*b*r*sin(theta/2)*cos(phi)^2 - 2*b*r*cos(theta/2)*sin(phi)^3 - 5*b*r*sin(theta/2)*cos(phi)^3 + b*r*sin(theta/2)*sin(phi)^2 - 3*r^2*cos(theta/2)^2*cos(phi)^2*sin(phi) + r^2*cos(theta/2)^3*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)*cos(phi)^2 - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)^3 + b*r*cos(phi)*sin(phi) - 2*r^2*sin(theta/2)^2*cos(phi)^2*sin(phi) + r^2*cos(theta/2)*sin(theta/2)*sin(phi)^2 + b*r*cos(theta/2)*cos(phi)*sin(phi) + r^2*cos(theta/2)*sin(theta/2)^2*cos(phi)^2*sin(phi) - 2*b*r*cos(theta/2)*cos(phi)^2*sin(phi) - 5*b*r*sin(theta/2)*cos(phi)*sin(phi)^2 - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2 + 2*b*r*cos(theta/2)^2*cos(phi)^2*sin(phi) + b*r*cos(theta/2)*sin(theta/2)*cos(phi)^3 + b*r*sin(theta/2)^2*cos(phi)^2*sin(phi) + b*r*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2)/((b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2 - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi))*(b*sin(phi) - r*sin(theta/2) - 2*r*sin(phi) + r*cos(theta/2)^2*sin(phi) + r*sin(theta/2)^2*sin(phi) + b*cos(theta/2)*sin(phi) - b*sin(theta/2)*cos(phi) - r*cos(theta/2)*sin(phi) + 3*r*sin(theta/2)*cos(phi)))
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                (b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + r*cos(phi) + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2)/(b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2 - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi))
% (b^2*sin(phi)^3 + 4*r^2*sin(phi)^3 - 3*r^2*cos(theta/2)^2*sin(phi)^3 + r^2*cos(theta/2)^3*sin(phi)^3 - 2*r^2*sin(theta/2)^2*sin(phi)^3 - r^2*sin(theta/2)^3*sin(phi)^2 + b^2*cos(phi)^2*sin(phi) + 4*r^2*cos(phi)^2*sin(phi) - 4*b*r*sin(phi)^3 + b^2*cos(theta/2)*sin(phi)^3 + b^2*sin(theta/2)*cos(phi)^3 - 2*r^2*sin(theta/2)*cos(phi)^2 + 6*r^2*sin(theta/2)*cos(phi)^3 - r^2*sin(theta/2)^2*sin(phi) + 2*b*r*cos(theta/2)^2*sin(phi)^3 + b*r*sin(theta/2)^2*sin(phi)^3 + 2*r^2*cos(theta/2)*cos(phi)*sin(phi) - 4*b*r*cos(phi)^2*sin(phi) - r^2*cos(theta/2)^2*sin(theta/2)*sin(phi)^2 + r^2*cos(theta/2)*sin(theta/2)^2*sin(phi)^3 + b^2*cos(theta/2)*cos(phi)^2*sin(phi) + b^2*sin(theta/2)*cos(phi)*sin(phi)^2 + r^2*cos(theta/2)^2*cos(phi)*sin(phi) - r^2*cos(theta/2)^3*cos(phi)*sin(phi) - r^2*cos(theta/2)*sin(theta/2)*cos(phi) + 6*r^2*sin(theta/2)*cos(phi)*sin(phi)^2 + 3*r^2*sin(theta/2)^2*cos(phi)*sin(phi) + b*r*sin(theta/2)*cos(phi)^2 - 2*b*r*cos(theta/2)*sin(phi)^3 - 5*b*r*sin(theta/2)*cos(phi)^3 - 3*r^2*cos(theta/2)^2*cos(phi)^2*sin(phi) + r^2*cos(theta/2)^3*cos(phi)^2*sin(phi) + 4*r^2*cos(theta/2)*sin(theta/2)*cos(phi)^2 - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)^3 - 2*r^2*sin(theta/2)^2*cos(phi)^2*sin(phi) + 2*r^2*cos(theta/2)*sin(theta/2)*sin(phi)^2 - b*r*cos(theta/2)*cos(phi)*sin(phi) + r^2*cos(theta/2)*sin(theta/2)^2*cos(phi)^2*sin(phi) - 2*b*r*cos(theta/2)*cos(phi)^2*sin(phi) - b*r*cos(theta/2)^2*cos(phi)*sin(phi) - 5*b*r*sin(theta/2)*cos(phi)*sin(phi)^2 - b*r*sin(theta/2)^2*cos(phi)*sin(phi) - 3*r^2*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2 - r^2*cos(theta/2)*sin(theta/2)^2*cos(phi)*sin(phi) + 2*b*r*cos(theta/2)^2*cos(phi)^2*sin(phi) - b*r*cos(theta/2)*sin(theta/2)*cos(phi)^2 + b*r*cos(theta/2)*sin(theta/2)*cos(phi)^3 + b*r*sin(theta/2)^2*cos(phi)^2*sin(phi) - b*r*cos(theta/2)*sin(theta/2)*sin(phi)^2 + b*r*cos(theta/2)*sin(theta/2)*cos(phi)*sin(phi)^2)/((b*sin(phi)^2 - 2*r*cos(phi)^2 - 2*r*sin(phi)^2 + b*cos(phi)^2 + r*cos(theta/2)*cos(phi)^2 + r*cos(theta/2)*sin(phi)^2 - r*cos(theta/2)*cos(phi) + r*sin(theta/2)*sin(phi))*(b*sin(phi) - r*sin(theta/2) - 2*r*sin(phi) + r*cos(theta/2)^2*sin(phi) + r*sin(theta/2)^2*sin(phi) + b*cos(theta/2)*sin(phi) - b*sin(theta/2)*cos(phi) - r*cos(theta/2)*sin(phi) + 3*r*sin(theta/2)*cos(phi)))
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               0
% ];
% gam1 = lastA(7);
% gam2 = lastA(8);
% gam3 = lastA(9);


%kappa = 0.05*10^9*pi*(0.05)^2/1;
kappa = 10; %arbitrary stiffness

gamma1 = gamma_1_soln*kappa;
gamma2 = gamma_2_soln*kappa;
gamma3 = gamma_3_soln*kappa;
gamma4 = gamma_4_soln*kappa;
Gamma = diag([gamma1 gamma2 gamma3 gamma4]);


%% initial FQ:dimension works
% FQ_0 = (W -(Q*psi_t+Y)*transpose(C_S)*Gamma*C_S)*psi;

%% construct M matrix in (Q double dot + Q*Xi)*M = F_Q
%try Q with 2 by 8
m = 0.1;%mass
J_tot = (1/4)*m*el^2;%moment of inertia, same for two spines??
J = J_tot/3;
M = diag([m,m,J,J,J,J,J,J]);

% %% solve Q by ode45
% tspan = linspace(0,6);
% Q0_st = [reshape(r1,[],1);reshape(r2,[],1);...
%             reshape(dB,[],1);zeros(24,1)];

%% solve for b1
%Q = [r1,r2,b1,b4] -> 2x4
%Q_stack is 8x1
%Q_st = stack[r1,r2,b1,b4,r1dot,r2dot,b1dot,b4dot] -> 16x1
tspan = 0:0.001:2;
%tspan = linspace(0,6);
R = [r1,r2];
solB = [b1,b4];
Q0_st = [reshape(R,[],1);reshape(solB,[],1);zeros(8,1)];

% Create function to be solved by the ODE by being run repeatedly
% Takes in t vector, current Q and Q at the next time step
% Q is a 16x1 column vector, with the first 8 being the r and b 1 and 2,
% the next 8 are zeroes, repp
odefun = @(t,Q_st)[Q_st(9:16);...
                    deri(W,psi,C_S,Y,Gamma,J, Q_st,m)];
options = odeset('RelTol',1e-13,'AbsTol',1e-13,'Refine',10);
[t,Q_sol] = ode45(odefun,tspan,Q0_st,options);
%[t,Q_sol] = ode45(odefun,tspan,Q0_st);

figure;
plot(t,Q_sol(:,1),t,Q_sol(:,2),t,Q_sol(:,3),t,Q_sol(:,4),'linew',1.5)
legend('$r_{1x}$','$r_{1y}$','$r_{2x}$','$r_{2y}$','interpreter','latex','Location','best');
title("Position vectors for both spines")
grid on;
xlabel('Time (s)')
ylabel('Position (m)')
set(gcf,'Color','w')
set(gca,'FontSize',14)
set(gca,'FontName','ComputerModern')
set(gca, 'box', 'off')

figure;
plot(t,rotn120*Q_sol(:,[5 6])',t,rotn120*Q_sol(:,[7 8])','linew',1.5)
legend('$b_{3x}$','$b_{3y}$','$b_{6x}$','$b_{6y}$','interpreter','latex','location','best');
title("Orientation vectors for both spines")
grid on;
xlabel('Time (s)')
ylabel('unitless')
set(gcf,'Color','w')
set(gca,'FontSize',14)
set(gca,'FontName','ComputerModern')
set(gca, 'box', 'off')
axis([tspan(1) tspan(end) -0.2 1.2])

%% function
function y = deri(W,psi,C_S,Y,Gamma,J, Q_st,m)
rotp120 = [cos(2*pi/3) -sin(2*pi/3); sin(2*pi/3) cos(2*pi/3)];
rotn120 = [cos(-2*pi/3) -sin(-2*pi/3); sin(-2*pi/3) cos(-2*pi/3)];
Q_all = [reshape(Q_st(1:6),2,[]),rotp120*Q_st(5:6),rotn120*Q_st(5:6),...
    Q_st(7:8),rotp120*Q_st(7:8),rotn120*Q_st(7:8)];


F_Q = (W - (Q_all*transpose(psi)+Y)*transpose(C_S)*Gamma*C_S)*psi;
fr1 = F_Q(:,1);
fr2 = F_Q(:,2);
fb1 = F_Q(:,3);
fb2 = F_Q(:,4);
fb3 = F_Q(:,5);
fb4 = F_Q(:,6);
fb5 = F_Q(:,7);
fb6 = F_Q(:,8);

rot3p120 = [cos(2*pi/3) -sin(2*pi/3) 0; sin(2*pi/3) cos(2*pi/3) 0; 0 0 1];
rot3n120 = [cos(-2*pi/3) -sin(-2*pi/3) 0; sin(-2*pi/3) cos(-2*pi/3) 0; 0 0 1];

b1 = [Q_st(5);Q_st(6)];
b1_m = [Q_st(5);Q_st(6);0];
b1dot = [Q_st(13),Q_st(14)];
tau_1m = cross(b1_m,[fb1;0])+cross(rot3p120*b1_m,[fb2;0])+cross(rot3n120*b1_m,[fb3;0]);
tau_1 = tau_1m(3);

b4 = [Q_st(7);Q_st(8)];
b4_m = [Q_st(7);Q_st(8);0];
b4dot = [Q_st(15),Q_st(16)];
tau_2m = cross(b4_m,[fb4;0])+cross(rot3p120*b4_m,[fb5;0])+cross(rot3n120*b4_m,[fb6;0]);
tau_2 = tau_2m(3);

cross1 = cross(b1_m,tau_1m);
cross4 = cross(b4_m,tau_2m);

b1doubledot = -((cross1(1:2))/(J*(norm(b1))^2))-b1*(norm(b1dot)^2)/(norm(b1)^2);
b4doubledot = -((cross4(1:2))/(J*(norm(b4))^2))-b4*(norm(b4dot)^2)/(norm(b4)^2);

r1doubledot = fr1/m;
r2doubledot = fr2/m;

y = [r1doubledot;r2doubledot;b1doubledot;b4doubledot];

end


