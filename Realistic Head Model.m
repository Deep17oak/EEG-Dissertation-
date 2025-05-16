% Realistic Head model with potential 
clc; clear;

% Load the realistic head model mesh
load('head_model.mat');  % This file  contain 'p' (points) and 't' (tetrahedral elements)

% Define Conductivities for Different Tissues 
conductivities = [0.33, 0.0041, 1.79, 0.33]; % [Scalp, Skull, Cerebrospinal Fluid, Brain]

% Create the Conductivity Tensor Matrix D
D = zeros(size(t, 1), 6);

% Assign conductivity values based on the tissue type
D(t(:, 5) == 1, [1, 4, 6]) = conductivities(1); % Scalp
D(t(:, 5) == 2, [1, 4, 6]) = conductivities(2); % Skull
D(t(:, 5) == 3, [1, 4, 6]) = conductivities(3); % Cerebrospinal Fluid
D(t(:, 5) == 4, [1, 4, 6]) = conductivities(4); % Brain

% Define Dipolar Moment
q = [1, 0, 0] * 1e-8; % Tangentially-oriented dipolar source [A.m]

% Compute Stiffness Matrix
% Use the femeg_stiffness function 
M_fo = femeg_stiffness(p, t, D);

% Compute Incomplete LU Factorization for Preconditioning
[L_fo, U_fo] = ilu(M_fo);

% Define Source Position
pos = [0.09, 0.06, 0.08]; 

% Compute Source Vector
% Use the femeg_indep_analyt function 
b_fo = femeg_indep_analyt(p, t, pos, q, D);

% Compute Infinite Medium Potential
% Use the femeg_uinf function 
uinf = femeg_uinf(p, pos, q, 0.61);

% Solve the Linear System
[u_n, flag_n] = qmr(M_fo, b_fo, 1e-10, 4000, L_fo, U_fo); % Using QMR solver
u_n = u_n + uinf; % Add the infinite medium potential

% Visualization
% Use the femeg_vis3d function 
figure, femeg_vis3d(p, t, u_n);
colorbar;
caxis([-1, 1] * .7e-6);