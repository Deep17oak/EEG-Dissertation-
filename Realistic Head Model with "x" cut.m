%% %% Load Realistic Head Model Cut in "x" direction for the visualization of the  dipole source
clc; clear;

% Load the realistic head model mesh
load('head_model.mat');   % 'head_model.mat' contains variables 'p' (points) and 't' (tetrahedral elements).

% Define Conductivities for Different Tissues
conductivities = [0.33, 0.0041, 1.79, 0.33]; % [Scalp, Skull, Cerebrospinal Fluid, Brain]


%Conductivity Tensor Matrix D
D = zeros(size(t, 1), 6);


% Assign conductivity values
D(t(:, 5) == 1, [1, 4, 6]) = conductivities(1); % Scalp
D(t(:, 5) == 2, [1, 4, 6]) = conductivities(2); % Skull
D(t(:, 5) == 3, [1, 4, 6]) = conductivities(3); % Cerebrospinal Fluid
D(t(:, 5) == 4, [1, 4, 6]) = conductivities(4); % Brain

% Define Dipole Moment
q = [1, 0, 0] * 1e-8; % Tangentially-oriented dipolar source [A.m]



% Selected Points on the Scalp for Reference Electrode Computation
ind_sc = find(t(:,5) == 1); 


% Compute Stiffness Matrix and Preconditioners
M_fo = femeg_stiffness(p, t, D);   % Adapt the femeg_stiffness function for realistic head model
[L_fo, U_fo] = ilu(M_fo);


% Define Source Position 
pos = [0.09, 0.06, 0.08];   

b_fo = femeg_indep_analyt(p, t, pos, q, D); % Source vector

%The last parameter is the background conductivity- 0.61 S/m
uinf = femeg_uinf(p, pos, q, 0.61);

% Solve the linear system for the electric potential using the Quasi-Minimal Residual (QMR) 
[u_n, flag_n] = qmr(M_fo, b_fo, 1e-10, 4000, L_fo, U_fo); 
u_n = u_n + uinf;

% Visualization
figure, femeg_vis3d(p, t, u_n,'p(:,1)>0.1');
colorbar;
caxis([-1, 1] * .7e-6); % Adjust color scale based on your results
