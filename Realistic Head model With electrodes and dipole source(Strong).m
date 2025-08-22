%% Load Realistic Head Model with electrodes and Dipole source(Strong)
clc; clear;  % Clear command window and workspace


% This section loads a precomputed head model mesh, assigns tissue conductivities,
load('head_model.mat');  % Load mesh data: points 'p', tetrahedra 't', electrode indices 'Ind_E'

% Define Conductivities for Different Tissues
conductivities = [0.33, 0.081, 1.79, 0.33];  % In S/m: [Scalp, Skull, CSF, Brain]

% Create the Conductivity Tensor Matrix D
% D is a matrix for anisotropic/isotropic conductivities per tetrahedron.
D = zeros(size(t, 1), 6);  

% Assign isotropic conductivity values based on tissue labels in t(:,5)
D(t(:,5) == 1, [1,4,6]) = conductivities(1);  % Scalp tissue
D(t(:,5) == 2, [1,4,6]) = conductivities(2);  % Skull tissue
D(t(:,5) == 3, [1,4,6]) = conductivities(3);  % Cerebrospinal Fluid (CSF)
D(t(:,5) == 4, [1,4,6]) = conductivities(4);  % Brain tissue

% Define Dipolar Moment
q = [1, 0, 0] * 1e-8;  % Tangential dipole moment in A.m 

% Electrode Locations
electrode_positions = p(Ind_E, :);  % Extract positions from mesh points

% Compute Stiffness Matrix and Preconditioners
M_fo = femeg_stiffness(p, t, D);  % Finite element stiffness matrix for EEG forward problem
[L_fo, U_fo] = ilu(M_fo);  % Incomplete LU factorization for preconditioning

%% Define Source Position and Compute Numerical Solution
pos = [0.09, 0.06, 0.08];  % Dipole source position in meters (x,y,z)
b_fo = femeg_indep_analyt(p, t, pos, q, D);  % Compute source term vector
uinf = femeg_uinf(p, pos, q, 0.61);  % Infinite medium potential (analytical)
[u_n, flag_n] = qmr(M_fo, b_fo, 1e-10, 4000, L_fo, U_fo);  % Solve linear system using QMR iterative solver
u_n = u_n + uinf;  % Add infinite medium correction
u_fsc = u_n;  % Final potential solution

% Visualization
figure;  % Create new figure
femeg_vis3d(p, t, u_fsc);  % Visualize 3D potential on head model
hold on;
% Add electrodes to the visualization
plot3(electrode_positions(:,1), electrode_positions(:,2), electrode_positions(:,3), ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Plot electrodes as red circles
hold off;
colorbar;  % Add colorbar for potential values
caxis([-1, 1] * .7e-6);  % Set color axis limits
caxis([-4,4]*1e-7);  % Override with tighter limits for better visualization
xlabel('X Direction','FontSize',20);  % Label axes
ylabel('Y Direction','FontSize',20);
zlabel('Z Direction','FontSize',20);