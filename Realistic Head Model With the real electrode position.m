% Realistic Head Model with the real electrode prosition 
clc, clear

% and 'Ind_E' (indices of electrode positions in 'p').
load('head_model.mat');

% Define Conductivities for Different Tissues (values in S/m)
conductivities = [0.33, 0.0041, 1.79, 0.33]; % [Scalp, Skull, Cerebrospinal Fluid, Brain]

% Create the Conductivity Tensor Matrix D
D = zeros(size(t, 1), 6);

% Assign conductivity values based on the tissue type (t(:,5) labels: 1=Scalp, 2=Skull, 3=CSF, 4=Brain)
D(t(:, 5) == 1, [1, 4, 6]) = conductivities(1); % Scalp
D(t(:, 5) == 2, [1, 4, 6]) = conductivities(2); % Skull
D(t(:, 5) == 3, [1, 4, 6]) = conductivities(3); % Cerebrospinal Fluid
D(t(:, 5) == 4, [1, 4, 6]) = conductivities(4); % Brain

% Define Dipolar Moment
% q is the dipole moment vector [qx, qy, qz] in A.m (Ampere-meter), here tangential (x-direction) with magnitude 1e-8.
q = [1, 0, 0] * 1e-8;

% Extract coordinates of electrodes from vertices 'p' using indices 'Ind_E'.
electrode_positions = p(Ind_E, :);

% Compute Stiffness Matrix and Preconditioners
% femeg_stiffness: Custom function to compute the FEM stiffness matrix from points, elements, and conductivities.
M_fo = femeg_stiffness(p, t, D);
% Incomplete LU decomposition for preconditioning the linear solver.
[L_fo, U_fo] = ilu(M_fo);

% Define Source Position and Compute Numerical Solution
% pos: XYZ coordinates of the dipole source (in meters, assuming model units).
pos = [0.09, 0.06, 0.08];
% femeg_indep_analyt: Custom function to compute the independent (source) vector b for the FEM system.
b_fo = femeg_indep_analyt(p, t, pos, q, D);
% femeg_uinf: Custom function to compute the potential in an infinite homogeneous medium (conductivity 0.61 S/m?).
uinf = femeg_uinf(p, pos, q, 0.61);
% Solve the sparse linear system M_fo * u = b_fo using Quasi-Minimal Residual (QMR) method,
% with tolerance 1e-10, max iterations 4000, and ILU preconditioners.
[u_n, flag_n] = qmr(M_fo, b_fo, 1e-10, 4000, L_fo, U_fo);
% Add the infinite medium potential to the solution (subtraction method).
u_n = u_n + uinf;
u_fsc = u_n;

% Extract Potential at Electrode Locations
% Potentials at electrode vertices.
electrode_potentials = u_fsc(Ind_E);

% Visualization
% Create a new figure and visualize the 3D potential distribution using custom function.
figure, femeg_vis3d(p, t, u_fsc);
hold on;

% Add electrodes to the visualization
plot3(electrode_positions(:,1), electrode_positions(:,2), electrode_positions(:,3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;

% Add colorbar and set color axis limits based on expected potential range.
colorbar;
caxis([-1, 1] * .7e-6); % Adjust color scale based on your results

% Display Electrode Potentials
disp('Electrode Potentials:');
disp(electrode_potentials);