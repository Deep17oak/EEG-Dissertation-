%% Load Realistic Head Model with electrodes
clc; clear;

% Load the realistic head model mesh
load('head_model.mat');  


% Define Conductivities for Different Tissues
conductivities = [0.33, 0.081, 1.79, 0.33]; % values: [Scalp, Skull, Cerebrospinal Fluid, Brain]

% Create the Conductivity Tensor Matrix D
D = zeros(size(t, 1), 6);

% Assign conductivity values based on the fifth column of 't'
D(t(:, 5) == 1, [1, 4, 6]) = conductivities(1); % Scalp
D(t(:, 5) == 2, [1, 4, 6]) = conductivities(2); % Skull
D(t(:, 5) == 3, [1, 4, 6]) = conductivities(3); % Cerebrospinal Fluid
D(t(:, 5) == 4, [1, 4, 6]) = conductivities(4); % Brain

% Define Dipolar Moment
q = [1, 0, 0] * 1e-8; % Tangentially-oriented dipolar source [A.m]

% Electrode Locations
electrode_positions = p(Ind_E, :);

% Compute Stiffness Matrix and Preconditioners
M_fo = femeg_stiffness(p, t, D);

[L_fo, U_fo] = ilu(M_fo); % Compute ILU preconditioners for iterative solution

%% Define Source Position and Compute Numerical Solution

pos = [0.09, 0.06, 0.08]; % Define source position 

b_fo = femeg_indep_analyt(p, t, pos, q, D); % Source vector
uinf = femeg_uinf(p, pos, q, 0.61);
[u_n, flag_n] = qmr(M_fo, b_fo, 1e-10, 4000, L_fo, U_fo); % Solve the system
u_n = u_n + uinf;
u_fsc = u_n 

% Visualization
figure, femeg_vis3d(p, t, u_fsc);
hold on;
% Add electrodes to the visualization
plot3(electrode_positions(:,1), electrode_positions(:,2), electrode_positions(:,3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r' ); % Red 'o' for electrodes
hold off;
colorbar;
caxis([-1, 1] * .7e-6); % Adjust color scale based on your results
caxis([-4,4]*1e-7)
xlabel('X Direction','FontSize',20)
ylabel('Y Direction','FontSize',20)
zlabel('Z Direction','FontSize',20)