% Initialization Spherical Head model with the depth on dipole position
%% Initialization
clc; clear;

% Load the mesh data for the spherical model
load Spherical_mesh;

%% Define one dipolar source
pos = [0, 0, 0.04];         % Position (m) along Z‐axis inside the sphere
q   = [0, 0, 1] * 1e-8; % Dipole moment (C·m) oriented along Z

    

%% Adjust scalp diffusivity and identify scalp nodes
% Set diffusivity for scalp elements (where t(:,5)==1)
ind_ele_scalp = t(:,5) == 1;
D(ind_ele_scalp, [1,4,6]) = 0.03;

% Mark scalp nodes for baseline subtraction
ind_node_scalp = sqrt(sum(p.^2, 2)) > r(1) - 1e-5;

%% Assemble FEM matrices and preconditioner
M_fo    = femeg_stiffness(p, t, D);
[L_fo, U_fo] = ilu(M_fo);

%% Build right‐hand side and analytical infinite‐medium term
b_fo = femeg_indep_analyt(p, t, pos, q, D);
uinf = femeg_uinf(p, pos, q, 0.41);

%% Solve for potential using QMR and add infinite‐medium term
tol    = 1e-10;
maxit  = 4000;
[u_num, flag] = qmr(M_fo, b_fo, tol, maxit, L_fo, U_fo);
u_total      = u_num + uinf;

% Reference to zero mean on scalp
u_fsc = u_total - mean(u_total(ind_node_scalp));

%% 3D Visualization
figure;
femeg_vis3d(p, t, u_fsc);
axis equal off;
colormap(parula);
colorbar('FontSize',12);
clim([-0.7, 0.7] * 1e-6);
title(sprintf('Dipole at Z=%.2f m, q=[%.0e %.0e %.0e]', pos(3), q(1), q(2), q(3)), 'FontSize',14);
