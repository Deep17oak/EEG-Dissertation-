% Spherical Head Model FEM with Dipole Source and 3D Visualization
%% Initialization
clc; clear;

% Load the mesh data for the spherical model
%  p – node coordinates (Nx3)
%  t – tetrahedral elements (Mx5: nodes)
%  D – conductivity/diffusion tensors (Mx6)
%  r – radii of scalp, skull outer, skull inner, brain shells
%  cond – conductivities for each shell
load Spherical_mesh;

%% Define one dipolar source
pos = [0, 0, 0.04];         % Position (m) of dipole inside sphere along Z-axis (4 cm depth)
q   = [0, 0, 1] * 1e-8;    % Dipole moment (C·m) oriented along Z

    
%% Adjust scalp diffusivity and identify scalp nodes
% Set diffusivity for scalp elements (where t(:,5)==1)
ind_ele_scalp = t(:,5) == 1;
D(ind_ele_scalp, [1,4,6]) = 0.03;

% Mark scalp nodes for baseline subtraction
ind_node_scalp = sqrt(sum(p.^2, 2)) > r(1) - 1e-5;

%Build global stiffness matrix M based on nodes p, elements t, and tensors D
M_fo    = femeg_stiffness(p, t, D);

% Compute incomplete LU factorization for preconditioning
[L_fo, U_fo] = ilu(M_fo);

% Assemble the source vector b from the analytic dipole contribution
b_fo = femeg_indep_analyt(p, t, pos, q, D);

%Compute the closed‐form infinite‐medium potential at each node using a homogeneous conductivity of 0.41 S/m
uinf = femeg_uinf(p, pos, q, 0.41);

% Solve for potential using QMR and add infinite‐medium term
tol    = 1e-10;  % Solver tolerance
maxit  = 4000;   % Maximum number of iterations
[u_num, flag] = qmr(M_fo, b_fo, tol, maxit, L_fo, U_fo);

% Combine numerical solution with analytic infinite‐medium part
u_total      = u_num + uinf;

% Reference to zero mean on scalp
u_fsc = u_total - mean(u_total(ind_node_scalp));

%% 3D Visualization of Scalp Potentials
figure;
femeg_vis3d(p, t, u_fsc);
axis equal off;
colormap(parula);
colorbar('FontSize',12);
clim([-0.7, 0.7] * 1e-6);

% Title with dipole depth and moment details
title(sprintf('Dipole at Z=%.2f m, q=[%.0e %.0e %.0e]', pos(3), q(1), q(2), q(3)), 'FontSize',14);
