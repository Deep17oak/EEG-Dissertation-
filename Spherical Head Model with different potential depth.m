
%% Spherical Head Model FEM Sweep: Multiple Dipole Orientations & Depths
clc; clear; close all;

% Load the spherical head mesh (assumes it defines p, t, r, D, etc.)
load Spherical_mesh;

%% Define source positions and dipole moments
positions = linspace(0.001, 0.06, 5);      % Generate 5 depths along Z between 1 mm and 60 mm


% Define three unit dipoles along X, Y, Z (scaled to 1e-8 A·m)
dipolar_moments = [       
    1, 0, 0;     % X‐oriented             
    0, 1, 0;     % Y‐oriented
    0, 0, 1      % Z‐oriented
    ] * 1e-8;

% Store counts for looping and subplot indexing
num_pos    = numel(positions);
num_dipole = size(dipolar_moments, 1);

%%  Adjust Scalp Conductivity (One‐Time)
ind_layers = (t(:,5) == 1);
D(ind_layers, [1,4,6]) = 0.03;  % Scalp element diagonal conductivities to 0.03 S/m

% Assemble global stiffness matrix M from nodes p, elements t, and tensors D
M_fo       = femeg_stiffness(p, t, D);
[L_fo, U_fo] = ilu(M_fo);

%% Precompute scalp‐node mask 
ind_scalp = sqrt(sum(p.^2,2)) > (r(1) - 1e-5);
if nnz(ind_scalp) == 0
    error('No scalp nodes found – check mesh radius r and tolerance.');
end

%% Solve & Plot for Each Dipole Orientation and Depth
figure;
for j = 1:num_dipole
    q = dipolar_moments(j, :);   % Current dipole orientation
    for i = 1:num_pos
        pos = [0, 0, positions(i)];  % Dipole location along Z
        
        % Build RHS and analytic infinite‐medium solution
        b_fo = femeg_indep_analyt(p, t, pos, q, D);
        uinf = femeg_uinf(p, pos, q, 0.41);  
        
        % Solve with QMR + preconditioner
        [u_n, flag_n] = qmr(M_fo, b_fo, 1e-10, 4000, L_fo, U_fo);
        if flag_n ~= 0
            warning('QMR failed to converge (flag = %d) for dipole %d at Z=%.4f m',flag_n, j, positions(i));
        end
        
        % Total potential and reference‐subtract on scalp
        u = u_n + uinf;
        u_ref_sub = u - mean(u(ind_scalp));
        
        % Subplot index: dipole along rows, position along columns
        subplot(num_dipole, num_pos, (j-1)*num_pos + i);
        femeg_vis3d(p, t, u_ref_sub);
        caxis([-0.7, 0.7] * 1e-6);
        colorbar;
        title(sprintf('q=[%g %g %g] at Z=%.3f m', q, positions(i)),'FontSize', 10);
    end
end

% Global figure formatting
colormap(parula);
axis equal off;
sgtitle('Scalp potentials for various dipole orientations & depths','FontSize', 14);

%femeg_vis3d: Renders a volumetric/isosurface plot of the scalar field u_fsc on the tetrahedral mesh.

%  Result: A 3×5 grid of 3D scalp potential maps showing how dipoles along X, Y, and Z at five different depths (1 mm–60 mm)
% produce distinct voltage distributions.