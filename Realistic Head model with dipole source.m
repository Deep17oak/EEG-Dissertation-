%% Creating Load Realistic Head Model with dipole source
clc; clear;
load('head_model.mat');  % Load the realistic head model mesh

% Define Conductivities for Different Tissues
conductivities = [0.33, 0.0041, 1.79, 0.33]; % Values= [Scalp,skull,CSF,Brain]

%creating the conductivity tensor Matrix D
D = zeros(size(t,1),6);

%  'head_model.mat' contains variables 'p' (points) and 't' (tetrahedral elements),
D(t(:,5)==1,[1 ,4, 6])=conductivities(1); %Scalp
D(t(:,5)==2,[1 ,4, 6])=conductivities(2); % Skull 
D(t(:,5)==3,[1 ,4, 6])=conductivities(3); % CSF
D(t(:,5)==4,[1 ,4, 6])=conductivities(4); % Brain 

%Defining Dipolar Moment
q = [1, 0, 0] * 1e-8;
 
% Select Pointson the Scalp for refenenceelectrode computation
ind_sc = find(t(:,5) == 1);

% The last column of 't' indicates tissue type, with 1 for scalp
ind_sc = find(t(:,5) == 1); 

% Compute Stiffness Matrix and Preconditioners
% Adapt the femeg_stiffness function for realistic head model nodes-p, elements-t, conductivity-D
M_fo = femeg_stiffness(p, t, D);
[L_fo, U_fo] = ilu(M_fo);

% Define Source Position and Compute Numerical Solution
pos = [0.09, 0.07, 0.09];  % Define source (dipole) position in head coordinates (meters)

b_fo = femeg_indep_analyt(p, t, pos, q, D); % femeg_indep_analyt computes contributions of dipole at pos with moment q

uinf = femeg_uinf(p, pos, q, 0.61); % Compute the analytical potential in an infinite homogeneous medium

[u_n, flag_n] = qmr(M_fo, b_fo, 1e-10, 4000, L_fo, U_fo); % Solve the system
u_n = u_n + uinf;
u_fsc = u_n

% Visualization: 3D plot of potential on head mesh
figure, femeg_vis3d(p, t, u_fsc);
colorbar ;
caxis([-1, 1] * .7e-6); % Adjust color scale 
caxis([-4,4]*1e-7)
caxis([-2,2]*1e-7)
xlabel('X Axis','FontSize',20)
ylabel('Y Axis','FontSize',20)
zlabel('Z Axis','FontSize',20)