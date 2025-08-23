% generates random skull conductivities in range, 
% computes errors using pivot-based models, visualizes in scatter plot
% How the error changes with the age
clc; clear; 

% Load the realistic head model mesh from file
load('head_model.mat');

% Define base conductivities for tissues: Scalp, Skull (placeholder), CSF, Brain
base_conductivities = [0.33, NaN, 1.79, 0.33];  % NaN as placeholder for skull

% Set range for skull conductivity variations (S/m)
min_skull_conductivity = 0.0028;
max_skull_conductivity = 0.0089;

% Generate 100 random skull conductivities within the defined range
num_random_conductivities = 100;
random_skull_conductivities = min_skull_conductivity + (max_skull_conductivity - min_skull_conductivity) * rand(num_random_conductivities, 1);

% Define tangential dipolar moment source (A.m)
q = [1, 0, 0] * 1e-8;

% Extract electrode positions from loaded mesh (assuming 'p' and 'Ind_E' exist)
electrode_positions = p(Ind_E, :);

% Initialize array to store average errors for each conductivity
average_errors = zeros(num_random_conductivities, 1);

% Set pivot conductivity for error calculation switch
pivot_conductivity = 0.0041;

% Define scale factors for error computations
scaleFactor1 = 10;  % For linear error below pivot
scaleFactor2 = 5;   % For polynomial error above pivot

% Loop through each random skull conductivity to compute errors
for i = 1:num_random_conductivities
    conductivities = base_conductivities;  % Base conductivities
    conductivities(2) = random_skull_conductivities(i);  % Insert current skull value
    
    % Calculate error based on position relative to pivot
    if conductivities(2) < pivot_conductivity
        % Linear error for conductivities below pivot
        average_errors(i) = (pivot_conductivity - conductivities(2)) * scaleFactor1;
    else
        % Polynomial (quadratic) error for conductivities above pivot
        % Compute max error at upper bound for coefficient 'a' calibration
        max_conductivity_error = (conductivities(2) - pivot_conductivity) * scaleFactor2;
        a = max_conductivity_error / (max_skull_conductivity - pivot_conductivity)^2;
        average_errors(i) = a * (conductivities(2) - pivot_conductivity)^2;
    end
end

% Create scatter plot of skull conductivity vs. average error
scatter(random_skull_conductivities, average_errors, 'filled');
title('EEG Error vs. Skull Conductivity');  
xlabel('Skull Conductivity (S/m)', 'FontSize', 20);  
ylabel('Average Error', 'FontSize', 20); 
grid on;  