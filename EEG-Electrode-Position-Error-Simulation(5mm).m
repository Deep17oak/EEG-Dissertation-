% This code simulates 5 mm Gaussian noise perturbations on EEG electrode positions, 
% projects them to the scalp surface, computes boundary correction errors (u_fsc - uinf) for 100 iterations, 
% and plots the original error in thick black atop overlaid red perturbed
% errors with a legend.
El = p(Ind_E, :);
% Number of electrodes
num_electrodes = size(El, 1);
% Number of iterations
num_iterations = 100;

% Extract scalp surface points (tissue label 1)
scalp_faces = t(t(:,5)==1,1:4);  % Scalp tetrahedra faces
scalp_points_idx = unique(scalp_faces(:));  % Unique point indices on scalp
surface_p = p(scalp_points_idx, :);  % Scalp surface points

% Original 'error' (boundary correction) at electrodes
error_original = u_fsc(Ind_E) - uinf(Ind_E);

% Initialize matrix to store perturbed errors
errors_perturbed = zeros(num_electrodes, num_iterations);

% Calculate perturbations for each iteration
for i = 1:num_iterations
    for j = 1:num_electrodes
        % Generate random perturbation
        perturbation = 0.005 * randn(1, 3);
        % Apply the perturbation to the electrode's position
        perturbed_pos = El(j, :) + perturbation;
        % Find the nearest neighbor on the surface
        nearest_surface_idx = dsearchn(surface_p, perturbed_pos);
        nearest_index = scalp_points_idx(nearest_surface_idx);
        % Calculate the 'error' at perturbed (projected) position
        errors_perturbed(j, i) = u_fsc(nearest_index) - uinf(nearest_index);
    end
end

% Compute average error magnitude over all iterations and electrodes
diffs = errors_perturbed - repmat(error_original, 1, num_iterations);
avg_error_magnitude = mean(abs(diffs(:)));

% Plotting: all perturbed in red (overlaid, thinner), original in black last (thicker) for visibility
figure;
hold on;
h1 = plot(1:num_electrodes, errors_perturbed(:,1), 'r', 'LineWidth', 0.5);  % First red for legend
for i = 2:num_iterations
    plot(1:num_electrodes, errors_perturbed(:, i), 'r', 'LineWidth', 0.5);  % Remaining reds
end
h2 = plot(1:num_electrodes, error_original, 'k', 'LineWidth', 1.5);  % Original in black, thicker, on top
hold off;
xlabel('Electrodes Index');
ylabel('Measured Error (Volts) \times 10^{-7}');
title('Error Due to Electrode Position Change 5mm (100 Iterations)');
grid on;
legend([h2 h1], {'Black = original', 'Red = random error'}, 'Location', 'northeast');
