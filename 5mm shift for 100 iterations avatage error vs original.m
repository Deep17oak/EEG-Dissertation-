% Standard deviation of 5 mm for 100 iterations: shift electrode positions, project to scalp surface, average error vs original.

El = p(Ind_E, :);  % Original electrode positions

% Extract scalp surface points (tissue label 1)
scalp_faces = t(t(:,5)==1,1:4);  % Scalp tetrahedra faces
scalp_points_idx = unique(scalp_faces(:));  % Unique point indices on scalp
surface_p = p(scalp_points_idx, :);  % Scalp surface points

% Initialize matrix for perturbed positions
El_1 = zeros(size(El,1), size(El,2), 100);

% Generate 100 sets of random variations
for i = 1:100
    perturbedEl = El + 0.005 * randn(size(El));  % Add 5 mm std dev noise
    nearestMeshPoints = dsearchn(surface_p, perturbedEl);  % Find nearest on surface
    projectedEl = surface_p(nearestMeshPoints, :);
    El_1(:,:,i) = projectedEl;
end

% Compute average error
differences = El_1 - repmat(El, [1 1 100]);  % Differences per set
errors = sqrt(sum(differences.^2, 2));  % Euclidean distances per electrode per set
avg_errors = mean(errors, 3);  % Average over 100 sets per electrode
overall_avg_error = mean(avg_errors(:));  % Overall average error
disp(['Overall average error: ', num2str(overall_avg_error), ' m']);

% Plotting (first 10 sets for clarity)
figure;
hold on;
scatter3(El(:,1), El(:,2), El(:,3), 'filled', 'r');  % Original in red
for i = 1:100
    projectedPositions = El_1(:,:,i);
    scatter3(projectedPositions(:,1), projectedPositions(:,2), projectedPositions(:,3), 'filled', 'b');  % Perturbed in blue
end
xlabel('X Axis','FontSize',20);
ylabel('Y Axis','FontSize',20);
zlabel('Z Axis','FontSize',20);
grid on;
hold off;