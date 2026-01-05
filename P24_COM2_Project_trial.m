% --- MATLAB Code: LTE Coverage with Spatially-Correlated Log-Normal Shadowing ---
% --- Main Sim File

% Parameter Realistic 3GPP Value (Suburban Macro)	Notes
% Pt (Tx Power)	20 to 40 (Linear Watts)	43 to 46 dBm is standard for Macro cells.
% Gt (Tx Gain)	63 (Linear, ~18 dBi)	Macro sectors use high-gain directional antennas.
% Gr (Rx Gain)	1 (Linear, 0 dBi)	Mobile phones have omni antennas (0 dBi gain).
% f_MHz		2100 (or 1800/800)	Band 1 (2100) is standard for capacity layers.
% L (Sys Loss)	2 to 4 (Linear, 3-6 dB)	Accounts for cable, connector, and body losses.
% MAPL_dB	140 to 148	132 is okay for high data rates; 148 is for edge voice coverage.
% pathlossExp	3.5 to 4.0	3.5 is standard suburban; 4.0 is urban.
% std_db	8.0	6-10 dB is standard for suburban shadowing.
% correlation_distance	0.05 (50m) to 0.1 (100m)	
% R_theoretical	Calculated	You must recalculate this if you change Pt/Gt/MAPL.

clear; clc; close all;

%% 1. Simulation and System Setup
n_users = 2500;           % More users for better resolution
simulation_area_km = 10.0; % Area size (+/- 3km box)

% System Parameters
f_MHz = 2100;  %3GPP LTE Band1
MAPL_dB = 132;            % Maximum Allowable Path Loss (Uplink Limited)
R_theoretical = 1.31;     % The "perfect" circle radius for reference

% Log-Normal Shadowing Parameters
pathlossExp = 3.5;        % Path loss exponent (typical for suburban: 3-5)
std_db = 8.0;             % Shadowing standard deviation (typical: 6-12 dB)
d0 = 0.1;                 % Reference distance in km (100m)
correlation_distance = 0.2; % Decorrelation distance in km (typical: 50-500m)

% Shadowing Model Selection
% 1 = Independent (original randn - most realistic, most conservative)
% 2 = Distance-correlated (rng based on distance - unrealistic)
% 3 = Spatially-correlated (best - models real-world spatial correlation)
shadowing_model = 3;

% Radio Parameters for Friis calculation
c = 3e8;                  % Speed of light (m/s)
lambda = c / (f_MHz * 1e6); % Wavelength in meters
Pt = 20;                   % Transmit power (normalized, since we work in dB)
Gt = 63;                    % Transmit antenna gain (linear)
Gr = 1;                  % Receive antenna gain (linear)
L = 2;                    % System loss factor (linear)

%% 2. Define Obstacles (Mountains) - Optional for visualization
% Format: [X_center(km), Y_center(km), Radius(km), Attenuation(dB)]
obstacles = [
    0.5,   0.5,   0.3,   3;  % Small Mountain
   -0.8,   0.4,   0.4,   6;  % Average Mountain
    0.2,  -0.7,   0.5,   10;  % Large Mountain
];
n_obstacles = size(obstacles, 1);

%% 3. Generate Users and Calculate Distance
% Generate random user coordinates
rng(42); % Fixed seed for reproducibility
user_x = (rand(n_users, 1) * 2 * simulation_area_km) - simulation_area_km;
user_y = (rand(n_users, 1) * 2 * simulation_area_km) - simulation_area_km;

% Calculate distance from tower (0,0)
d_km = sqrt(user_x.^2 + user_y.^2) + 0.0001; % Add epsilon to avoid log(0)

%% 4. Calculate Path Loss using Selected Shadowing Model

PL_total = zeros(n_users, 1);

if shadowing_model == 1
    %% MODEL 1: Independent Shadowing (Original - Most Realistic)
    fprintf('Using MODEL 1: Independent Shadowing (Original randn method)\n');
    
    for i = 1:n_users
        Pr = log_normal_shadowing_independent(Pt, Gt, Gr, lambda, L, ...
                                              pathlossExp, std_db, d0, d_km(i));
        PL_total(i) = -10 * log10(Pr + eps);
    end
    
elseif shadowing_model == 2
    %% MODEL 2: Distance-Correlated Shadowing (Unrealistic)
    fprintf('Using MODEL 2: Distance-Correlated Shadowing (rng method)\n');
    
    for i = 1:n_users
        Pr = log_normal_shadowing_distance_corr(Pt, Gt, Gr, lambda, L, ...
                                                pathlossExp, std_db, d0, d_km(i));
        PL_total(i) = -10 * log10(Pr + eps);
    end
    
elseif shadowing_model == 3
    %% MODEL 3: Spatially-Correlated Shadowing (Most Realistic)
    fprintf('Using MODEL 3: Spatially-Correlated Shadowing\n');
    
    % Generate spatially-correlated shadowing field
    shadowing_field_db = generate_correlated_shadowing(user_x, user_y, ...
                                                       std_db, correlation_distance);
    
    for i = 1:n_users
        % Calculate base received power (without shadowing)
        Pr0 = friis(Pt, Gt, Gr, lambda, L, d0);
        avg_db = -10.0 * pathlossExp * log10(d_km(i)/d0);
        
        % Apply spatially-correlated shadowing
        powerLoss_db = avg_db + shadowing_field_db(i);
        Pr = Pr0 * 10^(powerLoss_db/10);
        
        PL_total(i) = -10 * log10(Pr + eps);
    end
end

%% 5. Apply Additional Deterministic Obstacle Shadowing (Optional)
shadow_loss = zeros(n_users, 1);
tower_pos = [0, 0];

for i = 1:n_users
    user_pos = [user_x(i), user_y(i)];
    user_dist = d_km(i);
    
    vec_tu = user_pos - tower_pos;
    unit_vec_tu = vec_tu / norm(vec_tu);
    
    current_user_added_loss = 0;
    
    for j = 1:n_obstacles
        obs_center = obstacles(j, 1:2);
        obs_radius = obstacles(j, 3);
        obs_atten = obstacles(j, 4);
        
        vec_to = obs_center - tower_pos;
        projection_dist = dot(vec_to, unit_vec_tu);
        closest_point_on_line = tower_pos + projection_dist * unit_vec_tu;
        dist_to_line = norm(obs_center - closest_point_on_line);
        
        if (dist_to_line < obs_radius) && (projection_dist > 0) && ...
           (projection_dist < (user_dist - obs_radius*0.5))
             current_user_added_loss = max(current_user_added_loss, obs_atten);
        end
    end
    shadow_loss(i) = current_user_added_loss;
end

PL_total = PL_total + shadow_loss;

%% 6. Connectivity Check
is_connected = PL_total <= MAPL_dB;

% Calculate statistics
coverage_percentage = 100 * sum(is_connected) / n_users;
fprintf('\nCoverage Statistics:\n');
fprintf('  Connected Users: %d / %d (%.1f%%)\n', sum(is_connected), n_users, coverage_percentage);
fprintf('  Mean Path Loss: %.2f dB\n', mean(PL_total));
fprintf('  Std Path Loss: %.2f dB\n', std(PL_total));

%% 7. Visualization
figure('Position', [100, 100, 900, 800], 'Color', 'w');
hold on; grid on; axis equal;
box on;

% --- Draw Obstacles (Grey Circles) ---
theta_circ = linspace(0, 2*pi, 100);
for j = 1:n_obstacles
    ox = obstacles(j,1) + obstacles(j,3) * cos(theta_circ);
    oy = obstacles(j,2) + obstacles(j,3) * sin(theta_circ);
    fill(ox, oy, [0.6 0.6 0.6], 'EdgeColor', [0.4 0.4 0.4], ...
         'LineWidth', 2, 'FaceAlpha', 0.7);
    text(obstacles(j,1), obstacles(j,2), ['+' num2str(obstacles(j,4)) 'dB'], ...
        'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
end

% --- Plot Users ---
plot(user_x(~is_connected), user_y(~is_connected), 'r.', ...
     'MarkerSize', 6, 'DisplayName', 'No Service');
plot(user_x(is_connected), user_y(is_connected), 'g.', ...
     'MarkerSize', 6, 'DisplayName', 'Connected');

% --- Plot Tower and Reference Circle ---
plot(0, 0, 'k^', 'MarkerSize', 14, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'Base Station');
circ_x = R_theoretical * cos(theta_circ);
circ_y = R_theoretical * sin(theta_circ);
plot(circ_x, circ_y, 'b--', 'LineWidth', 2, ...
     'DisplayName', 'Theoretical Radius');

% --- Formatting ---
model_names = {'Independent', 'Distance-Correlated', 'Spatially-Correlated'};
title({sprintf('LTE Coverage with %s Shadowing (%.1f%% Coverage)', ...
       model_names{shadowing_model}, coverage_percentage); ...
       sprintf('n=%.1f, \\sigma=%.1f dB', pathlossExp, std_db)}, 'FontSize', 14);
xlabel('Distance East/West (km)', 'FontSize', 12);
ylabel('Distance North/South (km)', 'FontSize', 12);
legend('Location', 'northeastoutside', 'FontSize', 11);
xlim([-simulation_area_km simulation_area_km]);
ylim([-simulation_area_km simulation_area_km]);

set(gca, 'FontSize', 10);

hold off;

%% 8. Path Loss vs Distance Plot
figure('Position', [1050, 100, 700, 600], 'Color', 'w');
hold on; grid on; box on;

[d_sorted, idx] = sort(d_km);
PL_sorted = PL_total(idx);
connected_sorted = is_connected(idx);

plot(d_sorted(connected_sorted), PL_sorted(connected_sorted), 'g.', ...
     'MarkerSize', 6, 'DisplayName', 'Connected');
plot(d_sorted(~connected_sorted), PL_sorted(~connected_sorted), 'r.', ...
     'MarkerSize', 6, 'DisplayName', 'Disconnected');

yline(MAPL_dB, 'k--', 'LineWidth', 2, 'DisplayName', 'MAPL Threshold');

xlabel('Distance from Base Station (km)', 'FontSize', 12);
ylabel('Total Path Loss (dB)', 'FontSize', 12);
title(sprintf('Path Loss vs Distance - %s Shadowing', model_names{shadowing_model}), ...
      'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 10);

hold off;

%% ========== HELPER FUNCTIONS ==========

function Pr = log_normal_shadowing_independent(Pt, Gt, Gr, lambda, L, pathlossExp, std_db, d0, d)
% MODEL 1: Independent shadowing (original method - most realistic)
% Each user gets independent random shadowing

Pr0 = friis(Pt, Gt, Gr, lambda, L, d0);
avg_db = -10.0 * pathlossExp * log10(d/d0);

% Independent random variable for each call
powerLoss_db = avg_db + randn * std_db;

Pr = Pr0 * 10^(powerLoss_db/10);
end

function Pr = log_normal_shadowing_distance_corr(Pt, Gt, Gr, lambda, L, pathlossExp, std_db, d0, d)
% MODEL 2: Distance-correlated shadowing (unrealistic)
% Users at similar distances get similar shadowing

Pr0 = friis(Pt, Gt, Gr, lambda, L, d0);
avg_db = -10.0 * pathlossExp * log10(d/d0);

% Seed based on distance (causes correlation)
rng(round(d*10000), 'twister');
powerLoss_db = avg_db + randn * std_db;

Pr = Pr0 * 10^(powerLoss_db/10);
end

function shadowing_db = generate_correlated_shadowing(x, y, std_db, dcorr)
% MODEL 3: Generate spatially-correlated shadowing field
% Uses exponential correlation model: R(d) = exp(-d/dcorr)
%
% This is the most realistic model used in cellular network planning
% Nearby users have correlated shadowing, distant users are independent

n = length(x);

% Calculate distance matrix between all users
dist_matrix = zeros(n, n);
for i = 1:n
    for j = i:n
        d = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);
        dist_matrix(i,j) = d;
        dist_matrix(j,i) = d;
    end
end

% Create correlation matrix using exponential model
% R(d) = exp(-d/dcorr)
corr_matrix = exp(-dist_matrix / dcorr);

% Generate correlated Gaussian samples using Cholesky decomposition
try
    L = chol(corr_matrix, 'lower');
    % Generate independent samples
    independent_samples = randn(n, 1);
    % Create correlated samples
    correlated_samples = L * independent_samples;
    % Scale to desired standard deviation
    shadowing_db = std_db * correlated_samples;
catch
    % If Cholesky fails (shouldn't happen), fall back to independent
    warning('Cholesky decomposition failed, using independent shadowing');
    shadowing_db = std_db * randn(n, 1);
end

end

function Pr = friis(Pt, Gt, Gr, lambda, L, d)
% Friis free-space propagation model

d_m = d * 1000; % Convert km to meters
Pr = (Pt * Gt * Gr * lambda^2) / ((4 * pi * d_m)^2 * L);

end